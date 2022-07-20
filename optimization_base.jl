#########################################################################################
#= Import management
    Graph of includes
    
    import_management.jl --> graph.jl --> problemDataGenerator.jl --> optimization_base.jl  

=########################################################################################
include("problemDataGenerator.jl")

# Definition of global environment for Gurobi
# Used as unique environment for every JuMP model
#const GRB_ENV = Gurobi.Env()

############################################################################
# Creating variables
############################################################################
function variables_investment_flow(model, dataif, n_scenarios = 1)
    @variable(model, 0 <= invest_flow[e in data.network.edges])
end

function variables_investment_production(model, data)
    @variable(model, 0 <= invest_prod[n in 1:data.network.N; 
        data.has_production[n] == 1])
end

function variables_flow(model, data)
    @variable(model, flow[e in data.network.edges, t in 1:data.T])
end

function variables_absolute_value_flow(model, data)
    @variable(model, 0 <= flow_abs[e in data.network.edges, t in 1:data.T])
end 
    
function variables_production(model, data)
    @variable(model, 0 <= prod[n in 1:data.network.N, t in 1:data.T; 
        data.has_production[n] == 1])
end

function variables_unsupplied_energy(model, data)
    @variable(model, 0 <= unsupplied[i in 1:data.network.N, t in 1:data.T])
end

function variables_spilled_energy(model, data)    
    @variable(model, 0 <= spilled[i in 1:data.network.N, t in 1:data.T])
end

function variables_epigraph_on_scenario(model, data)
    @variable(model, 0 <= theta[s in 1:data.S])
end

function variables_epigraph_total(model, data)
    @variable(model, theta_sum)
end



############################################################################
# Creating constraints
# Requires reference to model and each used variable
############################################################################
function constraint_max_prod(model, prod, invest_prod, data, n_scenarios=1)
    if n_scenarios == 1
        @constraint(model, prod_max[n in 1:data.network.N, t in 1:data.T; 
            data.has_production[n] == 1],
            prod[n,t] <= invest_prod[n])
    end
end
    
function constraint_production_positive_gradient(model, invest_prod, prod, data, n_scenarios=1)
    if n_scenarios == 1
        @constraint(model, grad_positive[n in 1:data.network.N, t in 1:data.T; 
            data.has_production[n] == 1],
            prod[n,t] <= (t > 1 ? prod[n,t-1] : prod[n,data.T]) + data.grad_prod*invest_prod[n] )
    end
end

function constraint_production_negative_gradient(model, invest_prod, prod, data, n_scenarios = 1)
    @constraint(model, grad_negative[n in 1:data.network.N, t in 1:data.T; 
        data.has_production[n] == 1],
        prod[n,t] >= (t > 1 ? prod[n,t-1] : prod[n,data.T]) - data.grad_prod*invest_prod[n] )
end

function constraint_flow_max_positive(model, flow, invest_flow, data, n_scenarios = 1)
    @constraint(model, flow_max_positive[e in data.network.edges, t in 1:data.T], 
        flow[e,t] <= data.flow_init[e] + invest_flow[e])
end

function constraint_flow_max_negative(model, flow, invest_flow, data, n_scenarios = 1)
    @constraint(model, flow_max_negative[e in data.network.edges, t in 1:data.T], 
        -(invest_flow[e] + data.flow_init[e]) <= flow[e,t])
end

function constraint_flow_max(model, flow, invest_flow, data, n_scenarios = 1)
    constraint_flow_max_positive(model, flow, invest_flow, data, n_scenarios)
    constraint_flow_max_negative(model, flow, invest_flow, data, n_scenarios)
end

function constraint_absolute_flow_behavior(model, flow, flow_abs, n_scenarios = 1)
    # Absolute value of flow in cost
    @constraint(model, flow_abs_positive[e in data.network.edges, t in 1:data.T], 
        flow_abs[e,t] >= flow[e,t])
    @constraint(model, flow_abs_negative[e in data.network.edges, t in 1:data.T], 
        flow_abs[e,t] >= -flow[e,t])
end

function constraint_flow_conservation(model, prod, flow, unsupplied, spilled, data; n_scenarios = 1, which_scenario = -1)
    if n_scenarios == 1
        @constraint(model, flow_conservation[n in 1:data.network.N, t in 1:data.T], 
            sum(flow[e,t] for e in data.network.edges if e.to == n) - 
            sum(flow[e,t] for e in data.network.edges if e.from == n) + 
            (data.has_production[n] == 1 ? prod[n,t] : 0) + 
            unsupplied[n,t] - spilled[n,t] == data.scenario[which_scenario].demands[n,t]
            )
    end
end


function constraint_minimum_investment(model, invest_flow, invest_prod, data)
    @constraint(model, invest_cost,
        sum( [data.invest_flow_cost[e] * invest_flow[e] for e in data.network.edges]) +
        sum( data.invest_prod_cost[n] * invest_prod[n] for n in 1:data.network.N 
            if data.has_production[n] == 1) >= 0.0
        )
end
    
function constraint_total_epigraph(model, theta_sum, theta, data)
    @constraint(model, epigrah_cost, theta_sum == sum( data.probability[s] * theta[s] for s in 1:data.S) )
end

############################################################################
# Objectives
############################################################################
function objective_subproblem(model, unsupplied, prod, flow_abs, data; which_scenario=1)
    @objective(model, Min,
        # Sum on time steps
        sum(
        (
            # unsupplied costs
            sum( data.scenario[which_scenario].unsupplied_cost * unsupplied[n,t] 
                for n in 1:data.network.N ) +
            # production costs
            sum( data.scenario[which_scenario].prod_cost[n][t] * prod[n,t] 
                for n in 1:data.network.N 
                if data.has_production[n] == 1) +
            # flow cost
            sum( data.scenario[which_scenario].flow_cost[e] * flow_abs[e,t] 
                for e in data.network.edges )
        ) 
        for t in 1:data.T )
    )
end

function objective_master_problem(model, invest_flow, invest_prod, theta_sum, data)
    @objective(model, Min,
        sum( data.invest_flow_cost[e] * invest_flow[e] for e in data.network.edges ) +
        sum( data.invest_prod_cost[n] * invest_prod[n] for n in 1:data.network.N 
            if data.has_production[n] == 1) +
        theta_sum
    )
end

function set_invest_free_master_problem(master, data, perturbation)
    @objective(master.model, Min,
        sum( perturbation * master.invest_flow[e] for e in data.network.edges ) +
        sum( perturbation * master.invest_prod[n] for n in 1:data.network.N 
            if data.has_production[n] == 1) +
        master.theta_sum
    )
end

function set_initial_objective_master_problem(master, data)
    objective_master_problem(
        master.model, master.invest_flow, 
        master.invest_prod, master.theta_sum, 
        data
    )
end

############################################################################
# Solution methods and solutions getters
############################################################################
function solve(problem; silent_mode=false)
    silent_mode == true ? set_silent(problem.model) : nothing
    timer = @elapsed optimize!(problem.model)

    # Always unset silent before leaving
    unset_silent(problem.model)

    return timer
end

function get_objective_value(problem)
    objective_value(problem.model)
end


function print_solution(problem, data; null_tolerance=1e-6)
    invest_prod_sol = value.(problem.invest_prod)
    invest_flow_sol = value.(problem.invest_flow)
    for n in production_nodes(data)
        if invest_prod_sol[n] > null_tolerance
            @printf("Node %-15i%-10.3f\n", n, invest_prod_sol[n])
        end
    end
    for e in data.network.edges
        if invest_flow_sol[e] > null_tolerance
            @printf("%-20s%-10.2f\n", e, invest_flow_sol[e])
        end
    end
end




############################################################################
# Test
############################################################################
if PROGRAM_FILE == "graphGenerator.jl"
    N = 3
    graph_density = 20
    seed = 3
    time_graph = @elapsed network = create_network(N, graph_density, seed, plotGraph = false, drawGraph = true)


    seed >= 0 ? Random.seed!(seed) : nothing

    scenarios = 1
    time_steps = 3
    demand_range = 100:500
    prod_cost_range = 30:80
    unsupplied_cost = 1000
    epsilon_flow = 1.0
    grad_prod = 0.2
    invest_cost_range = 500:1000
    invest_prod_range = 500:1000
    flow_init_max = 5

    time_data = @elapsed data = investment_problem_data_generator(scenarios, network, time_steps, demand_range, 
    prod_cost_range, unsupplied_cost, epsilon_flow, flow_init_max, grad_prod, invest_cost_range, invest_prod_range)

    model = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GRB_ENV), "Threads" => 1, "TimeLimit" => 600))
    invest_flow = variables_investment_flow(model, data)
    invest_prod = variables_investment_production(model, data)

    flow = variables_flow(model, data)
    prod = variables_production(model, data)

    flow_abs = variables_absolute_value_flow(model, data)

    unsupplied = variables_unsupplied_energy(model, data)
    spilled = variables_spilled_energy(model, data)

    constraint_max_prod(model, prod, invest_prod, data)

    constraint_production_positive_gradient(model, prod, data)
    constraint_production_negative_gradient(model, prod, data)

    constraint_flow_max_positive(model, flow, invest_flow, data)
    constraint_flow_max_negative(model, flow, invest_flow, data)

    constraint_absolute_flow_behavior(model, flow, flow_abs)
    constraint_flow_conservation(model, prod, flow, unsupplied, spilled, data;  which_scenario = 1)

    objective_subproblem(model, unsupplied, prod, flow_abs, data, 1)
    print(model)

    optimize!(model)
end