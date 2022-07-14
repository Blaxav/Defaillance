#########################################################################################
# Import management
#########################################################################################
if "--install-packages" in ARGS
    import Pkg
    Pkg.add("Random")
    Pkg.add("Plots")
    Pkg.add("JuMP")
    Pkg.add("CPLEX")
    Pkg.add("Gurobi")
    Pkg.add("Dualization")
    Pkg.add("BilevelJuMP")
    Pkg.build("CPLEX")
end

try
    using Random
    using Plots
    using JuMP, CPLEX, BilevelJuMP
    using Printf
    using Gurobi
catch
    println("ERROR : Some recquired packages not found.")
    println("        Run with option --install-packages to install them")
    exit()
end

include("graphGenerator.jl")
include("problemDataGenerator.jl")

struct StochasticProblem
    model::Model
    invest_flow
    invest_prod
    flow
    prod
    unsupplied
    spilled
end

"""
function create_invest_optim_problem
    brief: Creates an stochastic investment optimization problem on a network
        It is possible to invest on every edge of the graph
"""
function create_invest_optim_problem(data)
    #model = Model(CPLEX.Optimizer)
    model = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GRB_ENV), "Threads" => 1, "TimeLimit" => 600))

    # invest variables
    @variable(model, 0 <= invest_flow[e in data.network.edges])
    @variable(model, 0 <= invest_prod[n in 1:data.network.N; 
        data.has_production[n] == 1])

    # Flow variables
    @variable(model, flow[s in 1:data.S, e in data.network.edges, t in 1:data.T])

    # Production variables
    @variable(model, 0 <= prod[s in 1:data.S, n in 1:data.network.N, t in 1:data.T; 
        data.scenario[s].has_production[n] == 1])
    
    # Production constraints
    @constraint(model, prod_max[s in 1:data.S, n in 1:data.network.N, t in 1:data.T; 
        data.has_production[n] == 1],
        prod[s,n,t] <= invest_prod[n])
    
    @constraint(model, grad_positive[s in 1:data.S, n in 1:data.network.N, t in 1:data.T; 
        data.has_production[n] == 1],
        prod[s,n,t] <= (t > 1 ? prod[s,n,t-1] : prod[s,n,data.T]) + data.scenario[s].grad_prod*invest_prod[n] )
    @constraint(model, grad_negative[s in 1:data.S, n in 1:data.network.N, t in 1:data.T; 
        data.has_production[n] == 1],
        prod[s,n,t] >= (t > 1 ? prod[s,n,t-1] : prod[s,n,data.T]) - data.scenario[s].grad_prod*invest_prod[n] )

    # Loss of load variables
    @variable(model, 0 <= unsupplied[s in 1:data.S, i in 1:data.network.N, t in 1:data.T])
    @variable(model, 0 <= spilled[s in 1:data.S, i in 1:data.network.N, t in 1:data.T])

    # Flow bounds
    @constraint(model, flow_max_positive[s in 1:data.S, e in data.network.edges, t in 1:data.T], 
        flow[s,e,t] <= data.scenario[s].flow_init[e] + invest_flow[e])
    @constraint(model, flow_max_negative[s in 1:data.S, e in data.network.edges, t in 1:data.T], 
        -(invest_flow[e] + data.scenario[s].flow_init[e]) <= flow[s,e,t])

    
    @constraint(model, flow_conservation[s in 1:data.S, n in 1:data.network.N, t in 1:data.T], 
        sum(flow[s,e,t] for e in data.network.edges if e.to == n) - 
        sum(flow[s,e,t] for e in data.network.edges if e.from == n) + 
        (data.scenario[s].has_production[n] == 1 ? prod[s,n,t] : 0) + 
        unsupplied[s,n,t] - spilled[s,n,t] == data.scenario[s].demands[n,t]
        )

    
    # Constraint to lead heuristic
    @constraint(model, invest_cost,
        sum( [data.invest_flow_cost[e] * invest_flow[e] for e in data.network.edges]) +
        sum( data.invest_prod_cost[n] * invest_prod[n] for n in 1:data.network.N 
            if data.has_production[n] == 1) >= 0.0
        )
    
    @objective(model, Min,
        sum( data.invest_flow_cost[e] * invest_flow[e] for e in data.network.edges ) +
        sum( data.invest_prod_cost[n] * invest_prod[n] for n in 1:data.network.N 
            if data.has_production[n] == 1) +
        # Sum on scenarios
        sum( data.probability[s] *
            (   
            # Sum on time steps
            sum( 
                (
                # unsupplied costs
                sum( data.scenario[s].unsupplied_cost * unsupplied[s,n,t] 
                    for n in 1:data.network.N ) +
                # production costs
                sum( data.scenario[s].prod_cost[n][t] * prod[s,n,t] 
                    for n in 1:data.network.N 
                    if data.scenario[s].has_production[n] == 1) +
                # flow cost
                sum( data.scenario[s].flow_cost[e] * flow[s,e,t] 
                    for e in data.network.edges )
                ) for t in 1:data.T
            )
            ) for s in 1:data.S
        )
    )

    return StochasticProblem(model, invest_flow, invest_prod, flow, prod, unsupplied, spilled)
end


"""
function counting_unsupplied_scenario
    brief: Computes the number of nodes which loss of load for a given scenario after optimization
"""
function counting_unsupplied_scenario(stoch_prob, scenario, epsilon, data)
    return sum( value(stoch_prob.unsupplied[scenario,n,t]) > epsilon ? 1 : 0 for n in 1:data.network.N, t in 1:data.T )
end



"""
function add_unsupplied_counter_constraint
    brief: Adds a binary variable on each node to say if it has unsupplied energy 
        and a constraints to limit the number of node with unsupplied energy to n_unsupplied
    args:
        stoch_prob : an instance of StochasticProblem
        max_unsupplied: maximum number of nodes with unsupplied enery in an optimal solution
        data : an instance of StochasticProblemData
    returns: references to new variables and new constraint
"""
function add_unsupplied_counter_constraint(stoch_prob, max_unsupplied, epsilon_cnt, data)
   
    @variable(stoch_prob.model, has_unsupplied[s in 1:data.S, i in 1:data.network.N, t in 1:data.T], Bin)

    # Variables behaviour
    @constraint(stoch_prob.model, unsupplied_to_zero[s in 1:data.S, n in 1:data.network.N, t in 1:data.T],
        has_unsupplied[s,n,t] <= (1/epsilon_cnt)*stoch_prob.unsupplied[s,n,t] )
    @constraint(stoch_prob.model, unsupplied_to_one[s in 1:data.S, n in 1:data.network.N, t in 1:data.T],
        2*data.scenario[s].demands[n,t]*has_unsupplied[s,n,t] >= stoch_prob.unsupplied[s,n,t] - epsilon_cnt )
    
    @constraint(stoch_prob.model, unsupplied_cnt, sum(data.probability .* sum( sum(has_unsupplied[:,n,t] for n = 1:data.network.N) for t in 1:data.T )) <= max_unsupplied )
    return has_unsupplied, unsupplied_cnt
end


"""
function add_unsupplied_counter_constraint
    brief: Adds a binary variable on each node to say if it has unsupplied energy
    args:
        stoch_prob : an instance of StochasticProblem
        max_unsupplied: maximum number of nodes with unsupplied enery in an optimal solution
        data : an instance of StochasticProblemData
    returns: references to new variables and new constraint
"""
function add_unsupplied_counter_variables(stoch_prob, max_unsupplied, epsilon_cnt, data)
   
    @variable(stoch_prob.model, has_unsupplied[s in 1:data.S, i in 1:data.network.N, t in 1:data.T], Bin)

    # Variables behaviour
    @constraint(stoch_prob.model, unsupplied_to_zero[s in 1:data.S, n in 1:data.network.N, t in 1:data.T],
        has_unsupplied[s,n,t] <= (1/epsilon_cnt)*stoch_prob.unsupplied[s,n,t] )
    @constraint(stoch_prob.model, unsupplied_to_one[s in 1:data.S, n in 1:data.network.N, t in 1:data.T],
        2*data.scenario[s].demands[n,t]*has_unsupplied[s,n,t] >= stoch_prob.unsupplied[s,n,t] - epsilon_cnt )
    
    return has_unsupplied
end


"""
function solve
    brief: solves a stochastic problem directly
    returns: time to solve
"""
function solve(stoch_prob, silent_mode)

    silent_mode == true ? set_silent(stoch_prob.model) : nothing
    timer = @elapsed optimize!(stoch_prob.model)

    # Always unset silent before leaving
    unset_silent(stoch_prob.model)

    return timer
end


"""
function investment_cost
    brief: computes the investment cost of the inner solution of stoch_prob
"""
function investment_cost(stoch_prob, data)
    sum( data.invest_flow_cost[e] * value(stoch_prob.invest_flow[e]) for e in data.network.edges) +
    sum( data.invest_prod_cost[n] * value(stoch_prob.invest_prod[n]) for n in 1:data.network.N 
            if data.has_production[n] == 1)
end


"""
function investment_heuristic
    brief : heuristic to force investment in order to set unsatisfied demands under
        a given value
"""
function investment_heuristic(stoch_prob, data, max_unsupplied, relative_gap, silent_mode, print_log)

    print_log == true ? @printf("%-20s%-20s%-20s%-20s%-20s%-20s%-15s%-20s%-20s\n", "Invest_min", "Invest_max", "Gap", "Constraint value", 
            "Invest cost", "Objective", "Unsupplied", "Optim time", "Counting time" ) : nothing

    # Initialization
    t_optim = @elapsed solve(stoch_prob, silent_mode)
    t_counting = @elapsed unsupplied_cnt = [counting_unsupplied_scenario(stoch_prob, s, 0.0, data) for s in 1:data.S]

    invest_min = 0.0
    invest_max = 100*objective_value(stoch_prob.model)

    #println(invest_min, "  ", invest_max)
    global best_obj = 0.0

    # Setting initial invest min and max
    if sum( data.probability .* unsupplied_cnt ) > max_unsupplied
        global invest_min = investment_cost(stoch_prob, data)
    else
        global invest_max = investment_cost(stoch_prob, data)
    end
    alpha = 0.5
    rhs = 0

    global best_unsup_val = 0.0
    #println(invest_min, "  ", invest_max)

    print_log == true ? @printf("%-20.6e%-20.6e%-20.2e%-20.6e%-20.6e%-20.6e%-15.2f%-20.6e%-20.6e\n", invest_min, invest_max, 
            (invest_max - invest_min)/invest_max, rhs, investment_cost(stoch_prob, data),
            objective_value(stoch_prob.model), sum( data.probability .* unsupplied_cnt ), 
            t_optim, t_counting ) : nothing

    while invest_max - invest_min > max(relative_gap*invest_max,1e-6)

        global rhs = (1-alpha)*invest_max + alpha*invest_min

        set_normalized_rhs(constraint_by_name(stoch_prob.model, "invest_cost"), rhs)
    
        global t_optim = @elapsed solve(stoch_prob, silent_mode)
        global t_counting = @elapsed global unsupplied_cnt = [ counting_unsupplied_scenario(stoch_prob, s, 0.0, data) for s in 1:data.S ]
        

        if sum( (data.probability .* unsupplied_cnt) ) > max_unsupplied
            global invest_min = rhs
        else
            global invest_max = rhs
            global best_obj = objective_value(stoch_prob.model)
            global best_unsup_val = sum( (data.probability .* unsupplied_cnt) )
        end

        print_log == true ? @printf("%-20.6e%-20.6e%-20.2e%-20.6e%-20.6e%-20.6e%-15.2f%-20.6e%-20.6e\n", invest_min, invest_max, 
            (invest_max - invest_min)/invest_max, rhs, investment_cost(stoch_prob, data),
            objective_value(stoch_prob.model), sum( data.probability .* unsupplied_cnt ), 
            t_optim, t_counting ) : nothing
    
    end

    #println("Unsup final : ", best_unsup_val)
    set_normalized_rhs(constraint_by_name(stoch_prob.model, "invest_cost"), 0)
    return invest_max, best_obj, best_unsup_val
end