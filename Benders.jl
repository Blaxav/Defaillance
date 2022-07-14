#########################################################################################
# Import management
#########################################################################################
if "--install-packages" in ARGS
    import Pkg
    Pkg.add("Random")
    Pkg.add("Plots")
    Pkg.add("JuMP")
    Pkg.add("CPLEX")
    Pkg.add("Dualization")
    Pkg.add("BilevelJuMP")
    Pkg.build("CPLEX")
end

try
    using Random
    using Plots
    using JuMP, CPLEX, BilevelJuMP
    using Printf
catch
    println("ERROR : Some recquired packages not found.")
    println("        Run with option --install-packages to install them")
    exit()
end

include("graphGenerator.jl")
include("problemDataGenerator.jl")
include("bilevelProblemGenerator.jl")
include("stochasticProblemGenerator.jl")
include("dataFromAntaresFormat.jl")


const GRB_ENV = Gurobi.Env()

#########################################################################################
# Benders problem creation
#########################################################################################
struct BendersMasterProblem
    model::Model
    invest_flow
    invest_prod
    theta_sum # expectation of epigraph variables
    theta # epigraph variable for each scenario
end


struct BendersSubroblem
    model::Model
    invest_flow
    invest_prod
    flow
    prod
    unsupplied
    spilled
end


struct BendersCountingSubroblem
    model::Model
    invest_flow
    invest_prod
    flow
    prod
    unsupplied
    spilled
    cost_constraint
    has_unsupplied
end


function create_master_benders_problem(data)
    model = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GRB_ENV), "Threads" => 1, "TimeLimit" => 600))
    #model = Model(CPLEX.Optimizer)


    # invest variables
    @variable(model, 0 <= invest_flow[e in data.network.edges])
    @variable(model, 0 <= invest_prod[n in 1:data.network.N; 
        data.has_production[n] == 1])
    
    # Constraint to lead heuristic
    @constraint(model, invest_cost,
        sum( [data.invest_flow_cost[e] * invest_flow[e] for e in data.network.edges]) +
        sum( data.invest_prod_cost[n] * invest_prod[n] for n in 1:data.network.N 
            if data.has_production[n] == 1) >= 0.0
        )

    # Epigraph variables
    @variable(model, theta_sum)
    @variable(model, 0 <= theta[s in 1:data.S])

    @constraint(model, epigrah_cost, theta_sum == sum( data.probability[s] * theta[s] for s in 1:data.S) )
    
    @objective(model, Min,
        sum( data.invest_flow_cost[e] * invest_flow[e] for e in data.network.edges ) +
        sum( data.invest_prod_cost[n] * invest_prod[n] for n in 1:data.network.N 
            if data.has_production[n] == 1) +
        theta_sum
    )

    return BendersMasterProblem(model, invest_flow, invest_prod, theta_sum, theta)
end


function create_benders_subproblem(data, s)
    model = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GRB_ENV), "Threads" => 1, "TimeLimit" => 600))
    #model = Model(CPLEX.Optimizer)

    # invest variables
    @variable(model, 0 <= invest_flow[e in data.network.edges])
    @variable(model, 0 <= invest_prod[n in 1:data.network.N; 
        data.has_production[n] == 1])

    # Flow variables
    @variable(model, flow[e in data.network.edges, t in 1:data.T])

    # Production variables
    @variable(model, 0 <= prod[n in 1:data.network.N, t in 1:data.T; 
        data.scenario[s].has_production[n] == 1])
    
    # Production constraints
    @constraint(model, prod_max[n in 1:data.network.N, t in 1:data.T; 
        data.has_production[n] == 1],
        prod[n,t] <= invest_prod[n])
    
    @constraint(model, grad_positive[n in 1:data.network.N, t in 1:data.T; 
        data.has_production[n] == 1],
        prod[n,t] <= (t > 1 ? prod[n,t-1] : prod[n,data.T]) + data.scenario[s].grad_prod*invest_prod[n] )
    @constraint(model, grad_negative[n in 1:data.network.N, t in 1:data.T; 
        data.has_production[n] == 1],
        prod[n,t] >= (t > 1 ? prod[n,t-1] : prod[n,data.T]) - data.scenario[s].grad_prod*invest_prod[n] )

    # Loss of load variables
    @variable(model, 0 <= unsupplied[i in 1:data.network.N, t in 1:data.T])
    @variable(model, 0 <= spilled[i in 1:data.network.N, t in 1:data.T])

    # Flow bounds
    @constraint(model, flow_max_positive[e in data.network.edges, t in 1:data.T], 
        flow[e,t] <= data.scenario[s].flow_init[e] + invest_flow[e])
    @constraint(model, flow_max_negative[e in data.network.edges, t in 1:data.T], 
        -(invest_flow[e] + data.scenario[s].flow_init[e]) <= flow[e,t])

    
    @constraint(model, flow_conservation[n in 1:data.network.N, t in 1:data.T], 
        sum(flow[e,t] for e in data.network.edges if e.to == n) - 
        sum(flow[e,t] for e in data.network.edges if e.from == n) + 
        (data.scenario[s].has_production[n] == 1 ? prod[n,t] : 0) + 
        unsupplied[n,t] - spilled[n,t] == data.scenario[s].demands[n,t]
        )

    
    @objective(model, Min,
        # Sum on time steps
        sum(
            (
            # unsupplied costs
            sum( data.scenario[s].unsupplied_cost * unsupplied[n,t] 
                for n in 1:data.network.N ) +
            # production costs
            sum( data.scenario[s].prod_cost[n][t] * prod[n,t] 
                for n in 1:data.network.N 
                if data.scenario[s].has_production[n] == 1) +
            # flow cost
            sum( data.scenario[s].flow_cost[e] * flow[e,t] 
                for e in data.network.edges )
            ) for t in 1:data.T
        )
    )

    return BendersSubroblem(model, invest_flow, invest_prod, flow, prod, unsupplied, spilled)
end

function counting_unsupplied_subproblem(stoch_prob, epsilon, data)
    return sum( value(stoch_prob.unsupplied[n,t]) > epsilon ? 1 : 0 for n in 1:data.network.N, t in 1:data.T )
end


function create_master_benders_problem(data)
    model = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GRB_ENV), "Threads" => 1, "TimeLimit" => 600))
    #model = Model(CPLEX.Optimizer)

    # invest variables
    @variable(model, 0 <= invest_flow[e in data.network.edges])
    @variable(model, 0 <= invest_prod[n in 1:data.network.N; 
        data.has_production[n] == 1])
    
    # Constraint to lead heuristic
    @constraint(model, invest_cost,
        sum( [data.invest_flow_cost[e] * invest_flow[e] for e in data.network.edges]) +
        sum( data.invest_prod_cost[n] * invest_prod[n] for n in 1:data.network.N 
            if data.has_production[n] == 1) >= 0.0
        )

    # Epigraph variables
    @variable(model, theta_sum)
    @variable(model, 0 <= theta[s in 1:data.S])

    @constraint(model, epigrah_cost, theta_sum == sum( data.probability[s] * theta[s] for s in 1:data.S) )
    
    @objective(model, Min,
        sum( data.invest_flow_cost[e] * invest_flow[e] for e in data.network.edges ) +
        sum( data.invest_prod_cost[n] * invest_prod[n] for n in 1:data.network.N 
            if data.has_production[n] == 1) +
        theta_sum
    )

    return BendersMasterProblem(model, invest_flow, invest_prod, theta_sum, theta)
end


function create_benders_subproblem_with_counting_unsupplied(data, s, epsilon_cnt)
    model = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GRB_ENV), "Threads" => 1, "TimeLimit" => 600))
    #model = Model(CPLEX.Optimizer)

    # invest variables
    @variable(model, 0 <= invest_flow[e in data.network.edges])
    @variable(model, 0 <= invest_prod[n in 1:data.network.N; 
        data.has_production[n] == 1])

    # Flow variables
    @variable(model, flow[e in data.network.edges, t in 1:data.T])

    # Production variables
    @variable(model, 0 <= prod[n in 1:data.network.N, t in 1:data.T; 
        data.scenario[s].has_production[n] == 1])
    
    # Production constraints
    @constraint(model, prod_max[n in 1:data.network.N, t in 1:data.T; 
        data.has_production[n] == 1],
        prod[n,t] <= invest_prod[n])
    
    @constraint(model, grad_positive[n in 1:data.network.N, t in 1:data.T; 
        data.has_production[n] == 1],
        prod[n,t] <= (t > 1 ? prod[n,t-1] : prod[n,data.T]) + data.scenario[s].grad_prod*invest_prod[n] )
    @constraint(model, grad_negative[n in 1:data.network.N, t in 1:data.T; 
        data.has_production[n] == 1],
        prod[n,t] >= (t > 1 ? prod[n,t-1] : prod[n,data.T]) - data.scenario[s].grad_prod*invest_prod[n] )

    # Loss of load variables
    @variable(model, 0 <= unsupplied[i in 1:data.network.N, t in 1:data.T])
    @variable(model, 0 <= spilled[i in 1:data.network.N, t in 1:data.T])

    # Variable saying if a node has unsupplied energy
    @variable(model, has_unsupplied[i in 1:data.network.N, t in 1:data.T], Bin)

    # Binary Variables behaviour
    @constraint(model, unsupplied_to_zero[n in 1:data.network.N, t in 1:data.T],
        has_unsupplied[n,t] <= (1/epsilon_cnt)*unsupplied[n,t] )
    @constraint(model, unsupplied_to_one[n in 1:data.network.N, t in 1:data.T],
        2*data.scenario[s].demands[n,t]*has_unsupplied[n,t] >= unsupplied[n,t] - epsilon_cnt )

    # Flow bounds
    @constraint(model, flow_max_positive[e in data.network.edges, t in 1:data.T], 
        flow[e,t] <= data.scenario[s].flow_init[e] + invest_flow[e])
    @constraint(model, flow_max_negative[e in data.network.edges, t in 1:data.T], 
        -(invest_flow[e] + data.scenario[s].flow_init[e]) <= flow[e,t])

    
    @constraint(model, flow_conservation[n in 1:data.network.N, t in 1:data.T], 
        sum(flow[e,t] for e in data.network.edges if e.to == n) - 
        sum(flow[e,t] for e in data.network.edges if e.from == n) + 
        (data.scenario[s].has_production[n] == 1 ? prod[n,t] : 0) + 
        unsupplied[n,t] - spilled[n,t] == data.scenario[s].demands[n,t]
        )


    # Objective cost constraint
    @constraint(model, cost_constraint, 
        sum(
            (
            # unsupplied costs
            sum( data.scenario[s].unsupplied_cost * unsupplied[n,t] 
                for n in 1:data.network.N ) +
            # production costs
            sum( data.scenario[s].prod_cost[n][t] * prod[n,t] 
                for n in 1:data.network.N 
                if data.scenario[s].has_production[n] == 1) +
            # flow cost
            sum( data.scenario[s].flow_cost[e] * flow[e,t] 
                for e in data.network.edges )
            ) for t in 1:data.T
        )
        <= 0
    )
    
    @objective(model, Min,
        # Sum on time steps
        sum(# Number of unsupplied nodes
            ( 
                sum(has_unsupplied[n,t] for n in 1:data.network.N )
            ) for t in 1:data.T
        )
    )

    return BendersCountingSubroblem(model, invest_flow, invest_prod, flow, prod, unsupplied, spilled, cost_constraint, has_unsupplied)
end




#########################################################################################
# Benders resolution
#########################################################################################
function get_master_value(bendersMaster::BendersMasterProblem)
    return objective_value(bendersMaster.model)
end

function get_master_solution(bendersMaster::BendersMasterProblem)
    candidate_flow = value.(bendersMaster.invest_flow)
    candidate_prod = value.(bendersMaster.invest_prod)
    return candidate_flow, candidate_prod
end

function fix_first_stage_candidate(subproblems, candidate_flow, candidate_prod, data)
    for s in 1:data.S
        # Can not broadcast fix on edges as they are modeled by a dictionary
        # Seems that broadcasting needs flat data
        for e in data.network.edges
            fix(subproblems[s].invest_flow[e], candidate_flow[e]; force=true)   
        end
        for n in 1:network.N 
            if data.scenario[1].has_production[n] == 1
                fix.(subproblems[s].invest_prod[n], candidate_prod[n]; force=true)
            end
        end
        #fix.(subproblems[s].invest_prod, candidate_prod; force=true)
    end
end

function investment_cost(invest_flow, invest_prod, data)
    sum( data.invest_flow_cost[e] * invest_flow[e] for e in data.network.edges) +
    sum( data.invest_prod_cost[n] * invest_prod[n] for n in 1:data.network.N 
            if data.has_production[n] == 1)
end


function master_solution_investment_cost(master, data)
    solution_invest_prod = value.(master.invest_prod)
    solution_invest_flow = value.(master.invest_flow)

    #println(solution_invest_prod)
    #println(solution_invest_flow)

    investment_cost(solution_invest_flow, solution_invest_prod, data)
end


function benders_sequential(master, subproblems, data, n_iteration_log, algo)

    global best_UB = 1e12
    LB = -1e6
    global best_invest_flow = Dict(zip(
        data.network.edges,
        [0.0 for e in data.network.edges]
        ))
    global separation_flow = Dict(zip(
            data.network.edges,
            [0.0 for e in data.network.edges]
            ))
    global best_invest_prod = Dict(zip(
        [i for i in 1:network.N if data.scenario[1].has_production[i] == 1],
        [0.0 for i in 1:sum(data.scenario[1].has_production)]
        ))
    global separation_prod = Dict(zip(
            [i for i in 1:network.N if data.scenario[1].has_production[i] == 1],
            [0.0 for i in 1:sum(data.scenario[1].has_production)]
            ))
    global iteration = 0
    global alpha = 0.5


    global grad_flow = Dict(zip(
            data.network.edges,
            [0.0 for e in data.network.edges]
            ))
    global grad_prod = Dict(zip(
        [i for i in 1:network.N if data.scenario[1].has_production[i] == 1],
        [0.0 for i in 1:sum(data.scenario[1].has_production)]
        ))

    global log_freq = n_iteration_log
    global stop = false

    #while best_UB - LB > 1e-6*best_UB
    while stop == false

        iteration += 1
        
        # 1. Solving master problem
        solve(master, true)
        
        # 2. Get master solution (investment candidates)
        candidate_flow, candidate_prod = get_master_solution(master)
        LB = objective_value(master.model)

        for e in data.network.edges
            global separation_flow[e] = alpha * candidate_flow[e] + (1-alpha) * best_invest_flow[e]
        end
        for n in 1:network.N 
            if data.scenario[1].has_production[n] == 1
                global separation_prod[n] = alpha * candidate_prod[n] + (1-alpha) * best_invest_prod[n]
            end
        end

        # 3. Fix investment candidates in subproblems
        fix_first_stage_candidate(subproblems, separation_flow, separation_prod, data)
    
        
        # 4. Solve Subproblems and get subprgradients
        local UB = investment_cost(separation_flow, separation_prod, data)
        for s in 1:data.S
            solve(subproblems[s], true)
            UB += data.probability[s] * objective_value(subproblems[s].model)
        end

        # Update UB
        if UB < best_UB
            global alpha = min(1.0, 1.2*alpha)
            global best_UB = UB
            global best_invest_prod = separation_prod
            global best_invest_flow = separation_flow
        else
            global alpha = max(0.1, 0.8*alpha)
        end


        # Log
        if iteration % log_freq == 0 #|| best_UB - LB <= 1e-6*best_UB
            @printf("%-10i%-20.6e%-20.6e%-15.3e%15.3f\n", iteration, LB, best_UB, (best_UB-LB)/best_UB, alpha)
        end

        
        # Stopping criterion
        if best_UB - LB <= 1e-6*best_UB
            global stop = true
        end


        # Get subgradients and build cut
        if stop == false
            local rhs = 0.0

            if algo == "monocut"
                for s in 1:data.S   
                    rhs += data.probability[s] * objective_value(subproblems[s].model)
    
                    if s == 1
                        global grad_flow = data.probability[s] * reduced_cost.(subproblems[s].invest_flow)
                        global grad_prod = data.probability[s] * reduced_cost.(subproblems[s].invest_prod)
                    else
                        for e in data.network.edges
                            global grad_flow[e] += data.probability[s] * reduced_cost(subproblems[s].invest_flow[e])
                        end
                        for n in 1:data.network.N
                            if data.scenario[s].has_production[n] == 1
                                global grad_prod[n] += data.probability[s] * reduced_cost(subproblems[s].invest_prod[n])
                            end
                        end
                    end
                end
            
                @constraint(master.model, master.theta_sum >= 
                rhs +
                sum( grad_flow[e] * (master.invest_flow[e] - separation_flow[e]) for e in data.network.edges) +
                sum( grad_prod[n] * (master.invest_prod[n] - separation_prod[n]) for n in 1:data.network.N if data.scenario[1].has_production[n] == 1)
                )
            elseif algo == "multicut"
                for s in 1:data.S   
                    rhs = objective_value(subproblems[s].model)
    
                    global grad_flow = reduced_cost.(subproblems[s].invest_flow)
                    global grad_prod = reduced_cost.(subproblems[s].invest_prod)

                    @constraint(master.model, master.theta[s] >= 
                        rhs +
                        sum( grad_flow[e] * (master.invest_flow[e] - separation_flow[e]) for e in data.network.edges) +
                        sum( grad_prod[n] * (master.invest_prod[n] - separation_prod[n]) for n in 1:data.network.N if data.scenario[1].has_production[n] == 1)
                        )
                end
            else
                println("Unknown algorithm ", algo)
                exit()
            end

            
        end
    end
end

"""
function investment_heuristic
    brief : heuristic to force investment in order to set unsatisfied demands under
        a given value
"""
function investment_heuristic_benders(benders_master, subproblems, data, max_unsupplied, relative_gap, silent_mode, print_log)

    print_log == true ? @printf("%-20s%-20s%-20s%-20s%-20s%-20s%-15s%-20s%-20s\n", "Invest_min", "Invest_max", "Gap", "Constraint value", 
            "Invest cost", "Objective", "Unsupplied", "Optim time", "Counting time" ) : nothing

    # Initialization
    t_optim = @elapsed benders_sequential(benders_master, subproblems, data, 50)

    t_counting = @elapsed unsupplied_cnt = [counting_unsupplied_subproblem(subproblems[s], 0.0001, data) for s in 1:data.S]

    invest_min = 0.0
    #invest_max = 100*objective_value(benders_master.model)
    invest_max = 1e12

    global best_obj = 0.0

    # Setting initial invest min and max
    if sum( data.probability .* unsupplied_cnt ) > max_unsupplied
        global invest_min = master_solution_investment_cost(benders_master, data)
    else
        global invest_max = master_solution_investment_cost(benders_master, data)
    end
    alpha = 0.5
    global rhs = 0.0


    print_log == true ? @printf("%-20.6e%-20.6e%-20.2e%-20.6e%-20.6e%-20.6e%-15.2f%-20.6e%-20.6e\n", invest_min, invest_max, 
            (invest_max - invest_min)/invest_max, rhs, investment_cost(benders_master, data),
            objective_value(benders_master.model), sum( data.probability .* unsupplied_cnt ), 
            t_optim, t_counting ) : nothing

    while invest_max - invest_min > max(relative_gap*invest_max,1e-6)

        global rhs = (1-alpha)*invest_max + alpha*invest_min

        set_normalized_rhs(constraint_by_name(benders_master.model, "invest_cost"), rhs)
    
        #global t_optim = @elapsed solve(stoch_prob, silent_mode)
        t_optim = @elapsed benders_sequential(benders_master, subproblems, data, 50)
        global t_counting = @elapsed global unsupplied_cnt = [counting_unsupplied_subproblem(subproblems[s], 0.0001, data) for s in 1:data.S]
        

        if sum( (data.probability .* unsupplied_cnt) ) > max_unsupplied
            global invest_min = rhs
        else
            global invest_max = rhs
            global best_obj = objective_value(benders_master.model)
        end

        print_log == true ? @printf("%-20.6e%-20.6e%-20.2e%-20.6e%-20.6e%-20.6e%-15.2f%-20.6e%-20.6e\n", invest_min, invest_max, 
            (invest_max - invest_min)/invest_max, rhs, master_solution_investment_cost(benders_master, data),
            objective_value(benders_master.model), sum( data.probability .* unsupplied_cnt ), 
            t_optim, t_counting ) : nothing
    
    end

    set_normalized_rhs(constraint_by_name(benders_master.model, "invest_cost"), 0)
    return invest_max, best_obj
end
#########################################################################################
#########################################################################################
#########################################################################################

#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################

#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################

#########################################################################################
# User options
#########################################################################################
N = 15
graph_density = 20
seed = 15
time_graph = @elapsed network = create_network(N, graph_density, seed, plotGraph = false, drawGraph = true)


seed >= 0 ? Random.seed!(seed) : nothing

scenarios = 2
time_steps = 5
demand_range = 100:1000
prod_cost_range = 10:40
unsupplied_cost = 500
epsilon_flow = 1.0
grad_prod = 0.2
invest_cost_range = 100:500
invest_prod_range = 100:500
flow_init_max = 500
time_data = @elapsed data = investment_problem_data_generator(scenarios, network, time_steps, demand_range, 
prod_cost_range, unsupplied_cost, epsilon_flow, flow_init_max, grad_prod, invest_cost_range, invest_prod_range)


# Benders problems creation
time_master_prob_creation = @elapsed benders_master = create_master_benders_problem(data)

println()
println("Subproblems creation")
println()

subproblems =  Vector{BendersSubroblem}(undef, scenarios)

epsilon_cnt = 0.001
counting_SPs =  Vector{BendersCountingSubroblem}(undef, scenarios)
for s in 1:data.S
    subproblems[s] = create_benders_subproblem(data, s)
    counting_SPs[s] = create_benders_subproblem_with_counting_unsupplied(data, s, epsilon_cnt)
    #println(subproblems[s].model)
end


benders_sequential(benders_master, subproblems, data, 1, "multicut")
exit()


max_unsupplied = 3
investment_heuristic_benders(benders_master, subproblems, data, max_unsupplied, 1e-6, true, true)

# Two phases resolution
# 1. LP relaxation
#unset_integer.(benders_master.invest_prod)

benders_sequential(benders_master, subproblems, data)

println(value.(benders_master.invest_prod))


# Evaluation du probleme stochastique
time_stoch_prob_creation = @elapsed stoch_prob = create_invest_optim_problem(data)

time_solve_stoch_prob = @elapsed solve(stoch_prob, false)

if termination_status(stoch_prob.model) == MOI.TIME_LIMIT
    println("TIME LIMIT reached")
    exit()
else
    val_stoch_prob = objective_value(stoch_prob.model)
    println("Problem value : ", val_stoch_prob)
    println("Optimal solution")
    println(value.(stoch_prob.invest_prod))
    for e in data.network.edges
        println(e, "  ", value(stoch_prob.invest_flow[e]))
    end
    println()
    println("*************************")
end

