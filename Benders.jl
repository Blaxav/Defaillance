include("optimization_base.jl")

################################################################
# Structures for master problem and subproblems
################################################################
struct BendersMasterProblem
    model::Model
    invest_flow
    invest_prod
    theta_sum # expectation of epigraph variables
    theta # epigraph variable for each scenario
    heuristic_constraint # Constraint on minimum cost investment for heuristic
end

# Structure for a subproblem in Benders decomposition
# cost_constraint and has_unsupplied can be set to "any"
# has_unsupplied : refers to binary variables saying if a node
#               has an unsupllied constraint or not at time t
# cost_constraint : constraint to set the cost in order to search 
#               only in optimal solutions (two-phase resolution)
struct BendersSubroblem
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

mutable struct HeuristicData
    invest_rhs::Float64
    UB_inv::Float64
    LB_inv::Float64
    epsilon::Float64
    alpha::Float64
    best_solution_cost::Float64
    best_flow_inv
    best_prod_inv
end


mutable struct CutData
    separation_flow
    separation_prod
    grad_flow
    grad_prod
    rhs
end


################################################################
# Problems constructors
################################################################
"""
create_benders_subproblem : Returns a subproblem of Benders decomposition
    data::StochasticProblemData     : Data of the problem
    s::Int64                        : Index of scenario
"""
function create_benders_subproblem(data, s)
    #model = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GRB_ENV), "Threads" => 1, "TimeLimit" => 600))
    model = Model(CPLEX.Optimizer)
    set_optimizer_attribute(model, "CPXPARAM_Threads", 1)

    # invest variables
    invest_flow = variables_investment_flow(model, data)
    invest_prod = variables_investment_production(model, data)

    # Flow variables
    flow = variables_flow(model, data)
    flow_abs = variables_absolute_value_flow(model, data)

    # Production variables
    prod = variables_production(model, data)
    
    # Production constraints
    constraint_max_prod(model, prod, invest_prod, data)
    
    constraint_production_positive_gradient(model, invest_prod, prod, data)
    constraint_production_negative_gradient(model, invest_prod, prod, data)
    # Loss of load variables
    unsupplied = variables_unsupplied_energy(model, data)
    spilled = variables_spilled_energy(model, data)

    # Flow bounds
    constraint_flow_max(model, flow, invest_flow, data)

    # Absolute value of flow in cost
    constraint_absolute_flow_behavior(model, flow, flow_abs)

    # Flow conservation
    constraint_flow_conservation(model, prod, flow, unsupplied, spilled, data; which_scenario = s)

    objective_subproblem(model, unsupplied, prod, flow_abs, data, which_scenario = s)

    return BendersSubroblem(model, invest_flow, invest_prod, flow, prod, unsupplied, spilled, undef, undef)
end


"""
create_master_benders_problem : Returns a Master problem of Benders decomposition
"""
function create_master_benders_problem(data)
    #model = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GRB_ENV), "Threads" => 1, "TimeLimit" => 600))
    model = Model(CPLEX.Optimizer)
    set_optimizer_attribute(model, "CPXPARAM_Threads", 1)

    # invest variables
    invest_flow = variables_investment_flow(model, data)
    invest_prod = variables_investment_production(model, data)
    
    # Constraint to lead heuristic
    invest_cost_ctr = constraint_minimum_investment(model, invest_flow, invest_prod, data)

    # Epigraph variables
    theta_sum = variables_epigraph_total(model, data)
    theta = variables_epigraph_on_scenario(model, data)
    constraint_total_epigraph(model, theta_sum, theta, data)
    
    objective_master_problem(model, invest_flow, invest_prod, theta_sum, data)

    return BendersMasterProblem(model, invest_flow, invest_prod, theta_sum, theta, invest_cost_ctr)
end


################################################################
# Benders algorithm functions
################################################################
function compute_separation_point(candidate_flow, candidate_prod, 
    separation_flow, separation_prod, alpha, data)

    for e in data.network.edges
        separation_flow[e] = alpha * candidate_flow[e] + (1-alpha) * best_invest_flow[e]
    end
    for n in 1:data.network.N 
        if data.scenario[1].has_production[n] == 1
            separation_prod[n] = alpha * candidate_prod[n] + (1-alpha) * best_invest_prod[n]
        end
    end
end

function compute_ub(subproblems, separation_flow, separation_prod, data; invest_free=false)
    UB = 0
    if invest_free == false
        UB += investment_cost(separation_flow, separation_prod, data)
    end
    for s in 1:data.S
        UB += data.probability[s] * get_objective_value(subproblems[s])
    end
    return UB
end

function update_best_solution(UB, best_UB, alpha, best_invest_flow, best_invest_prod)
    if UB < best_UB
        global alpha = min(1.0, 1.2*alpha)
        global best_UB = UB
        global best_invest_prod = separation_prod
        global best_invest_flow = separation_flow
    else
        global alpha = max(0.1, 0.8*alpha)
    end
end


function compute_monocut_gradient_and_rhs(cut_data, subproblems, data)
    cut_data.rhs = 0.0
    for s in 1:data.S   
        cut_data.rhs += data.probability[s] * get_objective_value(subproblems[s])

        if s == 1
            cut_data.grad_flow = data.probability[s] * reduced_cost.(subproblems[s].invest_flow)
            cut_data.grad_prod = data.probability[s] * reduced_cost.(subproblems[s].invest_prod)
        else
            for e in data.network.edges
                cut_data.grad_flow[e] += data.probability[s] * reduced_cost(subproblems[s].invest_flow[e])
            end
            for n in production_nodes(data)
                cut_data.grad_prod[n] += data.probability[s] * reduced_cost(subproblems[s].invest_prod[n])
            end
        end
    end
end

function compute_multicut_gradient_and_rhs(cut_data, subproblems, data, s)
    cut_data.rhs = get_objective_value(subproblems[s])
    cut_data.grad_flow = reduced_cost.(subproblems[s].invest_flow)
    cut_data.grad_prod = reduced_cost.(subproblems[s].invest_prod)
end

function add_cuts(master, subproblems, data, cut_data; algo="monocut")
    if algo == "monocut"
        compute_monocut_gradient_and_rhs(cut_data, subproblems, data)
    
        @constraint(master.model, master.theta_sum >= 
            cut_data.rhs +
            sum( cut_data.grad_flow[e] * (master.invest_flow[e] - separation_flow[e]) for e in data.network.edges) +
            sum( cut_data.grad_prod[n] * (master.invest_prod[n] - separation_prod[n]) for n in production_nodes(data))
        )
    elseif algo == "multicut"
        for s in 1:data.S   
            compute_multicut_gradient_and_rhs(cut_data, subproblems, data, s)

            @constraint(master.model, master.theta[s] >= 
            cut_data.rhs +
                sum( cut_data.grad_flow[e] * (master.invest_flow[e] - separation_flow[e]) for e in data.network.edges) +
                sum( cut_data.grad_prod[n] * (master.invest_prod[n] - separation_prod[n]) for n in production_nodes(data))
                )
        end
    else
        println("Unknown algorithm ", algo)
        exit()
    end
end



################################################################
# Benders algorithm
################################################################
function benders_sequential(master, subproblems, data, print_log, n_iteration_log, algo, invest_free, heuristic_data, check_heuristic)

    global best_UB = 1e20
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
    global alpha = 1.0

    local cut_data = CutData(
        undef, undef, 
        Dict(zip(
            data.network.edges,
            [0.0 for e in data.network.edges])),
        Dict(zip(
            [i for i in production_nodes(data)],
            [0.0 for i in 1:sum(data.has_production)])),
        0.0 )

    global log_freq = n_iteration_log
    global stop = false

    if print_log
        @printf("%-10s%-20s%-20s%-15s%15s\n", "ITER", "LB", "BEST UB", "GAP", "STEP SIZE")
    end

    while stop == false

        iteration += 1
        
        # 1. Solving master problem
        solve(master; silent_mode=true)
        
        # 2. Get master solution (investment candidates)
        candidate_flow, candidate_prod = get_master_solution(master)
        LB = get_objective_value(master)

        compute_separation_point(candidate_flow, candidate_prod, 
            separation_flow, separation_prod, alpha, data)

        
        separation_prod[1] = 382.0
        separation_prod[2] = 2200.50
        separation_prod[4] = 308.50

        for e in data.network.edges
            separation_flow[e] = 0.0
        end
        separation_flow[Edge(2,5)] = 865.0
        separation_flow[Edge(6,9)] = 299.0
        separation_flow[Edge(2,10)] = 663.50
        separation_flow[Edge(4,10)] = 71.50
        separation_flow[Edge(8,9)] = 201.0
        separation_flow[Edge(5,6)] = 612.0
        separation_flow[Edge(3,10)] = 390.0
        separation_flow[Edge(2,7)] = 396.0

        # 3. Fix investment candidates in subproblems
        fix_first_stage_candidate(subproblems, separation_flow, separation_prod, data)
    
        
        # 4. Solve Subproblems and get subprgradients
        true_unsup = 0.0
        SP_cnt =  Vector{BendersSubroblem}(undef, data.S)
        for s in 1:data.S
            SP_cnt[s] = create_benders_subproblem_with_counting(data, s)
        end
        fix_first_stage_candidate(SP_cnt, separation_flow, separation_prod, data)

        for s in 1:data.S
            solve(subproblems[s]; silent_mode=true)
            cost = objective_value(subproblems[s].model)
            println("Cost = ", cost, "  SP ", s)
            
            # Solving counting SPs
            set_normalized_rhs(SP_cnt[s].cost_constraint, cost + 1e-6)
            solve(SP_cnt[s]; silent_mode=false)
            true_unsup += data.probability[s] * objective_value(SP_cnt[s].model)
        end
        println("Min Unsupplied = ", true_unsup)
        local UB = compute_ub(subproblems, separation_flow, separation_prod, data; invest_free)


        unsupplied_cnt = sum([ data.probability[s] * counting_unsupplied_scenario(subproblems[s], 0.001, data) for s in 1:data.S ])
        println("Unsupplied = ", unsupplied_cnt)
        println("Val = ", UB)
        exit()

        # Check if investment solution is bilevel feasible
        if check_heuristic == true
            unsupplied_cnt = sum([ data.probability[s] * counting_unsupplied_scenario(subproblems[s], heuristic_data.epsilon, data) for s in 1:data.S ])
            if unsupplied_cnt <= heuristic_data.alpha
                if investment_cost(separation_flow, separation_prod, data) < heuristic_data.UB_inv
                    heuristic_data.UB_inv = investment_cost(separation_flow, separation_prod, data)
                end
                if UB < heuristic_data.best_solution_cost
                    heuristic_data.best_solution_cost = UB
                    heuristic_data.best_flow_inv = separation_flow
                    heuristic_data.best_prod_inv = separation_prod
                end
            end
        end


        # Update UB
        update_best_solution(UB, best_UB, alpha, best_invest_flow, best_invest_prod)
        
        # Stopping criterion
        if best_UB - LB <= 1e-6*best_UB
            global stop = true
        end

        # Get subgradients and build cut
        # It is important to add cuts only if the algorithm continue
        # as adding cuts erase all information about any optimal solution got
        # by solving the master problem
        if stop == false
            add_cuts(master, subproblems, data, cut_data; algo=algo)
        end

        # Log
        if print_log && ( stop == true || iteration % log_freq == 0 )
            @printf("%-10i%-20.6e%-20.6e%-15.3e%15.3f\n", iteration, LB, best_UB, (best_UB-LB)/best_UB, alpha)
        end
    end
end


function run_benders(options, data)
    time_master_prob_creation = @elapsed benders_master = create_master_benders_problem(data)
    subproblems =  Vector{BendersSubroblem}(undef, data.S)
    for s in 1:data.S
        subproblems[s] = create_benders_subproblem(data, s)
    end

    t_benders = @elapsed benders_sequential(benders_master, subproblems, data, true, 1, "multicut", false, undef, false)

    println()
    println("############################")
    println("Solution")
    println("############################")
    print_solution(benders_master, data; null_tolerance=1e-6)

    return t_benders
end





#=function create_benders_subproblem_with_counting_unsupplied(data, s, epsilon_cnt)
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
=#

function create_benders_subproblem_with_counting(data, s)
    #model = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GRB_ENV), "Threads" => 1, "TimeLimit" => 600))
    model = Model(CPLEX.Optimizer)
    set_optimizer_attribute(model, "CPXPARAM_Threads", 1)

    # invest variables
    invest_flow = variables_investment_flow(model, data)
    invest_prod = variables_investment_production(model, data)

    # Flow variables
    flow = variables_flow(model, data)
    flow_abs = variables_absolute_value_flow(model, data)

    # Production variables
    prod = variables_production(model, data)
    
    # Production constraints
    constraint_max_prod(model, prod, invest_prod, data)
    
    constraint_production_positive_gradient(model, invest_prod, prod, data)
    constraint_production_negative_gradient(model, invest_prod, prod, data)
    # Loss of load variables
    unsupplied = variables_unsupplied_energy(model, data)
    spilled = variables_spilled_energy(model, data)

    # Flow bounds
    constraint_flow_max(model, flow, invest_flow, data)

    # Absolute value of flow in cost
    constraint_absolute_flow_behavior(model, flow, flow_abs)

    # Flow conservation
    constraint_flow_conservation(model, prod, flow, unsupplied, spilled, data; which_scenario = s)


    # Variable saying if a node has unsupplied energy
    @variable(model, has_unsupplied[i in 1:data.network.N, t in 1:data.T], Bin)

    # Binary Variables behaviour
    epsilon_cnt = 1.0
    @constraint(model, unsupplied_to_zero[n in 1:data.network.N, t in 1:data.T],
        has_unsupplied[n,t] <= (1/epsilon_cnt)*unsupplied[n,t] )
    @constraint(model, unsupplied_to_one[n in 1:data.network.N, t in 1:data.T],
        100*data.scenario[s].demands[n,t]*has_unsupplied[n,t] >= unsupplied[n,t] - epsilon_cnt )

    @constraint(model, cost_constraint, 
        sum(
            (
            # unsupplied costs
            sum( data.scenario[s].unsupplied_cost * unsupplied[n,t] 
                for n in 1:data.network.N ) +
            # production costs
            sum( data.scenario[s].prod_cost[n][t] * prod[n,t] 
                for n in production_nodes(data)) +
            # flow cost
            sum( data.scenario[s].flow_cost[e] * flow_abs[e,t] 
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
    #objective_subproblem(model, unsupplied, prod, flow_abs, data, which_scenario = s)

    return BendersSubroblem(model, invest_flow, invest_prod, flow, prod, unsupplied, spilled, cost_constraint, has_unsupplied)
end





function counting_unsupplied_subproblem(stoch_prob, epsilon, data)
    return sum( value(stoch_prob.unsupplied[n,t]) > epsilon ? 1 : 0 for n in 1:data.network.N, t in 1:data.T )
end
"""
function counting_unsupplied_scenario
    brief: Computes the number of nodes which loss of load for a given scenario after optimization
"""
function counting_unsupplied_scenario(prob, epsilon, data)
    #return sum( value(prob.unsupplied[n,t]) > epsilon ? 1 : 0 for n in 1:data.network.N, t in 1:data.T )
    length( filter(x -> x > epsilon, value.(prob.unsupplied)) )
end

#########################################################################################
# Benders resolution -- OOOOLLLLLLDDDD
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

    investment_cost(solution_invest_flow, solution_invest_prod, data)
end


#=function benders_sequential(master, subproblems, data, print_log, n_iteration_log, algo, invest_free, heuristic_data, check_heuristic)

    global best_UB = 1e20
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
    global alpha = 1.0


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
        solve(master; silent_mode=true)
        
        # 2. Get master solution (investment candidates)
        candidate_flow, candidate_prod = get_master_solution(master)
        LB = get_objective_value(master)

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
        local UB = 0.0
        if invest_free == false
            UB += investment_cost(separation_flow, separation_prod, data)
        end
        for s in 1:data.S
            solve(subproblems[s]; silent_mode=true)
            UB += data.probability[s] * get_objective_value(subproblems[s])
        end


        # Check if investment solution is bilevel feasible
        if check_heuristic == true
            unsupplied_cnt = sum([ data.probability[s] * counting_unsupplied_scenario(subproblems[s], heuristic_data.epsilon, data) for s in 1:data.S ])
            if unsupplied_cnt <= heuristic_data.alpha
                #println("Bilevel feasible solution found at cost = ", UB)
                #println("Invest cost = ", investment_cost(separation_flow, separation_prod, data))
                #println("Total cost  = ", UB)
                if investment_cost(separation_flow, separation_prod, data) < heuristic_data.UB_inv
                    heuristic_data.UB_inv = investment_cost(separation_flow, separation_prod, data)
                end
                if UB < heuristic_data.best_solution_cost
                    heuristic_data.best_solution_cost = UB
                    heuristic_data.best_flow_inv = separation_flow
                    heuristic_data.best_prod_inv = separation_prod
                end
            end
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
        
        # Stopping criterion
        if best_UB - LB <= 1e-6*best_UB
            global stop = true
        end


        # Get subgradients and build cut
        if stop == false

            local rhs = 0.0

            if algo == "monocut"
                for s in 1:data.S   
                    rhs += data.probability[s] * get_objective_value(subproblems[s])
    
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
                    rhs = get_objective_value(subproblems[s])
    
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

        # Log
        if print_log && ( stop == true || iteration % log_freq == 0 )
            @printf("%-10i%-20.6e%-20.6e%-15.3e%15.3f\n", iteration, LB, best_UB, (best_UB-LB)/best_UB, alpha)
        end
    end
end
=#

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
if PROGRAM_FILE == "Benders.jl"
    N = 30
    graph_density = 10
    seed = 3
    time_graph = @elapsed network = create_network(N, graph_density, seed, plotGraph = false, drawGraph = true)


    seed >= 0 ? Random.seed!(seed) : nothing

    #=scenarios = 50
    time_steps = 24
    demand_range = 10:80
    prod_cost_range = 10:20
    unsupplied_cost = 300
    epsilon_flow = 0.1
    grad_prod = 0.3
    invest_cost_range = 100:500
    invest_prod_range = 50:200
    flow_init_max = 15=#
    scenarios = 10
    time_steps = 20
    demand_range = 100:500
    #prod_cost_range = 300:800
    prod_cost_range = 30:80
    unsupplied_cost = 1000
    epsilon_flow = 1.0
    grad_prod = 0.2
    invest_cost_range = 500:1000
    invest_prod_range = 500:1000
    flow_init_max = 200

    time_data = @elapsed data = investment_problem_data_generator(scenarios, network, time_steps, demand_range, 
    prod_cost_range, unsupplied_cost, epsilon_flow, flow_init_max, grad_prod, invest_cost_range, invest_prod_range)

    # Benders problems creation
    time_master_prob_creation = @elapsed benders_master = create_master_benders_problem(data)

    println()
    println("Subproblems creation")
    println()

    subproblems =  Vector{BendersSubroblem}(undef, scenarios)


    counting_SPs =  Vector{BendersSubroblem}(undef, scenarios)
    for s in 1:data.S
        subproblems[s] = create_benders_subproblem(data, s)
        #counting_SPs[s] = create_benders_subproblem_with_counting_unsupplied(data, s, epsilon_cnt)
        #println(subproblems[s].model)
    end


    # Heuristic Parameters
    global LB_inv = 0.0
    global UB_inv = 5e8 # a modifier with heuristic
    global heur_param = 0.2
    global Lambda_inv = 0.0
    global epsilon = 1e-3
    global max_unsupplied_ctr = 1.0
    heuristic_data = HeuristicData(Lambda_inv, UB_inv, LB_inv, epsilon, max_unsupplied_ctr, 1e20, undef, undef)

    ##### Feasibility check 
    #set_invest_free_master(benders_master, data, 1e-3)
    set_invest_free_master_problem(benders_master, data, 1e-3)
    #print(benders_master.model)

    println("#####################################")
    println("         Feasibility check           ")
    t_benders = @elapsed benders_sequential(benders_master, subproblems, data, true, 10, "multicut", true, heuristic_data, false)
    t_counting = @elapsed global unsupplied_cnt = sum([ data.probability[s] * counting_unsupplied_scenario(subproblems[s], epsilon, data) for s in 1:data.S ])

    println("Time feasibility check : ", t_benders + t_counting)
    if unsupplied_cnt > 0
        println("INFEASIBLE PROBLEM, minimum unsupplied  = ", unsupplied_cnt)
        exit()
    else
        UB_inv = master_solution_investment_cost(benders_master, data)
        println("Max invest cost = ", UB_inv)
    end
    println("#####################################")
    println()

    # Reset objective function to initial one
    set_initial_objective_master_problem(benders_master, data)



    t_benders = @elapsed benders_sequential(benders_master, subproblems, data, true, 10, "multicut", false, heuristic_data, true)

    println("Stochastic problem solution")
    println()
    invest_prod_sol = value.(benders_master.invest_prod)
    invest_flow_sol = value.(benders_master.invest_flow)
    for n in 1:data.network.N
    if data.has_production[n] == 1 &&  invest_prod_sol[n] > 1e-6
        @printf("%-10i%-10.3f\n", n, invest_prod_sol[n])
    end
    end
    for e in data.network.edges
        if invest_flow_sol[e] > 1e-6
            println(e, "   ", invest_flow_sol[e])
        end
    end

    t_counting = @elapsed global unsupplied_cnt = sum([ data.probability[s] * counting_unsupplied_scenario(subproblems[s], epsilon, data) for s in 1:data.S ])
    if unsupplied_cnt <= heuristic_data.alpha
        heuristic_data.UB_inv = heuristic_data.invest_rhs
    else
        heuristic_data.LB_inv = master_solution_investment_cost(benders_master, data)
    end

    @printf("%-20s%-20s%-20s%-20s%-20s%-20s%-20s%-15s%-20s%-20s\n", "Best sol", "Invest_min", "Invest_max", "Gap", "Constraint value", 
                "Invest cost", "Objective", "Unsupplied", "Optim time", "Counting time" )

    global invest_cost = master_solution_investment_cost(benders_master, data)
    global benders_val = get_objective_value(benders_master)
    @printf("%-20.6e%-20.6e%-20.6e%-20.2e%-20.6e%-20.6e%-20.6e%-15.2f%-20.6e%-20.6e\n", heuristic_data.best_solution_cost, heuristic_data.LB_inv, 
                heuristic_data.UB_inv, (heuristic_data.UB_inv - heuristic_data.LB_inv)/heuristic_data.UB_inv, 
                heuristic_data.invest_rhs, invest_cost, benders_val, unsupplied_cnt, t_benders, t_counting
            )

    while heuristic_data.UB_inv - heuristic_data.LB_inv > 1e-3*heuristic_data.UB_inv

        # Compute Lambda
        heuristic_data.invest_rhs = heur_param*heuristic_data.UB_inv + (1-heur_param)*heuristic_data.LB_inv

        # Modify RHS of investment cost constraint
        #set_normalized_rhs(constraint_by_name(benders_master.model, "invest_cost"), heuristic_data.invest_rhs)
        set_normalized_rhs(benders_master.heuristic_constraint, heuristic_data.invest_rhs)

        # Solve
        local t_benders = @elapsed benders_sequential(benders_master, subproblems, data, false, 10, "multicut", false, heuristic_data, false)
        global invest_cost = master_solution_investment_cost(benders_master, data)
        global benders_val = get_objective_value(benders_master)

        # Modify LB or UB
        local t_counting = @elapsed global unsupplied_cnt = sum([ data.probability[s] * counting_unsupplied_scenario(subproblems[s], epsilon, data) for s in 1:data.S ])
        if unsupplied_cnt <= heuristic_data.alpha
            if heuristic_data.UB_inv > master_solution_investment_cost(benders_master, data) 
                heuristic_data.UB_inv = master_solution_investment_cost(benders_master, data)
            end
            if get_objective_value(benders_master) < heuristic_data.best_solution_cost
                heuristic_data.best_solution_cost = get_objective_value(benders_master)
                #heuristic_data.best_invest_flow = value.(benders_master.invest_flow)
                #heuristic_data.best_invest_prod = value.(benders_master.invest_prod)
            end
        else
            heuristic_data.LB_inv = heuristic_data.invest_rhs
        end
        
        @printf("%-20.6e%-20.6e%-20.6e%-20.2e%-20.6e%-20.6e%-20.6e%-15.2f%-20.6e%-20.6e\n", heuristic_data.best_solution_cost, heuristic_data.LB_inv, 
                heuristic_data.UB_inv, (heuristic_data.UB_inv - heuristic_data.LB_inv)/heuristic_data.UB_inv, 
                heuristic_data.invest_rhs, invest_cost, benders_val, unsupplied_cnt, t_benders, t_counting
            )
        
    end

    println(value.(benders_master.invest_prod))
    println(value.(benders_master.invest_flow))


    bilev = create_bilevel_invest_problem(data, epsilon, max_unsupplied_ctr)
    solve(bilev, false)
end