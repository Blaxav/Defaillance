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

mutable struct BendersSolution
    flow
    prod
    val
    unsupplied
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


function create_all_subproblems(subproblems, counting_SP; create_counting=false)
    for s in 1:data.S
        subproblems[s] = create_benders_subproblem(data, s)
        if create_counting
            counting_SP[s] = create_benders_subproblem_with_counting(data, s, options.unsupplied_tolerance)
        end
    end
end

function create_mean_value_prob(data)
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
    
    #Flow conservation
    constraint_flow_conservation_expextation(model, prod, flow, unsupplied, spilled, data)
    
    objective_mean_value_prob(model, invest_flow, invest_prod, 
        unsupplied, prod, flow_abs, data)
    
    BendersSubroblem(model, invest_flow, invest_prod, flow, prod, unsupplied, spilled, undef, undef)
end


################################################################
# Benders algorithm functions
################################################################
function get_master_solution(bendersMaster, master_solution)
    master_solution.flow = value.(bendersMaster.invest_flow)
    master_solution.prod = value.(bendersMaster.invest_prod)
end

function compute_separation_point(best_sol, master_sol, separation_sol, algo, data)

    for e in data.network.edges
        separation_sol.flow[e] = algo.step_size * master_sol.flow[e] + (1-algo.step_size) * best_sol.flow[e]
    end
    for n in production_nodes(data)
        separation_sol.prod[n] = algo.step_size * master_sol.prod[n] + (1-algo.step_size) * best_sol.prod[n]
    end
end


function fix_first_stage_candidate(subproblems, separation_sol, data)
    for s in 1:data.S
        # Can not broadcast fix on edges as they are modeled by a dictionary
        # Seems that broadcasting needs flat data
        for e in data.network.edges
            fix(subproblems[s].invest_flow[e], separation_sol.flow[e]; force=true)   
        end
        for n in production_nodes(data)
            fix.(subproblems[s].invest_prod[n], separation_sol.prod[n]; force=true)
        end
        #fix.(subproblems[s].invest_prod, candidate_prod; force=true)
    end
end


function investment_cost(solution, data)
    sum( data.invest_flow_cost[e] * solution.flow[e] for e in data.network.edges) +
    sum( data.invest_prod_cost[n] * solution.prod[n] for n in production_nodes(data))
end


function compute_ub(subproblems, separation_sol, data; invest_free=false)
    separation_sol.val = 0.0
    if invest_free == false
        separation_sol.val += investment_cost(separation_sol, data)
    end
    for s in 1:data.S
        separation_sol.val += data.probability[s] * get_objective_value(subproblems[s])
    end
end

function update_best_solution(best_sol, separation_sol, algo, LB)
    if separation_sol.val < best_sol.val - algo.stab_center_tol * (best_sol.val - LB)
        algo.step_size = min(1.0, 1.2*algo.step_size)

        # Updating best solution
        best_sol.val = separation_sol.val
        best_sol.flow = separation_sol.flow
        best_sol.prod = separation_sol.prod
        best_sol.unsupplied = separation_sol.unsupplied
    else
        algo.step_size = max(0.1, 0.8*algo.step_size)
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

function add_cuts(master, subproblems, separation_sol, data, cut_data, algo)
    if algo.cut_aggregation == "monocut"
        compute_monocut_gradient_and_rhs(cut_data, subproblems, data)
    
        @constraint(master.model, master.theta_sum >= 
            cut_data.rhs +
            sum( cut_data.grad_flow[e] * (master.invest_flow[e] - separation_sol.flow[e]) for e in data.network.edges) +
            sum( cut_data.grad_prod[n] * (master.invest_prod[n] - separation_sol.prod[n]) for n in production_nodes(data))
        )
    elseif algo.cut_aggregation == "multicut"
        for s in 1:data.S   
            compute_multicut_gradient_and_rhs(cut_data, subproblems, data, s)

            @constraint(master.model, master.theta[s] >= 
            cut_data.rhs +
                sum( cut_data.grad_flow[e] * (master.invest_flow[e] - separation_sol.flow[e]) for e in data.network.edges) +
                sum( cut_data.grad_prod[n] * (master.invest_prod[n] - separation_sol.prod[n]) for n in production_nodes(data))
                )
        end
    else
        println("Unknown algorithm ", algo)
        exit()
    end
end

"""
function counting_unsupplied_scenario
    brief: Computes the number of nodes which loss of load for a given scenario after optimization
"""
function counting_unsupplied_scenario(prob, epsilon, data)
    #return sum( value(prob.unsupplied[n,t]) > epsilon ? 1 : 0 for n in 1:data.network.N, t in 1:data.T )
    length( filter(x -> x > epsilon, value.(prob.unsupplied)) )
end


function counting_unsupplied_total(subproblems, tolerance, data)
    total_unsupplied = 0.0
    for s in 1:data.S
        total_unsupplied += data.probability[s] * counting_unsupplied_scenario(subproblems[s], tolerance, data)
    end
    return total_unsupplied
end

################################################################
# Benders algorithm
################################################################
function benders_sequential(master, subproblems, data, algo, best_solution; 
    print_log=true, log_freq=1, invest_free=false, check_heuristic=false)

    separation_solution = BendersSolution(edge_dict(), prod_dict(), 0.0, 0.0)
    master_solution = BendersSolution(edge_dict(), prod_dict(), 0.0, 0.0)

    LB = 0.0
    iteration = 0
    cut_data = CutData(edge_dict(), prod_dict(), edge_dict(), prod_dict(), 0.0)
    stop = false

    if print_log
        @printf("%-10s%-20s%-20s%-15s%15s\n", "ITER", "LB", "BEST UB", "GAP", "STEP SIZE")
    end

    while stop == false

        iteration += 1
        
        # 1. Solving master problem
        solve(master; silent_mode=true)
        
        # 2. Get master solution (investment candidates)
        get_master_solution(master, master_solution)
        LB = get_objective_value(master)

        if iteration == 1
            # Always cut master solution at iteration 1 if 0.0 is not feasible
            compute_separation_point(master_solution, master_solution, separation_solution, algo, data)
        else 
            compute_separation_point(best_solution, master_solution, separation_solution, algo, data)
        end
        # 3. Fix investment candidates in subproblems
        fix_first_stage_candidate(subproblems, separation_solution, data)
    
        
        # 4. Solve Subproblems
        for s in 1:data.S
            solve(subproblems[s]; silent_mode=true)
        end
        compute_ub(subproblems, separation_solution, data; invest_free)

        # Update UB
        update_best_solution(best_solution, separation_solution, algo, LB)
        
        # Stopping criterion
        if best_solution.val - LB <= 1e-6*best_solution.val
            stop = true
        end

        # Get subgradients and build cut
        # It is important to add cuts only if the algorithm continue
        # as adding cuts erase all information about any optimal solution got
        # by solving the master problem
        if stop == false
            add_cuts(master, subproblems, separation_solution, data, cut_data, algo)
        end

        # Log
        if print_log && ( stop == true || iteration % log_freq == 0 )
            @printf("%-10i%-20.6e%-20.6e%-15.3e%15.3f\n", iteration, LB, best_solution.val, (best_solution.val-LB)/best_solution.val, algo.step_size)
        end
    end
end


function initialize_mean_value_solution(master, subproblems, best_solution, data)
    
    # 1. Define and solve mean_value_problem
    mean_value_prob = create_mean_value_prob(data)
    solve(mean_value_prob; silent_mode=true) 

    # 2. Set solution to initialize best_sol
    get_master_solution(mean_value_prob, best_solution)
    
    # Fix value in SPs and solve
    fix_first_stage_candidate(subproblems, best_solution, data)
    for s in 1:data.S
        solve(subproblems[s]; silent_mode=true)
    end
    # Set true value of best_sol
    compute_ub(subproblems, best_solution, data)

    # Add the cuts to initialize master problem
    cut_data = CutData(edge_dict(), prod_dict(), edge_dict(), prod_dict(), 0.0)
    add_cuts(master, subproblems, best_solution, data, cut_data, algo)
end

function run_benders(options, data, algo)
    time_master_prob_creation = @elapsed benders_master = create_master_benders_problem(data)
    subproblems =  Vector{BendersSubroblem}(undef, data.S)
    for s in 1:data.S
        subproblems[s] = create_benders_subproblem(data, s)
    end

    t_benders = 0.0

    best_solution = BendersSolution(edge_dict(), prod_dict(), 1e20, 0.0)
    if algo.init_mean_value_solution
        t_benders += @elapsed initialize_mean_value_solution(benders_master, subproblems, best_solution, data)
    end

    t_benders += @elapsed benders_sequential(benders_master, subproblems, data, algo, best_solution; 
        print_log=true, log_freq=1)

    println()
    println("############################")
    println("Solution")
    println("############################")
    print_solution(best_solution, data; null_tolerance=1e-6)

    total_unsupplied = 0.0
    for s in 1:data.S
        total_unsupplied += data.probability[s] * counting_unsupplied_scenario(subproblems[s], 0.001, data)
    end
    println()
    println("Total unsupplied = ", total_unsupplied)

    return t_benders
end



################################################################
# Counting Subproblems
################################################################
function create_benders_subproblem_with_counting(data, s, tolerance)
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
    @constraint(model, unsupplied_to_zero[n in 1:data.network.N, t in 1:data.T],
        has_unsupplied[n,t] <= (1/tolerance)*unsupplied[n,t] )
    @constraint(model, unsupplied_to_one[n in 1:data.network.N, t in 1:data.T],
        100*data.scenario[s].demands[n,t]*has_unsupplied[n,t] >= unsupplied[n,t] - tolerance )

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

    return BendersSubroblem(model, invest_flow, invest_prod, flow, prod, unsupplied, spilled, cost_constraint, has_unsupplied)
end
