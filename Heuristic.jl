include("Benders.jl")

mutable struct HeuristicData
    invest_rhs::Float64
    UB_inv::Float64
    LB_inv::Float64
    tolerance::Float64
    max_unsupplied::Float64
    best_sol::BendersSolution
    found_solution_ite::Bool
    min_unsupplied::Float64
    counting_time::Float64
    total_time::Float64
    benders_iteration::Int64
end


function feasibility_check(master, subproblems, heuristic_data, options, data)

    # Set investment free objective function (investment cost set to perturbation value)
    set_invest_free_master_problem(master, data; perturbation = 0.0)

    best_solution = BendersSolution(edge_dict(), prod_dict(), 1e20, 0.0)
    t_feasibility = @elapsed benders_sequential(master, subproblems, data, algo, best_solution; 
        print_log=false, log_freq=1, invest_free=true)

    # Counting unsupplied 
    total_unsupplied = counting_unsupplied_total(subproblems, options.unsupplied_tolerance, data)

    if total_unsupplied > 0.0
        println("Infeasible instance. Minimum unsupplied = ", total_unsupplied)
        exit()
    else
        println("Feasilibity check passed.")
        heuristic_data.UB_inv = investment_cost(best_solution, data)
        heuristic_data.best_sol.flow = best_solution.flow
        heuristic_data.best_sol.prod = best_solution.prod
        heuristic_data.best_sol.val  = best_solution.val + investment_cost(best_solution, data)
        println("Max investment cost = ", heuristic_data.UB_inv)
        println()
    end

    # Reset original objective value
    set_initial_objective_master_problem(master, data)

    return t_feasibility
end



function counting_unsupplied__solution(h_data, separation_solution, subproblems, counting_SP, algo, options, data)

    rand_unsupplied_value = counting_unsupplied_total(subproblems, options.unsupplied_tolerance, data)

    if algo.heuristic_strategy == "Min" && 
        rand_unsupplied_value > options.max_unsupplied
        # Solve auxiliary SP
        separation_solution.unsupplied = 0.0
        fix_first_stage_candidate(counting_SP, separation_solution, data)
        for s in 1:data.S
            opt_val = get_objective_value(subproblems[s])
            set_normalized_rhs(counting_SP[s].cost_constraint, opt_val + 1e-3)

            solve(counting_SP[s]; silent_mode=true)
            separation_solution.unsupplied += data.probability[s] * get_objective_value(counting_SP[s])

            if separation_solution.unsupplied > options.max_unsupplied
                break
            end
        end
    else
        separation_solution.unsupplied = rand_unsupplied_value
    end
    
    h_data.min_unsupplied = min(rand_unsupplied_value, h_data.min_unsupplied)
end

function update_h_data_best_sol(solution, options, h_data, data)
    if solution.unsupplied <= options.max_unsupplied
        h_data.found_solution_ite = true
        if solution.val < h_data.best_sol.val
            h_data.best_sol.val = solution.val
            h_data.best_sol.flow = solution.flow
            h_data.best_sol.prod = solution.prod
        end
        if investment_cost(solution, data) < h_data.UB_inv
            h_data.UB_inv = investment_cost(solution, data)
        end
    end
end


function counting_benders(master, subproblems, counting_SP,
    data, algo, best_solution, h_data, max_unsupplied; 
    print_log=true, log_freq=1, invest_free=false, 
    count_frequency="Opt", count_strategy="Rand")

    separation_solution = BendersSolution(edge_dict(), prod_dict(), 0.0, 0.0)
    master_solution = BendersSolution(edge_dict(), prod_dict(), 0.0, 0.0)

    LB = 0.0
    iteration = 0
    cut_data = CutData(edge_dict(), prod_dict(), edge_dict(), prod_dict(), 0.0)
    stop = false
    h_data.counting_time    = 0.0
    h_data.min_unsupplied   = 0.0
    h_data.min_unsupplied = 1e20

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

        # 3. Fix investment candidates in subproblems and solve SPs
        fix_first_stage_candidate(subproblems, separation_solution, data)
        for s in 1:data.S
            solve(subproblems[s]; silent_mode=true)
        end
        compute_ub(subproblems, separation_solution, data; invest_free)

        # Compute value and unsupplied number
        if algo.heuristic_frequency == "All" && separation_solution.val < h_data.best_sol.val
            
            h_data.counting_time += @elapsed counting_unsupplied__solution(h_data, separation_solution, subproblems, counting_SP, algo, options, data)           
            update_h_data_best_sol(separation_solution, options, h_data, data)
        end

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

        # Test time limit
        if h_data.total_time + h_data.counting_time >= algo.time_limit
            println("Time limit exceeded. Time : ", h_data.total_time + h_data.counting_time)
            println("Best solution value = ", h_data.best_sol.val)
            println("############################")
            println("Solution")
            println("############################")
            print_solution(h_data.best_sol, data)
            exit()
        end
    end


    # Check unsupplied of benders_best_solution
    # As solution in SPs is not from the same solution, we need to solve them again
    fix_first_stage_candidate(subproblems, best_solution, data)
    for s in 1:data.S
        solve(subproblems[s]; silent_mode=true)
    end
    compute_ub(subproblems, best_solution, data; invest_free)

    # Then update heuristic data
    h_data.counting_time += @elapsed counting_unsupplied__solution(h_data, best_solution, subproblems, counting_SP, algo, options, data)           
    update_h_data_best_sol(best_solution, options, h_data, data)

    # If feasible for opt, then set RHS as UB_inv (IMPORTANT !!)
    if h_data.found_solution_ite
        h_data.UB_inv = h_data.invest_rhs
    end


    h_data.benders_iteration = iteration
end



function run_heuristic(options, data, algo)

    time_master_prob_creation = @elapsed master = create_master_benders_problem(data)
    subproblems =  Vector{BendersSubroblem}(undef, data.S)
    counting_SP =  Vector{BendersSubroblem}(undef, data.S)
    time_subproblems_creation = @elapsed create_all_subproblems(subproblems, counting_SP; create_counting=true)

    println("Time master problem creation : ", time_master_prob_creation)
    println("Time subproblem creation     : ", time_subproblems_creation)
    println()

    # Instanciate heuristic data
    h_data = HeuristicData(
        0.0, 1e20, 0.0, options.unsupplied_tolerance, options.max_unsupplied,
        BendersSolution(edge_dict(), prod_dict(), 1e20, 0.0), false,
        0.0, 0.0, 0.0, 0
    )

    # Feasibility check
    t_feasibility_check = feasibility_check(master, subproblems, h_data, options, data)

    iteration = 0
    # Initializing with stochastic solution
    @printf("%-10s%-20s%-20s%-20s%-12s%-20s%-20s%-20s%-15s%-15s%-15s%-15s\n", "Ite", "Best sol", 
            "Invest_min", "Invest_max", "Gap", "Constraint value", 
            "Invest cost", "Objective", "Unsupplied", "Optim time", "Count time", "Total time" )
    while (h_data.UB_inv - h_data.LB_inv)/h_data.UB_inv > 1e-6

        iteration += 1
        h_data.found_solution_ite = false

        # Solve
        benders_best_solution = BendersSolution(edge_dict(), prod_dict(), 1e20, 0.0)
        t_benders = @elapsed counting_benders(master, subproblems, counting_SP,
                        data, algo, benders_best_solution, h_data, options.max_unsupplied; 
                        print_log=false, log_freq=1, invest_free=false, 
                        count_frequency="Opt", count_strategy="Rand")
        
        invest_cost = investment_cost(benders_best_solution, data) 
        benders_val = benders_best_solution.val
        h_data.total_time += t_benders

        if h_data.found_solution_ite == false
            h_data.LB_inv = max(invest_cost, h_data.LB_inv)
        end   

        @printf("%-10i%-20.6e%-20.6e%-20.6e%-12.2e%-20.6e%-20.6e%-20.6e%-15.2f%-15.2e%-15.2e%-15.2e\n", iteration, 
                h_data.best_sol.val, h_data.LB_inv, 
                h_data.UB_inv, (h_data.UB_inv - h_data.LB_inv)/h_data.UB_inv, h_data.invest_rhs, 
                invest_cost, benders_val, h_data.min_unsupplied, t_benders, 
                h_data.counting_time, h_data.total_time
        )

        # Update investment constraint
        h_level = 0.5
        h_data.invest_rhs = h_level * h_data.UB_inv + (1-h_level) * h_data.LB_inv
        if h_data.UB_inv / h_data.LB_inv > 100
            h_data.invest_rhs = 10*h_data.LB_inv
        end
        set_normalized_rhs(master.heuristic_constraint, h_data.invest_rhs)

        if h_data.total_time >= algo.time_limit
            println("Time limit exceeded. Time : ", h_data.total_time + h_data.counting_time)
            println("Best solution value = ", h_data.best_sol.val)
            println("############################")
            println("Solution")
            println("############################")
            print_solution(h_data.best_sol, data)
            exit()
        end
    end

    println()
    println("Solution found")
    println("Best solution value = ", h_data.best_sol.val)
    println("############################")
    println("Solution")
    println("############################")
    print_solution(h_data.best_sol, data)

    return h_data.total_time
end