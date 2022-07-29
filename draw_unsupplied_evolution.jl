include("Heuristic.jl") # includes "Benders.jl"
include("Bilevel.jl")
using Distributions


struct Unsupplied
    invest
    min
    rand
    max
end


#################################################
# Functions
#################################################
function set_max_unsupplied_objective(model, has_unsupplied, data)
    @objective(model, Max,
        # Sum on time steps
        sum(# Number of unsupplied nodes
            ( 
                sum(has_unsupplied[n,t] for n in 1:data.network.N )
            ) for t in 1:data.T
        )
    )
end

function set_min_unsupplied_objective(model, has_unsupplied, data)
    @objective(model, Min,
        # Sum on time steps
        sum(# Number of unsupplied nodes
            ( 
                sum(has_unsupplied[n,t] for n in 1:data.network.N )
            ) for t in 1:data.T
        )
    )
end


function compute_std_sample(mean, min_val, proba)
    d = Normal(0, 1)
    return sigma = (min_val - mean) / quantile(d, 1 - proba) 
end
  

function initial_heuristic(h_data, master, subproblems, counting_SP, data, algo)

    iteration = 0
    @printf("%-10s%-20s%-20s%-20s%-12s%-20s%-20s%-20s%-15s%-15s%-15s%-15s\n", "Ite", "Best sol", 
            "Invest_min", "Invest_max", "Gap", "Constraint value", 
            "Invest cost", "Objective", "Unsupplied", "Optim time", "Count time", "Total time" )

    min_value = 0.0
    max_value = h_data.UB_inv


    while (h_data.UB_inv - h_data.LB_inv)/h_data.UB_inv > 1e-6

        iteration += 1
        h_data.found_solution_ite = false

        # Solve
        benders_best_solution = BendersSolution(edge_dict(), prod_dict(), 1e20, 0.0)
        t_benders = @elapsed counting_benders(master, subproblems, counting_SP,
                data, algo, benders_best_solution, h_data, options.max_unsupplied; 
                print_log=false, log_freq=1, invest_free=false)
        
        invest_cost = investment_cost(benders_best_solution, data)
        
        # Get stochastic value as minimum invest val
        if iteration == 1
            min_value = invest_cost
        end

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
            if h_data.invest_rhs == 0.0
                h_data.invest_rhs = 0.1 * h_data.UB_inv
            end
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

    return min_value
end

function compute_unsupplied_evolution(vect_invest, vect_rand, vect_min, vect_max, 
    master,subproblems, counting_SP, h_data, data, algo, options)

    unsupplied_variables = zeros(Float64, data.network.N, data.T)
    iteration = 0

    @printf("%-10s%-15s%-15s%-15s%-15s\n", "It", "Val", "Min", "Rand", "Max")
    for invest_val in vect_invest
        iteration += 1

        # Set sampled value as investment RHS
        set_normalized_rhs(master.heuristic_constraint, invest_val)
        
        # Solve Benders
        benders_best_solution = BendersSolution(edge_dict(), prod_dict(), 1e20, 0.0)
        t_benders = @elapsed counting_benders(master, subproblems, counting_SP,
                data, algo, benders_best_solution, h_data, options.max_unsupplied; 
                print_log=false, log_freq=1, invest_free=false)

        
        # Check unsupplied of benders_best_solution
        fix_first_stage_candidate(subproblems, benders_best_solution, data)
        for s in 1:data.S
            solve(subproblems[s]; silent_mode=true)
        end
        compute_ub(subproblems, benders_best_solution, data; invest_free=false)

        # Counting Rand, Min and Max val
        rand_unsupplied = counting_unsupplied_total(subproblems, options.unsupplied_tolerance, 
                                data, unsupplied_variables)

        # Solve auxiliary SP
        min_unsupplied = 0.0
        max_unsupplied = 0.0
        absolute_perturbation = 1e-6
        relative_perturbation = 1e-8
        fix_first_stage_candidate(counting_SP, benders_best_solution, data)
        for s in 1:data.S
            opt_val = get_objective_value(subproblems[s])
            set_normalized_rhs(counting_SP[s].cost_constraint, (1 + relative_perturbation) * opt_val + absolute_perturbation)

            set_min_unsupplied_objective(counting_SP[s].model, counting_SP[s].has_unsupplied, data)
            #println("Min ", s)
            solve(counting_SP[s]; silent_mode=true)
            min_unsupplied += data.probability[s] * get_objective_value(counting_SP[s])

            set_max_unsupplied_objective(counting_SP[s].model, counting_SP[s].has_unsupplied, data)
            #println("Max ", s)
            solve(counting_SP[s]; silent_mode=true)
            max_unsupplied += data.probability[s] * get_objective_value(counting_SP[s])
        end
        @printf("%-10i%-15.3e%-15.3f%-15.3f%-15.3f\n", iteration, invest_val, min_unsupplied, rand_unsupplied, max_unsupplied)

        #if first_reached == 0.0 && min_unsupplied <= 3.0
        #    first_reached = invest_val
        #end

        push!(vect_rand, rand_unsupplied)
        push!(vect_min, min_unsupplied)
        push!(vect_max, max_unsupplied)
    end
end


function main()

    time_master_prob_creation = @elapsed master = create_master_benders_problem(data)
    subproblems =  Vector{BendersSubroblem}(undef, data.S)
    counting_SP =  Vector{BendersSubroblem}(undef, data.S)
    time_subproblems_creation = @elapsed create_all_subproblems(subproblems, counting_SP; create_counting=true)


    println("Time master problem creation : ", time_master_prob_creation)
    println("Time subproblem creation     : ", time_subproblems_creation)
    println()
    print_problem_summary(master, subproblems, data)

    # Instanciate heuristic data
    h_data = HeuristicData(
        0.0, 1e20, 0.0, options.unsupplied_tolerance, options.max_unsupplied,
        BendersSolution(edge_dict(), prod_dict(), 1e20, 0.0), false,
        0.0, 0.0, 0.0, 0, zeros(Float64, data.network.N, data.T)
    )

    #=print_feasibility_log = false
    perturbation = 0.0
    println()
    println("Feasibility check")
    t_feasibility_check = feasibility_check(master, subproblems, h_data, 
        options, data, print_feasibility_log, perturbation)=#

    min_value = initial_heuristic(h_data, master, subproblems, counting_SP, data, algo)

    println("Heuristic finished")
    mean = investment_cost(h_data.best_sol, data)
    println("Mean value of sampling   = ", mean)
    println("Minimum investment value = ", min_value)
    
    n_samples = 500
    coef_high = 1.05
    proba = 0.01
    println("Sampling ", n_samples, " random investments")
    println("Computing standard deviation such that the probability to be higher than ", coef_high, " times the expected value is lower than ", proba)
    println("")
    sigma = compute_std_sample(mean, coef_high*mean, proba)
    println("Sigma = ", sigma)

    d = Normal(mean, sigma)
    td = truncated(d, min_value, 100*min_value)

    println()
    
    
    sampled_invest = sort(rand(td, n_samples))
    
    #vect_invest = sampled_invest

    # Adding a linear sample from min to max
    N = 100
    vect_invest = [min_value + (k/N)*(1.3*mean - min_value) for k in 0:N]
    vect_invest = sort([vect_invest; sampled_invest])

    

    rand_unsupplied = 0.0
    min_unsupplied = 0.0
    max_unsupplied = 0.0

    vect_min    = Vector{Float64}()
    vect_rand   = Vector{Float64}()
    vect_max    = Vector{Float64}()
    compute_unsupplied_evolution(vect_invest, vect_rand, vect_min, vect_max, 
        master,subproblems, counting_SP, h_data, data, algo, options)

    
    graph_plot_final = plot(vect_invest, vect_rand, labels="Random")
    plot!(vect_invest, vect_min, labels="Minimum")
    plot!(vect_invest, vect_max, labels="Maximum")
    
    display(graph_plot_final)
    println("Press a key to continue")
    readline()
end


options_path = ARGS[1]
options = read_option_file(options_path)
print_options(options)

algo = create_algo(options)

println("Unsupplied evolution")
println()



#################################################
# Call Main function
#################################################
time_graph = @elapsed network = create_network(options, plotGraph = false, drawGraph = false)
#seed >= 0 ? Random.seed!(seed) : nothing
time_data = @elapsed data = investment_problem_data_generator(options, network)
println()
println("Graph generation time = ", time_graph)
println("Problem data generation time = ", time_data)
println()




main()
