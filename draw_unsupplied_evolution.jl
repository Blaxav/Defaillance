include("Heuristic.jl") # includes "Benders.jl"
include("Bilevel.jl")

#################################################
# Read option file
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




#################################################
# Read option file
#################################################
options_path = ARGS[1]
options = read_option_file(options_path)
print_options(options)

algo = create_algo(options)

println("Unsupplied evolution")
println()

#################################################
# Generate random datas
#################################################
time_graph = @elapsed network = create_network(options, plotGraph = false, drawGraph = false)
#seed >= 0 ? Random.seed!(seed) : nothing
time_data = @elapsed data = investment_problem_data_generator(options, network)
println()
println("Graph generation time = ", time_graph)
println("Problem data generation time = ", time_data)
println()

time_master_prob_creation = @elapsed master = create_master_benders_problem(data)
subproblems =  Vector{BendersSubroblem}(undef, data.S)
counting_SP =  Vector{BendersSubroblem}(undef, data.S)

println("Create counting subproblems : ", create_counting)
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

print_feasibility_log = false
perturbation = 0.0
println()
println("Feasibility check")
t_feasibility_check = feasibility_check(master, subproblems, h_data, 
    options, data, print_feasibility_log, perturbation)


iteration = 0
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
            print_log=false, log_freq=1, invest_free=false)
    
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