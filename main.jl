include("Benders.jl")
include("Bilevel.jl")

#################################################
# Read option file
#################################################
options_path = ARGS[1]
options = read_option_file(options_path)
print_options(options)

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

if options.algorithm == "bilevel"

    t_creation = @elapsed bilev = create_bilevel_invest_problem(data; unsupplied_tolerance=1e-3, max_unsupplied=3)
    t_solve = @elapsed solve(bilev; silent_mode=false)
    println()
    println("############################")
    println("Solution")
    println("############################")
    print_solution(bilev, data; null_tolerance=1e-6)

    println()
    println("Creation problem time = ", t_creation)
    println("Solution time = ", t_solve)

elseif options.algorithm == "stochastic"
elseif options.algorithm == "benders"
    total_time = @elapsed t_benders = run_benders(options, data)
    println()
    println("Solution time = ", t_benders)
    println("Total time = ", total_time)
elseif options.algorithm == "heuristic"
else
    println("Unknown algortihm ", options.algorithm)
end



