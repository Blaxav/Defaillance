include("Benders.jl")

#################################################
# Read option file
#################################################
options = read_option_file("options.txt")



#################################################
# Generate random datas
#################################################
time_graph = @elapsed network = create_network(options, plotGraph = false, drawGraph = false)
#seed >= 0 ? Random.seed!(seed) : nothing
time_data = @elapsed data = investment_problem_data_generator(options, network)


if options.algorithm == "bilevel_sos1"
elseif options.algorithm == "stochastic"
elseif options.algorithm == "benders"
    total_time = @elapsed run_benders(options, data)
elseif options.algorithm == "heuristic"
else
    println("Unknown algortihm ", options.algorithm)
end

println("Graph generation time = ", time_graph)
println("Problem data generation time = ", time_data)
println("Total time = ", total_time)