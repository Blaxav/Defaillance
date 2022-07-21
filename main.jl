include("Benders.jl")
include("Bilevel.jl")

#if "--compile" in ARGS
##    include("small_exec_compiler.jl")
#end 

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
    
    #=fix(bilev.invest_prod[1], 382.0; force=true)
    fix(bilev.invest_prod[2], 2237.498; force=true)
    fix(bilev.invest_prod[4], 374.50; force=true)

    for e in data.network.edges
        fix(bilev.invest_flow[e], 0.0; force=true)
    end
    fix(bilev.invest_flow[Edge(2,5)], 1064.0; force=true)
    fix(bilev.invest_flow[Edge(6,9)], 380.0; force=true)
    fix(bilev.invest_flow[Edge(2,10)], 592.50; force=true)
    fix(bilev.invest_flow[Edge(4,10)], 5.50; force=true)
    fix(bilev.invest_flow[Edge(8,9)], 201.0; force=true)
    fix(bilev.invest_flow[Edge(5,6)], 715.0; force=true)
    fix(bilev.invest_flow[Edge(3,10)], 381.0; force=true)
    fix(bilev.invest_flow[Edge(2,7)], 396.0; force=true)

    t_solve = @elapsed solve(bilev; silent_mode=false)

    println()
    println()
    println()

    unfix(bilev.invest_prod[1])
    unfix(bilev.invest_prod[2])
    unfix(bilev.invest_prod[4])

    for e in data.network.edges
        unfix(bilev.invest_flow[e])
    end=#
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



