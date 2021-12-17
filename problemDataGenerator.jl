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

##################################################################################################
# Network investment problem class
##################################################################################################
"""
struct NetworkFlowProblemData
    brief: contains all the data to build a network flow problem on one scenario
"""
struct NetworkFlowProblemData
    network::Network
    demands::Vector{Float64}
    has_production::Vector{Int}
    unsupplied_cost::Float64
    prod_cost::Dict{Int,Float64}
    epsilon_flow::Float64
end



"""
function sample_network_data
    brief: Generates a random vector of NetworkFlowProblemData for each scenario
"""
function sample_network_data(scenarios, network, demand_range, 
    prod_cost_range, unsupplied_cost, epsilon_flow)

    data_flow = Vector{NetworkFlowProblemData}(undef, scenarios)
    for s in 1:scenarios

        demands = rand(demand_range, network.N)

        # Sampling production nodes
        has_production = rand(0:1, network.N)
        prod_cost = Dict(zip(
            [i for i in 1:N if has_production[i] == 1],
            rand(prod_cost_range, sum(has_production))
            ))

        data_flow[s] = NetworkFlowProblemData(network, demands, 
            has_production, unsupplied_cost, prod_cost, epsilon_flow)
    end
    return data_flow
end



##################################################################################################
# Main
##################################################################################################
if PROGRAM_FILE == "problemGenerator.jl"
    N = 15
    graph_density = 10
    seed = 0
    print("Generating graph ")
    @time network = create_network(N, graph_density, seed, plotGraph = false, drawGraph = true)


    seed >= 0 ? Random.seed!(seed) : nothing

    scenarios = 5
    proba = [0.6, 0.1, 0.1, 0.1, 0.1]
    print("Sample data      ")
    @time data_flow = sample_network_data(scenarios, network, 1:50, 1:10, 1000, 0.1)

    invest_flow_cost = Dict(zip(network.edges, rand(50:100, length(network.edges))))

    print("Create bilevel model     ")
    @time bilev, bi_invest, bi_flow, bi_prod, bi_unsupplied = crate_kkt_bilevel_problem(network, scenarios, proba, data_flow, invest_flow_cost)
    set_silent(bilev)
    optimize!(bilev)

    println("Obj bilevel : ", objective_value(bilev))
    println("Invest solution cost = ", sum([ invest_flow_cost[i,j] * value(bi_invest[(i,j)]) for (i,j) in network.edges]))    
    for (i,j) in network.edges
        if value(bi_invest[(i,j)]) != 0.0
            println("    invest ", (i,j), " = ", value(bi_invest[(i,j)]) )
        end
    end
    println("unsuppliedaillance")
    unsupplied_cnt = [ sum([ value(bi_unsupplied[s,n]) > 0 ? 1 : 0 for n in 1:network.N]) for s in 1:scenarios]
    for s in 1:scenarios
        println("    Scenario ", s, " n_unsupplied = ", unsupplied_cnt[s])
    end
    println("    unsupplied totale = ", sum( (proba .* unsupplied_cnt) ))
    println()
    

    ########################################
    # Problem formulation
    ########################################
    println()
    print("Create model     ")
    @time model = crate_invest_optim_problem(network, scenarios, proba, data_flow, invest_flow_cost)

    set_silent(model)
    print("Solving model    ")
    @time optimize!(model)
    
    ########################################
    # Solution analysis
    ########################################
    println()
    println("Obj : ", objective_value(model))

    println("Invest solution cost = ", sum([ invest_flow_cost[i,j] * value(variable_by_name(model,"invest_flow[$i,$j]")) for (i,j) in network.edges]))    
    for (i,j) in network.edges
        if value(variable_by_name(model,"invest_flow[$i,$j]")) != 0.0
            println("    invest ", (i,j), " = ", value(variable_by_name(model,"invest_flow[$i,$j]")) )
        end
    end
    println("unsuppliedaillance")
    unsupplied_cnt = [counting_unsupplied_scenario(model, s) for s in 1:scenarios]
    for s in 1:scenarios
        println("    Scenario ", s, " n_unsupplied = ", unsupplied_cnt[s])
    end
    println("    unsupplied totale = ", sum( (proba .* unsupplied_cnt) ))

    exit()

    ########################################
    # Heuristic
    ########################################
    println()
    println("Investment heuristic")
    @printf("%-15s%-15s%-15s%-10s\n", "Invest min", "Invest max", "Obj", "unsupplied count")

    invest_min = 0.0
    invest_max = 1e5
    while invest_max - invest_min > 1
        if sum( (proba .* unsupplied_cnt) ) > 0
            global invest_min = sum([ invest_flow_cost[i,j] * value(variable_by_name(model,"invest_flow[$i,$j]")) for (i,j) in network.edges])
        else
            global invest_max = sum([ invest_flow_cost[i,j] * value(variable_by_name(model,"invest_flow[$i,$j]")) for (i,j) in network.edges])
        end

        set_normalized_rhs(constraint_by_name(model, "invest_cost"), (invest_max + invest_min) / 2)

        optimize!(model)
        global unsupplied_cnt = [sum([ value(variable_by_name(model,"unsupplied[$s,$n]")) > 0 ? 1 : 0 for n in 1:network.N]) for s in 1:scenarios]
        @printf("%-15.3f%-15.3f%-15.3f%-10.3f\n", invest_min, invest_max, objective_value(model), sum( (proba .* unsupplied_cnt) ))
    end
end