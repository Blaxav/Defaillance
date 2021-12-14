#########################################################################################
# Import management
#########################################################################################
if "--install-packages" in ARGS
    import Pkg
    Pkg.add("Random")
    Pkg.add("Plots")
    Pkg.add("JuMP")
    Pkg.add("CPLEX")
    Pkg.build("CPLEX")
end

try
    using Random
    using Plots
    using JuMP, CPLEX
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
struct NetworkFlowProblemData
    network::Network
    demands::Vector{Float64}
    has_production::Vector{Int}
    def_cost::Int
    prod_cost::Dict{Int,Float64}
    epsilon_flow::Float64
end

function sample_network_data(scenarios, network, demand_range, 
    prod_cost_range, def_cost, epsilon_flow)

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
            has_production, def_cost, prod_cost, epsilon_flow)
    end
    return data_flow
end

function crate_optim_problem(network, scenarios, proba, data_flow, invest_flow_cost)
    model = Model(CPLEX.Optimizer)

    # invest flow variables
    @variable(model, 0 <= invest_flow[i = 1:network.N, j = i:network.N; (i,j) 
        in network.edges || (j,i) in network.edges])

    # Flow variables
    @variable(model, flot[s = 1:scenarios, i = 1:network.N, j = i:network.N; (i,j) 
        in network.edges || (j,i) in network.edges])

    # Production variables
    @variable(model, 0 <= prod[s = 1:scenarios, n in 1:network.N; 
        data_flow[s].has_production[n] == 1])

    # Loss of load variables
    @variable(model, 0 <= def[s = 1:scenarios, i in 1:network.N])

    # Flow bounds
    invest_init = 5
    @constraint(model, flow_max_positive[s = 1:scenarios, (i,j) in network.edges], 
        flot[s,i,j] <= invest_init + invest_flow[i,j])
    @constraint(model, flow_max_negative[s = 1:scenarios, (i,j) in network.edges], 
        -(invest_flow[i,j] + invest_init) <= flot[s,i,j])

    @constraint(model, flow_conservation[s = 1:scenarios, n in 1:network.N], 
        sum(flot[s,n,i] for i in 1:network.N if (n,i) in network.edges) - 
        sum(flot[s,i,n] for i in 1:network.N if (i,n) in network.edges) + 
        (data_flow[s].has_production[n] == 1 ? prod[s,n] : 0) + 
        def[s,n] == data_flow[s].demands[n]
        )
    
    @constraint(model, invest_cost,
        sum( [invest_flow_cost[i,j] * invest_flow[i,j] for (i,j) in network.edges]) >= 0.0
        )

    @objective(model, Min, 
        sum(proba[s] * prod[s,n] * data_flow[s].prod_cost[n] for s in 1:scenarios, n in 1:network.N if data_flow[s].has_production[n] == 1) +
        sum(sum((proba[s] * data_flow[s].def_cost) .* def[s,:]) for s in 1:scenarios) +
        sum(proba[s] * data_flow[s].epsilon_flow * flot[s,i,j] for s in 1:scenarios, (i,j) in network.edges) + 
        sum(invest_flow_cost[i,j] * invest_flow[i,j] for (i,j) in network.edges)
        )
    
    return model
end

##################################################################################################
# Main
##################################################################################################
if PROGRAM_FILE == "problemGenerator.jl"
    N = 50
    graph_density = 10
    seed = 0
    print("Generating graph")
    @time network = create_network(N, graph_density, seed, plotGraph = false, drawGraph = true)


    seed >= 0 ? Random.seed!(seed) : nothing

    scenarios = 2
    print("Sample data     ")
    @time data_flow = sample_network_data(scenarios, network, 1:50, 1:10, 1000, 0.1)

    invest_flow_cost = Dict(zip(network.edges, rand(50:100, length(network.edges))))

    ##################################################################################################
    # Problem formulation
    ##################################################################################################
    print("Create model    ")
    proba = [0.99, 0.01]
    @time model = crate_optim_problem(network, scenarios, proba, data_flow, invest_flow_cost)

    set_silent(model)
    print("Solving model   ")
    @time optimize!(model)
    

    # Solution analysis
    println()
    println("Obj : ", objective_value(model))

    println("Invest solution cost = ", sum([ invest_flow_cost[i,j] * value(variable_by_name(model,"invest_flow[$i,$j]")) for (i,j) in network.edges]))    
    for (i,j) in network.edges
        if value(variable_by_name(model,"invest_flow[$i,$j]")) != 0.0
            println("    invest ", (i,j), " = ", value(variable_by_name(model,"invest_flow[$i,$j]")) )
        end
    end
    println("Defaillance")
    def_cnt = [sum([ value(variable_by_name(model,"def[$s,$n]")) > 0 ? 1 : 0 for n in 1:network.N]) for s in 1:scenarios]
    for s in 1:scenarios
        println("    Scenario ", s, " n_def = ", def_cnt[s])
    end
    println("    Def totale = ", sum( (proba .* def_cnt) ))

    println()
    println("Investment heuristic")
    @printf("%-15s%-15s%-15s%-10s\n", "Invest min", "Invest max", "Obj", "Def count")

    invest_min = 0.0
    invest_max = 1e5
    while invest_max - invest_min > 1
        if sum( (proba .* def_cnt) ) > 0
            global invest_min = sum([ invest_flow_cost[i,j] * value(variable_by_name(model,"invest_flow[$i,$j]")) for (i,j) in network.edges])
        else
            global invest_max = sum([ invest_flow_cost[i,j] * value(variable_by_name(model,"invest_flow[$i,$j]")) for (i,j) in network.edges])
        end

        #println("Invest range : ", invest_min, " - ", invest_max)
        set_normalized_rhs(constraint_by_name(model, "invest_cost"), (invest_max + invest_min) / 2)

        optimize!(model)
        global def_cnt = [sum([ value(variable_by_name(model,"def[$s,$n]")) > 0 ? 1 : 0 for n in 1:network.N]) for s in 1:scenarios]
        @printf("%-15.3f%-15.3f%-15.3f%-10.3f\n", invest_min, invest_max, objective_value(model), sum( (proba .* def_cnt) ))
        #println("    Obj : ", objective_value(model))
        #println("    Def totale = ", sum( (proba .* def_cnt) ))
    end
end