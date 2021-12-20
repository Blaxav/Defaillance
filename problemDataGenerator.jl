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
struct StochasticProblemData
    brief: contains a vector of NetworkFlowProblemData associated with probabilities
"""
struct StochasticProblemData
    scenario::Vector{NetworkFlowProblemData}
    probability::Vector{Float64}
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

function generate_probabilities(n_scenarios)
    proba = rand(1:100, n_scenarios)
    return (1/sum(proba)) .* proba
end