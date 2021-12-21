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
    demands::Matrix{Float64}
    has_production::Vector{Int}
    unsupplied_cost::Float64
    prod_cost::Dict{Int,Vector{Float64}}
    epsilon_flow::Float64
end


"""
struct StochasticProblemData
    brief: contains a vector of NetworkFlowProblemData associated with probabilities
"""
struct StochasticProblemData
    S::Int # Number of scenarios
    scenario::Vector{NetworkFlowProblemData}
    probability::Vector{Float64}
    network::Network
    T::Int # Time steps
    invest_flow_cost::Dict{Edge, Float64}
end


"""
function sample_network_data
    brief: Generates a random vector of NetworkFlowProblemData for each scenario
"""
function sample_network_data(scenarios, network, time_steps, demand_range, 
    prod_cost_range, unsupplied_cost, epsilon_flow)

    # Same production for each scenario
    has_production = rand(0:1, network.N)

    data_flow = Vector{NetworkFlowProblemData}(undef, scenarios)
    for s in 1:scenarios

        demands = rand(demand_range, network.N, time_steps)

        # Sampling production nodes    
        prod_cost = Dict(zip(
            [i for i in 1:N if has_production[i] == 1],
            [rand(prod_cost_range, time_steps) for i in 1:sum(has_production)]
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

function investment_problem_data_generator(scenarios, network, time_steps, demand_range, 
    prod_cost_range, unsupplied_cost, epsilon_flow, invest_cost_range)
    
    data_flow = sample_network_data(scenarios, network, time_steps, demand_range, 
    prod_cost_range, unsupplied_cost, epsilon_flow)

    proba = generate_probabilities(scenarios)

    invest_flow_cost = Dict(zip(network.edges, rand(invest_cost_range, network.n_edges)))

    return StochasticProblemData(scenarios, data_flow, proba, network, time_steps, invest_flow_cost)
end