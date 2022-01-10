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
    grad_prod::Float64
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
    has_production::Vector{Int}
    T::Int # Time steps
    invest_flow_cost::Dict{Edge, Float64}
    invest_prod_cost::Dict{Int, Float64}
end


"""
function sample_network_data
    brief: Generates a random vector of NetworkFlowProblemData for each scenario
"""
function sample_network_data(scenarios, network, time_steps, demand_range, 
    prod_cost_range, unsupplied_cost, epsilon_flow, grad_prod, has_production)

    data_flow = Vector{NetworkFlowProblemData}(undef, scenarios)
    for s in 1:scenarios

        demands = rand(demand_range, network.N, time_steps)

        # Sampling production nodes    
        prod_cost = Dict(zip(
            [i for i in 1:N if has_production[i] == 1],
            [rand(prod_cost_range, time_steps) for i in 1:sum(has_production)]
            ))

        data_flow[s] = NetworkFlowProblemData(network, demands, 
            has_production, unsupplied_cost, prod_cost, epsilon_flow, grad_prod)
    end
    return data_flow
end

function generate_probabilities(n_scenarios)
    proba = rand(1:100, n_scenarios)
    return (1/sum(proba)) .* proba
end

function investment_problem_data_generator(scenarios, network, time_steps, demand_range, 
    prod_cost_range, unsupplied_cost, epsilon_flow, grad_prod, invest_cost_range, invest_prod_range)

    # Samplng probabilities
    proba = generate_probabilities(scenarios)

    # Sampling production nodes
    has_production = rand(0:1, network.N)

    # Network flow data
    data_flow = sample_network_data(scenarios, network, time_steps, demand_range, 
    prod_cost_range, unsupplied_cost, epsilon_flow, grad_prod, has_production)

    # Sampling investment costs
    invest_flow_cost = Dict(zip(network.edges, rand(invest_cost_range, network.n_edges)))
    
    prod_nodes = [i for i in 1:network.N if data_flow[1].has_production[i] == 1]
    invest_prod_cost = Dict(zip(prod_nodes, rand(invest_prod_range, length(prod_nodes) )))
    
    return StochasticProblemData(scenarios, data_flow, proba, network, has_production, time_steps, 
        invest_flow_cost, invest_prod_cost)
end