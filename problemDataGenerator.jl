#########################################################################################
# Import management
#########################################################################################
#=if "--install-packages" in ARGS
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
end=#

include("graph.jl")

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
    flow_cost::Dict{Edge, Float64}
    flow_init::Dict{Edge, Int}
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
    flow_init::Dict{Edge, Int}
    grad_prod::Float64
end


"""
function sample_network_data
    brief: Generates a random vector of NetworkFlowProblemData for each scenario
"""
function sample_network_data(scenarios, network, time_steps, demand_range, 
    prod_cost_range, unsupplied_cost, flow_cost_range, flow_init, grad_prod, has_production)

    data_flow = Vector{NetworkFlowProblemData}(undef, scenarios)
    for s in 1:scenarios

        demands = rand(demand_range, network.N, time_steps)

        # Sampling production nodes    
        prod_cost = Dict(zip(
            [i for i in 1:network.N if has_production[i] == 1],
            [rand(prod_cost_range, time_steps) for i in 1:sum(has_production)]
            ))

        flow_cost = Dict(zip(
            network.edges,
            rand(flow_cost_range, network.n_edges)
            ))
        
        #=flow_init = Dict(zip(
            network.edges,
            rand(0:flow_init_max, network.n_edges)
            ))
        =#

        data_flow[s] = NetworkFlowProblemData(network, demands, 
            has_production, unsupplied_cost, prod_cost, flow_cost, 
            flow_init, grad_prod)
    end
    return data_flow
end


function generate_probabilities(n_scenarios)
    proba = rand(1:100, n_scenarios)
    return (1/sum(proba)) .* proba
end



function investment_problem_data_generator(scenarios, network, time_steps, demand_range, 
    prod_cost_range, unsupplied_cost, flow_cost_range, flow_init_max, grad_prod, invest_flow_range, invest_prod_range)

    # Samplng probabilities
    proba = generate_probabilities(scenarios)

    # Sampling production nodes
    has_production = rand(0:1, network.N)
    if sum(has_production) == 0
        has_production[1] = 1
    end

    flow_init = Dict(zip(
        network.edges,
        rand(0:flow_init_max, network.n_edges)
        ))

    # Network flow data
    data_flow = sample_network_data(scenarios, network, time_steps, demand_range, 
    prod_cost_range, unsupplied_cost, flow_cost_range, flow_init, grad_prod, has_production)

    # Sampling investment costs
    invest_flow_cost = Dict(zip(network.edges, rand(invest_flow_range, network.n_edges)))
    
    prod_nodes = [i for i in 1:network.N if data_flow[1].has_production[i] == 1]
    invest_prod_cost = Dict(zip(prod_nodes, rand(invest_prod_range, length(prod_nodes) )))


    
    return StochasticProblemData(scenarios, data_flow, proba, network, has_production, time_steps, 
        invest_flow_cost, invest_prod_cost, flow_init, grad_prod)
end


function investment_problem_data_generator(options, network)
    investment_problem_data_generator(
        options.scenarios, network, options.time_steps, options.demand_range, 
        options.prod_cost_range, options.unsupplied_cost, options.flow_cost_range, 
        options.flow_init_max, options.grad_prod, options.invest_flow_range, 
        options.invest_prod_range)
end


function production_nodes(data)
    ( n for n in 1:data.network.N if data.has_production[n] == 1 )
end

function edge_dict()
    Dict(zip(data.network.edges,[0.0 for e in data.network.edges]))
end

function prod_dict()
    Dict(zip(
            [i for i in production_nodes(data)],
            [0.0 for i in 1:sum(data.has_production)])
        )
end