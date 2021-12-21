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
include("problemDataGenerator.jl")

struct StochasticProblem
    model::Model
    invest_flow
    invest_prod
    flow
    prod
    unsupplied
end

"""
function create_invest_optim_problem
    brief: Creates an stochastic investment optimization problem on a network
        It is possible to invest on every edge of the graph
"""
function create_invest_optim_problem(data)
    model = Model(CPLEX.Optimizer)

    # invest flow variables
    @variable(model, 0 <= invest_flow[e in data.network.edges])

    # Flow variables
    @variable(model, flow[s in 1:data.S, e in data.network.edges, t in 1:data.T])

    # Production variables
    @variable(model, 0 <= prod[s in 1:data.S, n in 1:data.network.N, t in 1:data.T; 
        data.scenario[s].has_production[n] == 1])

    # Loss of load variables
    @variable(model, 0 <= unsupplied[s in 1:data.S, i in 1:data.network.N, t in 1:data.T])


    # Flow bounds
    invest_init = 5
    @constraint(model, flow_max_positive[s in 1:data.S, e in data.network.edges, t in 1:data.T], 
        flow[s,e,t] <= invest_init + invest_flow[e])
    @constraint(model, flow_max_negative[s in 1:data.S, e in data.network.edges, t in 1:data.T], 
        -(invest_flow[e] + invest_init) <= flow[s,e,t])

    
    @constraint(model, flow_conservation[s in 1:data.S, n in 1:data.network.N, t in 1:data.T], 
        sum(flow[s,e,t] for e in data.network.edges if e.to == n) - 
        sum(flow[s,e,t] for e in data.network.edges if e.from == n) + 
        (data.scenario[s].has_production[n] == 1 ? prod[s,n,t] : 0) + 
        unsupplied[s,n,t] == data.scenario[s].demands[n,t]
        )
    
    # Constraint to lead heuristic
    @constraint(model, invest_cost,
        sum( [data.invest_flow_cost[e] * invest_flow[e] for e in data.network.edges]) >= 0.0
        )
    
    @objective(model, Min,
        sum( data.invest_flow_cost[e] * invest_flow[e] for e in data.network.edges ) +
        # Sum on scenarios
        sum( data.probability[s] *
            (   
            # Sum on time steps
            sum( 
                (
                # unsupplied costs
                sum( data.scenario[s].unsupplied_cost * unsupplied[s,n,t] 
                    for n in 1:data.network.N ) +
                # production costs
                sum( data.scenario[s].prod_cost[n][t] * prod[s,n,t] 
                    for n in 1:data.network.N 
                    if data.scenario[s].has_production[n] == 1) +
                # flow cost
                sum( data.scenario[s].epsilon_flow * flow[s,e,t] 
                    for e in data.network.edges )
                ) for t in 1:data.T
            )
            ) for s in 1:data.S
        )
    )
    
    return StochasticProblem(model, invest_flow, [], flow, prod, unsupplied)
end

"""
function counting_unsupplied_scenario
    brief: Computes the number of nodes which loss of load for a given scenario after optimization
"""
function counting_unsupplied_scenario(stoch_prob, scenario, epsilon, data)
    return sum( value(stoch_prob.unsupplied[scenario,n,t]) > epsilon ? 1 : 0 for n in 1:data.network.N, t in 1:data.T)
end

"""
function add_unsupplied_counter_constraint
    brief: Adds a binary variable on each node to say if it has unsupplied energy 
        and a constraints to limit the number of node with unsupplied energy to n_unsupplied
    args:
        stoch_prob : an instance of StochasticProblem
        n_unsupplied: maximum number of nodes with unsupplied enery in an optimal solution
        network: an instance of Network
        proba: vector of prabilities 
    returns: references to new variables and new constraint
"""
function add_unsupplied_counter_constraint(stoch_prob, n_unsupplied, network, scenarios, proba, epsilon_cnt, networkData)
    @variable(stoch_prob.model, has_unsupplied[s = 1:scenarios, i in 1:network.N], Bin)

    # Variables behaviour
    @constraint(stoch_prob.model, unsupplied_to_zero[s = 1:scenarios, n = 1:network.N],
        has_unsupplied[s,n] <= (1/epsilon_cnt)*stoch_prob.unsupplied[s,n] )
    @constraint(stoch_prob.model, unsupplied_to_one[s = 1:scenarios, n = 1:network.N],
        2*networkData[s].demands[n]*has_unsupplied[s,n] >= stoch_prob.unsupplied[s,n] - epsilon_cnt )

    @constraint(stoch_prob.model, unsupplied_cnt, sum(proba .* sum(has_unsupplied[:,n] for n = 1:network.N) ) <= n_unsupplied )
    return has_unsupplied, unsupplied_cnt
end

"""
function solve
    brief: solves a stochastic problem directly
    returns: time to solve
"""
function solve(stoch_prob, silent_mode)

    silent_mode == true ? set_silent(stoch_prob.model) : nothing
    timer = @elapsed optimize!(stoch_prob.model)
    return timer
end