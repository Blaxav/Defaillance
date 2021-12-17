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
function create_invest_optim_problem(network, scenarios, proba, data_flow, invest_flow_cost)
    model = Model(CPLEX.Optimizer)

    # invest flow variables
    @variable(model, 0 <= invest_flow[(i,j) in network.edges])

    # Flow variables
    @variable(model, flow[s = 1:scenarios, (i,j) in network.edges])

    # Production variables
    @variable(model, 0 <= prod[s = 1:scenarios, n in 1:network.N; 
        data_flow[s].has_production[n] == 1])

    # Loss of load variables
    @variable(model, 0 <= unsupplied[s = 1:scenarios, i in 1:network.N])

    # Flow bounds
    invest_init = 5
    @constraint(model, flow_max_positive[s = 1:scenarios, (i,j) in network.edges], 
        flow[s,(i,j)] <= invest_init + invest_flow[(i,j)])
    @constraint(model, flow_max_negative[s = 1:scenarios, (i,j) in network.edges], 
        -(invest_flow[(i,j)] + invest_init) <= flow[s,(i,j)])

    @constraint(model, flow_conservation[s = 1:scenarios, n in 1:network.N], 
        sum(flow[s,(n,i)] for i in 1:network.N if (n,i) in network.edges) - 
        sum(flow[s,(i,n)] for i in 1:network.N if (i,n) in network.edges) + 
        (data_flow[s].has_production[n] == 1 ? prod[s,n] : 0) + 
        unsupplied[s,n] == data_flow[s].demands[n]
        )
    
    # Constraint to lead heuristic
    @constraint(model, invest_cost,
        sum( [invest_flow_cost[(i,j)] * invest_flow[(i,j)] for (i,j) in network.edges]) >= 0.0
        )

    @objective(model, Min, 
        sum(proba[s] * prod[s,n] * data_flow[s].prod_cost[n] for s in 1:scenarios, n in 1:network.N if data_flow[s].has_production[n] == 1) +
        sum(sum((proba[s] * data_flow[s].unsupplied_cost) .* unsupplied[s,:]) for s in 1:scenarios) +
        sum(proba[s] * data_flow[s].epsilon_flow * flow[s,(i,j)] for s in 1:scenarios, (i,j) in network.edges) + 
        sum(invest_flow_cost[(i,j)] * invest_flow[(i,j)] for (i,j) in network.edges)
        )
    
    return StochasticProblem(model, invest_flow, [], flow, prod, unsupplied)
end

"""
function counting_unsupplied_scenario
    brief: Computes the number of nodes which loss of load for a given scenario after optimization
"""
function counting_unsupplied_scenario(stoch_prob, scenario)
    return sum([ value(stoch_prob.unsupplied[scenario,n]) > 0 ? 1 : 0 for n in 1:network.N])
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
    @constraint(stoch_prob.model, unsupplied_to_zero[s = 1:scenarios, n = 1:network.n],
        has_unsupplied[s,n] <= (1/epsilon_cnt)*stoch_prob.unsupplied[s,n] )
    @constraint(stoch_prob.model, unsupplied_to_one[s = 1:scenarios, n = 1:network.n],
        2*networkData.demands[s,n]*has_unsupplied[s,n] >= stoch_prob.unsupplied[s,n] * epsilon_cnt )

    @constraint(stoch_prob.model, unsupplied_cnt, sum(proba .* [sum(stoch_prob.has_unsupplied[:,n] for n = 1:network.N)]) <= n_unsupplied )
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