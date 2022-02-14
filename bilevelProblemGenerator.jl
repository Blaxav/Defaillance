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

struct BilevProblem
    model::BilevelModel
    invest_flow
    invest_prod
    flow
    prod
    unsupplied
    spilled
    has_unsupplied
end

"""
function create_bilevel_invest_problem
    brief: Creates a stochastic bilevel investment optimization problem on a network
        It is possible to invest on every edge of the graph
"""
function create_bilevel_invest_problem(data, epsilon_cnt, max_unsupplied)

    model = BilevelModel(Gurobi.Optimizer, mode = BilevelJuMP.SOS1Mode())

    # invest variables
    @variable(Upper(model), 0 <= invest_flow[e in data.network.edges])
    @variable(Upper(model), 0 <= invest_prod[n in 1:data.network.N; 
        data.has_production[n] == 1])

    # Flow variables
    @variable(Lower(model), flow[s in 1:data.S, e in data.network.edges, t in 1:data.T])

    # Production variables
    @variable(Lower(model), 0 <= prod[s in 1:data.S, n in 1:data.network.N, t in 1:data.T; 
        data.scenario[s].has_production[n] == 1])
    
    @constraint(Lower(model), prod_max[s in 1:data.S, n in 1:data.network.N, t in 1:data.T; 
        data.scenario[s].has_production[n] == 1],
        prod[s,n,t] <= invest_prod[n])
    
    @constraint(Lower(model), grad_positive[s in 1:data.S, n in 1:data.network.N, t in 1:data.T; 
        data.has_production[n] == 1],
        prod[s,n,t] <= (t > 1 ? prod[s,n,t-1] : prod[s,n,data.T]) + data.scenario[s].grad_prod*invest_prod[n] )
    @constraint(Lower(model), grad_negative[s in 1:data.S, n in 1:data.network.N, t in 1:data.T; 
        data.has_production[n] == 1],
        prod[s,n,t] >= (t > 1 ? prod[s,n,t-1] : prod[s,n,data.T]) - data.scenario[s].grad_prod*invest_prod[n] )

    # Loss of load variables
    @variable(Lower(model), 0 <= unsupplied[s in 1:data.S, i in 1:data.network.N, t in 1:data.T])
    @variable(Lower(model), 0 <= spilled[s in 1:data.S, i in 1:data.network.N, t in 1:data.T])

    # Flow bounds
    invest_init = 5
    @constraint(Lower(model), flow_max_positive[s in 1:data.S, e in data.network.edges, t in 1:data.T], 
        flow[s,e,t] <= invest_init + invest_flow[e])
    @constraint(Lower(model), flow_max_negative[s in 1:data.S, e in data.network.edges, t in 1:data.T], 
        -(invest_flow[e] + invest_init) <= flow[s,e,t])

    
    @constraint(Lower(model), flow_conservation[s in 1:data.S, n in 1:data.network.N, t in 1:data.T], 
        sum(flow[s,e,t] for e in data.network.edges if e.to == n) - 
        sum(flow[s,e,t] for e in data.network.edges if e.from == n) + 
        (data.scenario[s].has_production[n] == 1 ? prod[s,n,t] : 0) + 
        unsupplied[s,n,t] - spilled[s,n,t] == data.scenario[s].demands[n,t]
        )
    
    @variable(Upper(model), has_unsupplied[s in 1:data.S, i in 1:data.network.N, t in 1:data.T], Bin)

    # Variables behaviour
    @constraint(Upper(model), unsupplied_to_zero[s in 1:data.S, n in 1:data.network.N, t in 1:data.T],
        has_unsupplied[s,n,t] <= (1/epsilon_cnt)*unsupplied[s,n,t] )
    @constraint(Upper(model), unsupplied_to_one[s in 1:data.S, n in 1:data.network.N, t in 1:data.T],
        2*data.scenario[s].demands[n,t]*has_unsupplied[s,n,t] >= unsupplied[s,n,t] - epsilon_cnt )
    
    @constraint(Upper(model), unsupplied_cnt, sum(data.probability .* sum( sum(has_unsupplied[:,n,t] for n = 1:data.network.N) for t in 1:data.T )) <= max_unsupplied )
    
    @objective(Lower(model), Min,
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

    @objective(Upper(model), Min,
        sum( data.invest_flow_cost[e] * invest_flow[e] for e in data.network.edges ) +
        sum( data.invest_prod_cost[n] * invest_prod[n] for n in 1:data.network.N 
            if data.has_production[n] == 1) +
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

    return BilevProblem(model,invest_flow, invest_prod, flow, prod, unsupplied, spilled, has_unsupplied)
end


function count_bilev_unsupplied_nodes(bilev_problem, scenarios, network)
    return [ sum([ value(bilev_problem.unsupplied[s,n]) > 0 ? 1 : 0 for n in 1:network.N]) for s in 1:scenarios]
end

"""
function solve
    brief: solves a bilevel problem with SOS1 method of MIP solvers
    returns: time to solve
"""
function solve(bilev_prob, silent_mode)

    silent_mode == true ? set_silent(bilev_prob.model) : nothing
    timer = @elapsed optimize!(bilev_prob.model)
    return timer
end

"""
function investment_cost
    brief: computes the investment cost of the inner solution of stoch_prob
"""
function investment_cost(bilev_prob, data)
    sum( data.invest_flow_cost[e] * value(bilev_prob.invest_flow[e]) for e in data.network.edges)
end