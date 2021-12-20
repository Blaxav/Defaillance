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
    has_unsupplied
end

"""
function create_bilevel_invest_problem
    brief: Creates a stochastic bilevel investment optimization problem on a network
        It is possible to invest on every edge of the graph
"""
function create_bilevel_invest_problem(network, scenarios, proba, data_flow, invest_flow_cost,
    epsilon_cnt, n_unsupplied)

    model = BilevelModel(CPLEX.Optimizer, mode = BilevelJuMP.SOS1Mode())

    ## Upper level variables : invest flow variables
    @variable(Upper(model), 0 <= invest_flow[(i,j) in network.edges])

    ## Lower level variables
    # Flow variables
    @variable(Lower(model), flow[s = 1:scenarios, (i,j) in network.edges])
    # Production variables
    @variable(Lower(model), 0 <= prod[s = 1:scenarios, n in 1:network.N; 
        data_flow[s].has_production[n] == 1])
    # Loss of load variables
    @variable(Lower(model), 0 <= unsupplied[s = 1:scenarios, i in 1:network.N])

    # Counting unsupplied
    @variable(Upper(model), has_unsupplied[s = 1:scenarios, i in 1:network.N], Bin)


    ## Lower level constraints
    ## ----------------
    # Flow bounds
    invest_init = 5
    @constraint(Lower(model), flow_max_positive[s = 1:scenarios, (i,j) in network.edges], 
        flow[s,(i,j)] <= invest_init + invest_flow[(i,j)])
    @constraint(Lower(model), flow_max_negative[s = 1:scenarios, (i,j) in network.edges], 
        -(invest_init + invest_flow[(i,j)]) <= flow[s,(i,j)])
    # Flow conservation
    @constraint(Lower(model), flow_conservation[s = 1:scenarios, n in 1:network.N], 
        sum(flow[s,(n,i)] for i in 1:network.N if (n,i) in network.edges) - 
        sum(flow[s,(i,n)] for i in 1:network.N if (i,n) in network.edges) + 
        (data_flow[s].has_production[n] == 1 ? prod[s,n] : 0) + 
        unsupplied[s,n] == data_flow[s].demands[n]
        )
    
    # Counting variables behaviour
    @constraint(Upper(model), unsupplied_to_zero[s = 1:scenarios, n = 1:network.N],
        has_unsupplied[s,n] <= (1/epsilon_cnt)*unsupplied[s,n] )
    @constraint(Upper(model), unsupplied_to_one[s = 1:scenarios, n = 1:network.N],
        2*data_flow[s].demands[n]*has_unsupplied[s,n] >= unsupplied[s,n] - epsilon_cnt )

    @constraint(Upper(model), unsupplied_cnt, sum(proba .* sum(has_unsupplied[:,n] for n = 1:network.N) ) <= n_unsupplied )
    
    # Lower level obj : Objective without first stage cost to dualize 
    @objective(Lower(model), Min, 
        sum(proba[s] * prod[s,n] * data_flow[s].prod_cost[n] for s in 1:scenarios, n in 1:network.N if data_flow[s].has_production[n] == 1) +
        sum(sum((proba[s] * data_flow[s].unsupplied_cost) .* unsupplied[s,:]) for s in 1:scenarios) +
        sum(proba[s] * data_flow[s].epsilon_flow * flow[s,(i,j)] for s in 1:scenarios, (i,j) in network.edges)
        )
    
    # Upper level obj : Objective with first stage cost to dualize 
    @objective(Upper(model), Min, 
    sum(proba[s] * prod[s,n] * data_flow[s].prod_cost[n] for s in 1:scenarios, n in 1:network.N if data_flow[s].has_production[n] == 1) +
    sum(sum((proba[s] * data_flow[s].unsupplied_cost) .* unsupplied[s,:]) for s in 1:scenarios) +
    sum(proba[s] * data_flow[s].epsilon_flow * flow[s,(i,j)] for s in 1:scenarios, (i,j) in network.edges) +
    sum(invest_flow_cost[(i,j)] * invest_flow[(i,j)] for (i,j) in network.edges)
    )

    return BilevProblem(model,invest_flow,[], flow, prod, unsupplied, has_unsupplied)
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