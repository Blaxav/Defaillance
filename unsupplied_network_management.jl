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
include("bilevelProblemGenerator.jl")
include("stochasticProblemGenerator.jl")

#########################################################################################
# User options
#########################################################################################

N = 50
graph_density = 10
seed = 0
print("Generating graph ")
@time network = create_network(N, graph_density, seed, plotGraph = false, drawGraph = true)


seed >= 0 ? Random.seed!(seed) : nothing

scenarios = 10
proba = (1/scenarios) .* ones(scenarios)
print(proba)
print("Sample data      ")
@time data_flow = sample_network_data(scenarios, network, 1:50, 1:10, 1000, 0.1)

invest_flow_cost = Dict(zip(network.edges, rand(50:100, length(network.edges))))

########################################
# Bilevel problem
########################################
print("Create bilevel model     ")
@time bilev = create_bilevel_invest_problem(network, scenarios, proba, data_flow, invest_flow_cost)

solve(bilev, false)

println("Obj bilevel : ", objective_value(bilev.model))
println("Invest solution cost = ", sum([ invest_flow_cost[i,j] * value(bilev.invest_flow[(i,j)]) for (i,j) in network.edges]))    
for (i,j) in network.edges
    if value(bilev.invest_flow[(i,j)]) != 0.0
        println("    invest ", (i,j), " = ", value(bilev.invest_flow[(i,j)]) )
    end
end
println("unsupplied")
unsupplied_cnt = [ sum([ value(bilev.unsupplied[s,n]) > 0 ? 1 : 0 for n in 1:network.N]) for s in 1:scenarios]
for s in 1:scenarios
    println("    Scenario ", s, " n_unsupplied = ", unsupplied_cnt[s])
end
println("    unsupplied totale = ", sum( (proba .* unsupplied_cnt) ))
println()


########################################
# Stochastic problem
########################################
println()
print("Create model     ")
@time stoch_prob = create_invest_optim_problem(network, scenarios, proba, data_flow, invest_flow_cost)


print("Solving model    ")
solve(stoch_prob, false)

println()
println("Obj : ", objective_value(stoch_prob.model))

println("Invest solution cost = ", sum([ invest_flow_cost[i,j] * value(stoch_prob.invest_flow[(i,j)]) for (i,j) in network.edges]))    
for (i,j) in network.edges
    if value(stoch_prob.invest_flow[(i,j)]) != 0.0
        println("    invest ", (i,j), " = ", value(stoch_prob.invest_flow[(i,j)]) )
    end
end
println("unsupplied")
unsupplied_cnt = [counting_unsupplied_scenario(stoch_prob, s) for s in 1:scenarios]
for s in 1:scenarios
    println("    Scenario ", s, " n_unsupplied = ", unsupplied_cnt[s])
end
println("    unsupplied totale = ", sum( (proba .* unsupplied_cnt) ))