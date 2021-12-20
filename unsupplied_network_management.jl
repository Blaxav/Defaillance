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

N = 20
graph_density = 30
seed = 0
print("Generating graph ")
@time network = create_network(N, graph_density, seed, plotGraph = false, drawGraph = true)


seed >= 0 ? Random.seed!(seed) : nothing

scenarios = 10
proba = generate_probabilities(scenarios)
println("Proba = ", proba)
print("Sample data      ")
@time data_flow = sample_network_data(scenarios, network, 1:50, 1:10, 1000, 0.1)

invest_flow_cost = Dict(zip(network.edges, rand(50:100, length(network.edges))))

println()
########################################
# Bilevel problem
########################################
println("#########################################################################################")
print("Create bilevel model     ")
epsilon_cnt = 0.01
n_unsupplied = 0.0
@time bilev = create_bilevel_invest_problem(network, scenarios, proba, data_flow, invest_flow_cost, epsilon_cnt, n_unsupplied)

print("Solving bilevel model     ")
@time solve(bilev, true)

println("Obj bilevel : ", objective_value(bilev.model))
println("Invest solution cost = ", sum([ invest_flow_cost[i,j] * value(bilev.invest_flow[(i,j)]) for (i,j) in network.edges]))    
println("unsupplied")
unsupplied_cnt = [ sum([ value(bilev.unsupplied[s,n]) > epsilon_cnt + 1e-6 ? 1 : 0 for n in 1:network.N]) for s in 1:scenarios]
for s in 1:scenarios
    println("    Scenario ", s, " n_unsupplied = ", unsupplied_cnt[s])
end
println("    unsupplied totale = ", sum( (proba .* unsupplied_cnt) ))
println()

########################################
# Stochastic problem
########################################
println()
println("#########################################################################################")
print("Create stochastic model     ")
@time stoch_prob = create_invest_optim_problem(network, scenarios, proba, data_flow, invest_flow_cost)


print("Solving model    ")
@time solve(stoch_prob, true)

println()
println("Obj : ", objective_value(stoch_prob.model))
println("Invest solution cost = ", sum([ invest_flow_cost[i,j] * value(stoch_prob.invest_flow[(i,j)]) for (i,j) in network.edges]))  
println("unsupplied")
unsupplied_cnt = [counting_unsupplied_scenario(stoch_prob, s) for s in 1:scenarios]
println("    unsupplied totale = ", sum( (proba .* unsupplied_cnt) ))


########################################
# Heuristic
########################################
println()
println("Investment heuristic")
@printf("%-15s%-15s%-15s%-15s%-10s\n", "Invest min", "Invest max", "Invest ctr","Obj", "unsupplied count")

invest_min = 0.0
invest_max = 1e6
max_unsupplied = 0

rhs = 0

while invest_max - invest_min > 1e-3
    if sum( (proba .* unsupplied_cnt) ) > max_unsupplied
        global invest_min = rhs
    else
        global invest_max = rhs
    end

    global rhs = (invest_max + invest_min) / 2
    set_normalized_rhs(constraint_by_name(stoch_prob.model, "invest_cost"), rhs)

    optimize!(stoch_prob.model)
    global unsupplied_cnt = [sum([ value(variable_by_name(stoch_prob.model,"unsupplied[$s,$n]")) > 0 ? 1 : 0 for n in 1:network.N]) for s in 1:scenarios]
    @printf("%-15.3f%-15.3f%-15.3f%-15.3f%-10.3f\n", invest_min, invest_max, rhs, objective_value(stoch_prob.model), sum( (proba .* unsupplied_cnt) ))
end

set_normalized_rhs(constraint_by_name(stoch_prob.model, "invest_cost"), 0)

########################################
# N hours counting constraint
########################################
has_unsupplied, unsupplied_cnt = add_unsupplied_counter_constraint(stoch_prob, max_unsupplied, network, scenarios, proba, 0.01, data_flow)
println()
println("Adding unsupplied constraint")
print("Solving stochastic model with binary variables ")
@time solve(stoch_prob, true)
println()
println("Obj : ", objective_value(stoch_prob.model))

println("Invest solution cost = ", sum([ invest_flow_cost[i,j] * value(stoch_prob.invest_flow[(i,j)]) for (i,j) in network.edges]))    
println("has_unsupplied")
has_unsupplied_cnt = [sum([ value(has_unsupplied[s,n]) for n in 1:network.N]) for s in 1:scenarios]
println("    unsupplied totale = ", sum( (proba .* has_unsupplied_cnt) ))