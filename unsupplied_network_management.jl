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
include("dataFromAntaresFormat.jl")

#########################################################################################
# User options
#########################################################################################
generate_graph_from_antares(".")
exit()

N = 30
graph_density = 10
seed = 0
print("Generating graph ")
@time network = create_network(N, graph_density, seed, plotGraph = false, drawGraph = true)


seed >= 0 ? Random.seed!(seed) : nothing

scenarios = 10
time_steps = 168
demand_range = 50:200
prod_cost_range = 50:200
unsupplied_cost = 10000
epsilon_flow = 0.1
grad_prod = 0.7
invest_cost_range = 10000:80000
invest_prod_range = 10000:80000
print("Sample data      ")
@time data = investment_problem_data_generator(scenarios, network, time_steps, demand_range, 
prod_cost_range, unsupplied_cost, epsilon_flow, grad_prod, invest_cost_range, invest_prod_range)

println()

########################################
# Stochastic problem
########################################
println()
println("#########################################################################################")
print("Create stochastic model     ")
@time stoch_prob = create_invest_optim_problem(data)


print("Solving model    ")
@time solve(stoch_prob, true)

println()
println("Obj : ", objective_value(stoch_prob.model))
println("Invest solution cost = ", investment_cost(stoch_prob, data))  
println("unsupplied")
unsupplied_cnt = [counting_unsupplied_scenario(stoch_prob, s, 0.0, data) for s in 1:data.S]
println("    unsupplied totale = ", sum( data.probability .* unsupplied_cnt ))

########################################
# Heuristic
########################################
println()
println("Investment heuristic")
@printf("%-15s%-15s%-15s%-15s%-10s\n", "Invest min", "Invest max", "Invest ctr","Obj", "unsupplied count")

invest_min = 0.0
invest_max = 5*objective_value(stoch_prob.model)
max_unsupplied = 0
rhs = investment_cost(stoch_prob, data)

while invest_max - invest_min > 1e-6*invest_max
    if sum( (data.probability .* unsupplied_cnt) ) > max_unsupplied
        global invest_min = rhs
    else
        global invest_max = rhs
    end

    global rhs = (invest_max + invest_min) / 2
    set_normalized_rhs(constraint_by_name(stoch_prob.model, "invest_cost"), rhs)

    optimize!(stoch_prob.model)
    global unsupplied_cnt = [ counting_unsupplied_scenario(stoch_prob, s, 0.0, data) for s in 1:data.S ]
    @printf("%-15.3f%-15.3f%-15.3f%-15.3f%-10.3f\n", invest_min, invest_max, rhs, 
        objective_value(stoch_prob.model), sum( data.probability .* unsupplied_cnt ) )

end

set_normalized_rhs(constraint_by_name(stoch_prob.model, "invest_cost"), 0)

########################################
# N hours counting constraint
########################################
println("Adding unsupplied constraint")
epsilon_cnt = 0.001
@time has_unsupplied, unsupplied_cnt = add_unsupplied_counter_constraint(stoch_prob, max_unsupplied, epsilon_cnt, data)
println()

print("Solving stochastic model with binary variables ")
@time solve(stoch_prob, true)
println()
println("Obj : ", objective_value(stoch_prob.model))

println("Invest solution cost = ", investment_cost(stoch_prob, data) )    
println("has_unsupplied")
has_unsupplied_cnt = [sum( value(has_unsupplied[s,n,t]) for n in 1:data.network.N, t in 1:data.T) for s in 1:scenarios]
println("    unsupplied totale = ", sum( (data.probability .* has_unsupplied_cnt) ))
println()

########################################
# Bilevel problem
########################################
println("#########################################################################################")
print("Create bilevel model     ")
@time bilev = create_bilevel_invest_problem(data, epsilon_cnt, max_unsupplied)

print("Solving bilevel model     ")
@time solve(bilev, false)

println("Obj bilevel : ", objective_value(bilev.model))
println("Invest solution cost = ", investment_cost(bilev, data))    
println("unsupplied")
unsupplied_cnt = [ sum([ value(bilev.has_unsupplied[s,n,t]) for n in 1:network.N, t in 1:data.T]) for s in 1:scenarios]
for s in 1:scenarios
    println("    Scenario ", s, " n_unsupplied = ", unsupplied_cnt[s])
end
println("    unsupplied totale = ", sum( (data.probability .* unsupplied_cnt) ))
println()