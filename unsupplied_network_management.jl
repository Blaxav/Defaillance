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
println()
N = 20
graph_density = 10
seed = 1
print("Generating graph ")
@time network = create_network(N, graph_density, seed, plotGraph = false, drawGraph = true)


seed >= 0 ? Random.seed!(seed) : nothing

scenarios = 1
time_steps = 24
demand_range = 50:400
prod_cost_range = 200:800
unsupplied_cost = 1000
epsilon_flow = 0.1
grad_prod = 0.1
invest_cost_range = 1000:8000
invest_prod_range = 1000:8000
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
@time unsup = value.(stoch_prob.unsupplied)
@time unsupplied_cnt = [counting_unsupplied_scenario(stoch_prob, s, 0.0, data) for s in 1:data.S]
println("    unsupplied totale = ", sum( data.probability .* unsupplied_cnt ))


########################################
# Heuristic
########################################
println()
println()
println("#########################################################################################")
println("Heuristic ")
max_unsupplied = 3
@time investment_heuristic(stoch_prob, data, max_unsupplied, 1e-6, true, false)

########################################
# N hours counting constraint
########################################
println()
println("#########################################################################################")
print("Solving stochastic model with binary variables ")
epsilon_cnt = 0.0001
has_unsupplied, unsupplied_cnt_cstr = add_unsupplied_counter_constraint(stoch_prob, max_unsupplied, epsilon_cnt, data)

@time solve(stoch_prob, true)

println()
println("Obj : ", objective_value(stoch_prob.model))

println("Invest solution cost = ", investment_cost(stoch_prob, data) )
println("has_unsupplied")
has_unsupplied_cnt = [sum( value(has_unsupplied[s,n,t]) for n in 1:data.network.N, t in 1:data.T) for s in 1:scenarios]
println("    unsupplied totale = ", sum( (data.probability .* has_unsupplied_cnt) ))


########################################
# Bilevel problem
########################################
println()
println("#########################################################################################")
print("Create bilevel model     ")
@time bilev = create_bilevel_invest_problem(data, epsilon_cnt, max_unsupplied)

print("Solving bilevel model     ")
@time solve(bilev, true)

println("Obj bilevel : ", objective_value(bilev.model))
println("Invest solution cost = ", investment_cost(bilev, data))    
println("unsupplied")
unsupplied_cnt = [ sum([ value(bilev.has_unsupplied[s,n,t]) for n in 1:network.N, t in 1:data.T]) for s in 1:scenarios]
for s in 1:scenarios
    println("    Scenario ", s, " n_unsupplied = ", unsupplied_cnt[s])
end
println("    unsupplied totale = ", sum( (data.probability .* unsupplied_cnt) ))
println()

#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
exit()
#########################################################################################
# Reading Antares Study
#########################################################################################
println("#########################################################################################")
println("#########################################################################################")
println("#########################################################################################")
println()
path = "C:/Users/xavbl/these/instances/OSMOSE_20210528"
antaGraph, node_to_int = generate_graph_from_antares(path)
graph_plot_final = plot(antaGraph.positions, seriestype = :scatter, showaxis = false, grid = false, ticks = false, legend = false)
for e in antaGraph.edges
    plot!([antaGraph.positions[e.from][1], antaGraph.positions[e.to][1]], 
    [antaGraph.positions[e.from][2], antaGraph.positions[e.to][2]])
end

#display(graph_plot_final)
#println("Press a key to continue")
#readline()

initial_capacities = get_init_link_capacities_from_antares(path, antaGraph, node_to_int)
loads = get_load_chronicles_from_antares(path, antaGraph, node_to_int)

scenarios = 2
time_steps = 10
demand_range = 50:200
prod_cost_range = 50:200
unsupplied_cost = 10000
epsilon_flow = 0.1
grad_prod = 0.7
invest_cost_range = 10000:80000
invest_prod_range = 10000:80000
print("Sample data      ")
@time data = investment_problem_data_generator(scenarios, antaGraph, time_steps, demand_range, 
prod_cost_range, unsupplied_cost, epsilon_flow, grad_prod, invest_cost_range, invest_prod_range)

# Taking known demands after sampling
"""for s in 1:scenarios
    for (id, load) in loads
        #println(load[1][1+(s-1)*168:s*168])
        for t in 1:time_steps
            data.scenario[s].demands[id,t] = load[1][(s-1)*168 + t]
        end
    end
end"""

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
@time unsupplied_cnt = [counting_unsupplied_scenario(stoch_prob, s, 0.0, data) for s in 1:data.S]
@time println("    unsupplied totale = ", sum( data.probability .* unsupplied_cnt ))
@time unsup = value.(stoch_prob.unsupplied)