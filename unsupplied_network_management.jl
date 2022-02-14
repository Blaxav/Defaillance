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
N = 15
graph_density = 30
#seed = 1

for seed in 1:50
    println()
    println("Seed ", seed)
    time_graph = @elapsed network = create_network(N, graph_density, seed, plotGraph = false, drawGraph = true)


    seed >= 0 ? Random.seed!(seed) : nothing

    scenarios = 1
    time_steps = 5
    demand_range = 50:400
    prod_cost_range = 100:200
    unsupplied_cost = 120
    epsilon_flow = 0.1
    grad_prod = 0.01
    invest_cost_range = 300:1200
    invest_prod_range = 300:1200
    time_data = @elapsed data = investment_problem_data_generator(scenarios, network, time_steps, demand_range, 
    prod_cost_range, unsupplied_cost, epsilon_flow, grad_prod, invest_cost_range, invest_prod_range)

    ########################################
    # Stochastic problem
    ########################################
    time_stoch_prob_creation = @elapsed stoch_prob = create_invest_optim_problem(data)

    time_solve_stoch_prob = @elapsed solve(stoch_prob, true)

    val_stoch_prob = objective_value(stoch_prob.model)
    invest_stoch_prob = investment_cost(stoch_prob, data)

    unsupplied_cnt = [counting_unsupplied_scenario(stoch_prob, s, 0.0, data) for s in 1:data.S]
    unsup_stoch_prob = sum( data.probability .* unsupplied_cnt )

    println()
    @printf("%-20s%-20s%-15s%-15s\n", "Invest cost", "Total cost", "Unsupplied", "Time")
    @printf("%-20.4e%-20.4e%-15.2f%-15.3f\n", invest_stoch_prob, val_stoch_prob, unsup_stoch_prob, time_solve_stoch_prob)

    ########################################
    # Heuristic
    ########################################
    max_unsupplied = 3
    time_heuristic = @elapsed invest_heuristic, val_heuristic = investment_heuristic(stoch_prob, data, max_unsupplied, 1e-6, true, false)

    @printf("%-20.4e%-20.4e%-15.2f%-15.3f\n", invest_heuristic, val_heuristic, max_unsupplied, time_heuristic)

    ########################################
    # N hours counting constraint
    ########################################
    epsilon_cnt = 0.0001
    has_unsupplied, unsupplied_cnt_cstr = add_unsupplied_counter_constraint(stoch_prob, max_unsupplied, epsilon_cnt, data)

    time_mip = @elapsed solve(stoch_prob, true)

    val_mip =  objective_value(stoch_prob.model)
    invest_mip = investment_cost(stoch_prob, data)

    has_unsupplied_cnt = [sum( value(has_unsupplied[s,n,t]) for n in 1:data.network.N, t in 1:data.T) for s in 1:scenarios]
    unsup_mip = sum( (data.probability .* has_unsupplied_cnt) )

    @printf("%-20.4e%-20.4e%-15.2f%-15.3f\n", invest_mip, val_mip, unsup_mip, time_mip)

    ########################################
    # Bilevel problem
    ########################################
    time_bilev_prob_creation = @elapsed bilev = create_bilevel_invest_problem(data, epsilon_cnt, max_unsupplied)
    time_bilev = @elapsed solve(bilev, true)

    val_bilev = objective_value(bilev.model)
    invest_bilev = investment_cost(bilev, data)

    unsupplied_cnt = [ sum([ value(bilev.has_unsupplied[s,n,t]) for n in 1:network.N, t in 1:data.T]) for s in 1:scenarios]
    unsup_bilev = sum( (data.probability .* unsupplied_cnt) )

    @printf("%-20.4e%-20.4e%-15.2f%-15.3f\n", invest_bilev, val_bilev, unsup_bilev, time_bilev)

    @printf("Bilevel-MIP Gap = %-.6e\n" , val_bilev - val_mip)
    if val_bilev - val_mip > 1e-2 * val_mip
        println("####################################")
        println("Instance found for seed = ", seed)
        println("####################################")
    end
end






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