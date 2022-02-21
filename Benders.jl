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
# Benders problem creation
#########################################################################################
struct BendersMasterProblem
    model::Model
    invest_flow
    invest_prod
    theta_sum # expectation of epigraph variables
    theta # epigraph variable for each scenario
end


struct BendersSubroblem
    model::Model
    invest_flow
    invest_prod
    flow
    prod
    unsupplied
    spilled
end


function create_master_benders_problem(data)
    model = Model(Gurobi.Optimizer)

    # invest variables
    @variable(model, 0 <= invest_flow[e in data.network.edges])
    @variable(model, 0 <= invest_prod[n in 1:data.network.N; 
        data.has_production[n] == 1])
    
    # Constraint to lead heuristic
    @constraint(model, invest_cost,
        sum( [data.invest_flow_cost[e] * invest_flow[e] for e in data.network.edges]) +
        sum( data.invest_prod_cost[n] * invest_prod[n] for n in 1:data.network.N 
            if data.has_production[n] == 1) >= 0.0
        )

    # Epigraph variables
    @variable(model, theta_sum)
    @variable(model, 0 <= theta[s in 1:data.S])

    @constraint(model, epigrah_cost, theta_sum == sum( data.probability[s] * theta[s] for s in 1:data.S) )
    
    @objective(model, Min,
        sum( data.invest_flow_cost[e] * invest_flow[e] for e in data.network.edges ) +
        sum( data.invest_prod_cost[n] * invest_prod[n] for n in 1:data.network.N 
            if data.has_production[n] == 1) +
        theta_sum
    )

    return BendersMasterProblem(model, invest_flow, invest_prod, theta_sum, theta)
end


function create_benders_subproblem(data, s)
    model = Model(Gurobi.Optimizer)

    # invest variables
    @variable(model, 0 <= invest_flow[e in data.network.edges])
    @variable(model, 0 <= invest_prod[n in 1:data.network.N; 
        data.has_production[n] == 1])

    # Flow variables
    @variable(model, flow[e in data.network.edges, t in 1:data.T])

    # Production variables
    @variable(model, 0 <= prod[n in 1:data.network.N, t in 1:data.T; 
        data.scenario[s].has_production[n] == 1])
    
    # Production constraints
    @constraint(model, prod_max[n in 1:data.network.N, t in 1:data.T; 
        data.has_production[n] == 1],
        prod[n,t] <= invest_prod[n])
    
    @constraint(model, grad_positive[n in 1:data.network.N, t in 1:data.T; 
        data.has_production[n] == 1],
        prod[n,t] <= (t > 1 ? prod[n,t-1] : prod[n,data.T]) + data.scenario[s].grad_prod*invest_prod[n] )
    @constraint(model, grad_negative[n in 1:data.network.N, t in 1:data.T; 
        data.has_production[n] == 1],
        prod[n,t] >= (t > 1 ? prod[n,t-1] : prod[n,data.T]) - data.scenario[s].grad_prod*invest_prod[n] )

    # Loss of load variables
    @variable(model, 0 <= unsupplied[i in 1:data.network.N, t in 1:data.T])
    @variable(model, 0 <= spilled[i in 1:data.network.N, t in 1:data.T])

    # Flow bounds
    invest_init = 5
    @constraint(model, flow_max_positive[e in data.network.edges, t in 1:data.T], 
        flow[e,t] <= invest_init + invest_flow[e])
    @constraint(model, flow_max_negative[e in data.network.edges, t in 1:data.T], 
        -(invest_flow[e] + invest_init) <= flow[e,t])

    
    @constraint(model, flow_conservation[n in 1:data.network.N, t in 1:data.T], 
        sum(flow[e,t] for e in data.network.edges if e.to == n) - 
        sum(flow[e,t] for e in data.network.edges if e.from == n) + 
        (data.scenario[s].has_production[n] == 1 ? prod[n,t] : 0) + 
        unsupplied[n,t] - spilled[n,t] == data.scenario[s].demands[n,t]
        )

    
    @objective(model, Min,
        # Sum on time steps
        sum(
            (
            # unsupplied costs
            sum( data.scenario[s].unsupplied_cost * unsupplied[n,t] 
                for n in 1:data.network.N ) +
            # production costs
            sum( data.scenario[s].prod_cost[n][t] * prod[n,t] 
                for n in 1:data.network.N 
                if data.scenario[s].has_production[n] == 1) +
            # flow cost
            sum( data.scenario[s].flow_cost[e] * flow[e,t] 
                for e in data.network.edges )
            ) for t in 1:data.T
        )
    )

    return BendersSubroblem(model, invest_flow, invest_prod, flow, prod, unsupplied, spilled)
end



#########################################################################################
# Benders resolution
#########################################################################################
function get_master_value(bendersMaster::BendersMasterProblem)
    return objective_value(bendersMaster.model)
end

function get_master_solution(bendersMaster::BendersMasterProblem)
    candidate_flow = value.(bendersMaster.invest_flow)
    candidate_prod = value.(bendersMaster.invest_prod)
    return candidate_flow, candidate_prod
end

function fix_first_stage_candidate(subproblems, candidate_flow, candidate_prod, data)
    for s in 1:data.S
        # Can not broadcast fix on edges as they are modeled by a dictionary
        # Seems that broadcasting needs flat data
        for e in data.network.edges
            fix(subproblems[s].invest_flow[e], candidate_flow[e]; force=true)   
        end
        for n in 1:network.N 
            if data.scenario[1].has_production[n] == 1
                fix.(subproblems[s].invest_prod[n], candidate_prod[n]; force=true)
            end
        end
        #fix.(subproblems[s].invest_prod, candidate_prod; force=true)
    end
end

function investment_cost(invest_flow, invest_prod, data)
    sum( data.invest_flow_cost[e] * invest_flow[e] for e in data.network.edges) +
    sum( data.invest_prod_cost[n] * invest_prod[n] for n in 1:data.network.N 
            if data.has_production[n] == 1)
end


function benders_sequential(master, subproblems, data)

    global best_UB = 1e12
    global LB = -1e6
    global best_invest_flow = Dict(zip(
        data.network.edges,
        [0.0 for e in data.network.edges]
        ))
    global separation_flow = Dict(zip(
            data.network.edges,
            [0.0 for e in data.network.edges]
            ))
    global best_invest_prod = Dict(zip(
        [i for i in 1:network.N if data.scenario[1].has_production[i] == 1],
        [0.0 for i in 1:sum(data.scenario[1].has_production)]
        ))
    global separation_prod = Dict(zip(
            [i for i in 1:network.N if data.scenario[1].has_production[i] == 1],
            [0.0 for i in 1:sum(data.scenario[1].has_production)]
            ))
    global iteration = 0
    global alpha = 0.5


    global grad_flow = Dict(zip(
            data.network.edges,
            [0.0 for e in data.network.edges]
            ))
    global grad_prod = Dict(zip(
        [i for i in 1:network.N if data.scenario[1].has_production[i] == 1],
        [0.0 for i in 1:sum(data.scenario[1].has_production)]
        ))

    
    while best_UB - LB > 1e-6
    
        global iteration += 1
    
        # 1. Solving master problem
        solve(master, true)
    
        # 2. Get master solution (investment candidates)
        candidate_flow, candidate_prod = get_master_solution(master)
    
        for e in data.network.edges
            global separation_flow[e] = alpha * candidate_flow[e] + (1-alpha) * best_invest_flow[e]
        end
        for n in 1:network.N 
            if data.scenario[1].has_production[n] == 1
                global separation_prod[n] = alpha * candidate_prod[n] + (1-alpha) * best_invest_prod[n]
            end
        end
        #println(separation_flow)
        #println(separation_prod)
    
        # 3. Fix investment candidates in subproblems
        fix_first_stage_candidate(subproblems, separation_flow, separation_prod, data)
    
        local UB = investment_cost(separation_flow, separation_prod, data)
        local rhs = 0.0
        #println(UB)
        # 4. Solve Subproblems and get subprgradients
        for s in 1:data.S
            solve(subproblems[s], true)
    
            UB += data.probability[s] * objective_value(subproblems[s].model)
            rhs += data.probability[s] * objective_value(subproblems[s].model)
            #println(UB)
            if s == 1
                global grad_flow = data.probability[s] * reduced_cost.(subproblems[s].invest_flow)
                global grad_prod = data.probability[s] * reduced_cost.(subproblems[s].invest_prod)
            else
                for e in data.network.edges
                    global grad_flow[e] += data.probability[s] * reduced_cost(subproblems[s].invest_flow[e])
                end
                for n in 1:data.network.N
                    if data.scenario[s].has_production[n] == 1
                        global grad_prod[n] += data.probability[s] * reduced_cost(subproblems[s].invest_prod[n])
                    end
                end
            end
        end

        # Build cut
        @constraint(benders_master.model, benders_master.theta_sum >= 
        rhs +
        sum(grad_flow[e] * (benders_master.invest_flow[e] - separation_flow[e]) for e in data.network.edges) +
        sum(grad_prod[n] * (benders_master.invest_prod[n] - separation_prod[n]) for n in 1:data.network.N if data.scenario[1].has_production[n] == 1)
        )
    
        if UB < best_UB
            global alpha = min(1.0, 1.2*alpha)
            global best_UB = UB
            global best_invest_prod = separation_prod
            global best_invest_flow = separation_flow
        else
            global alpha = max(0.1, 0.8*alpha)
        end
    
        global LB = objective_value(master.model)
        if iteration % 10 == 0 || best_UB - LB <= 1e-6
            @printf("%-10i%-20.6e%-20.6e%-15.3e%15.3f\n", iteration, LB, best_UB, best_UB-LB, alpha)
        end
    end
end

function benders_lazy_constraint_callback(cb_data)

    global iter_num
    global best_bound
    iter_num += 1
    #println("Iteration number = ", iter_num)

    # Get master problem solution
    invest_prod_current = callback_value.(Ref(cb_data), benders_master.invest_prod)
    invest_flow_current = Dict{Edge, Float64}()
    for e in data.network.edges
        invest_flow_current[e] = callback_value(cb_data, benders_master.invest_flow[e])
    end
    theta_current = callback_value.(Ref(cb_data), benders_master.theta)
    theta_total_current = callback_value(cb_data, benders_master.theta_sum)


    # Fix first stage
    fix_first_stage_candidate(subproblems, invest_flow_current, invest_prod_current, data)


    local UB = investment_cost(invest_flow_current, invest_prod_current, data)
    local rhs = 0.0

    # 4. Solve Subproblems and get subprgradients
    for s in 1:data.S
        solve(subproblems[s], true)

        UB += data.probability[s] * objective_value(subproblems[s].model)
        
        rhs += data.probability[s] * objective_value(subproblems[s].model)
        #println(UB)
        if s == 1
            global grad_flow = data.probability[s] * reduced_cost.(subproblems[s].invest_flow)
            global grad_prod = data.probability[s] * reduced_cost.(subproblems[s].invest_prod)
        else
            for e in data.network.edges
                global grad_flow[e] += data.probability[s] * reduced_cost(subproblems[s].invest_flow[e])
            end
            for n in 1:data.network.N
                if data.scenario[s].has_production[n] == 1
                    global grad_prod[n] += data.probability[s] * reduced_cost(subproblems[s].invest_prod[n])
                end
            end
        end
    end




    if UB < best_bound
        global best_bound = UB
        println("BOUND FOUND ! ", best_bound)
    #    println("test bound : ", MOI.get(benders_master.model, MOI.ObjectiveBound()))
    #    MOI.set(benders_master.model, MOI.ObjectiveBound(), best_bound)
    end

    if rhs > theta_total_current + 1e-3
        # Build cut
        new_optimality_cons = @build_constraint(
            benders_master.theta_sum >= 
            rhs +
            sum(grad_flow[e] * (benders_master.invest_flow[e] - invest_flow_current[e]) for e in data.network.edges) +
            sum(grad_prod[n] * (benders_master.invest_prod[n] - invest_prod_current[n]) for n in 1:data.network.N if data.scenario[1].has_production[n] == 1)
            )
        MOI.submit(
            benders_master.model,
            MOI.LazyConstraint(cb_data),
            new_optimality_cons,
        )

    else
        println("No cut here !")
    end
    
end


function benders_user_cut_callback(cb_data)
    global user_iter
    user_iter += 1
    #println("Iteration number = ", user_iter)

    # Get master problem solution
    #invest_prod_current = callback_value.(Ref(cb_data), benders_master.invest_prod)
    invest_prod_current = Dict{Int, Float64}()
    for n in 1:data.network.N
        if data.scenario[1].has_production[n] == 1
            invest_prod_current[n] = round(Int, callback_value(cb_data, benders_master.invest_prod[n]))
        end
    end
    #println("User : ", invest_prod_current)

    invest_flow_current = Dict{Edge, Float64}()
    for e in data.network.edges
        invest_flow_current[e] = callback_value(cb_data, benders_master.invest_flow[e])
    end
    theta_current = callback_value.(Ref(cb_data), benders_master.theta)
    theta_total_current = callback_value(cb_data, benders_master.theta_sum)

    # Fix first stage
    fix_first_stage_candidate(subproblems, invest_flow_current, invest_prod_current, data)

    #println("Current integer sol : ",   invest_prod_current)

    local UB = investment_cost(invest_flow_current, invest_prod_current, data)
    local rhs = 0.0

    # 4. Solve Subproblems and get subprgradients
    for s in 1:data.S
        solve(subproblems[s], true)

        UB += data.probability[s] * objective_value(subproblems[s].model)
        
        rhs += data.probability[s] * objective_value(subproblems[s].model)
        #println(UB)
        if s == 1
            global grad_flow = data.probability[s] * reduced_cost.(subproblems[s].invest_flow)
            global grad_prod = data.probability[s] * reduced_cost.(subproblems[s].invest_prod)
        else
            for e in data.network.edges
                global grad_flow[e] += data.probability[s] * reduced_cost(subproblems[s].invest_flow[e])
            end
            for n in 1:data.network.N
                if data.scenario[s].has_production[n] == 1
                    global grad_prod[n] += data.probability[s] * reduced_cost(subproblems[s].invest_prod[n])
                end
            end
        end
    end

    #println("User UB = ", UB)

    # Build cut
    if rhs > theta_total_current + 1e3
        new_optimality_cons = @build_constraint(
            benders_master.theta_sum >= 
            rhs +
            sum(grad_flow[e] * (benders_master.invest_flow[e] - invest_flow_current[e]) for e in data.network.edges) +
            sum(grad_prod[n] * (benders_master.invest_prod[n] - invest_prod_current[n]) for n in 1:data.network.N if data.scenario[1].has_production[n] == 1)
            )
        MOI.submit(
            benders_master.model,
            MOI.UserCut(cb_data),
            new_optimality_cons,
        )
    else
        println("No cut here !")
    end
    
end

function my_heuristic_callback(cb_data)

    #println("Submit heuristic solution")
    invest_prod_current = callback_value.(Ref(cb_data), benders_master.invest_prod)
    
    invest_heuristic = [ round(Int, v) for v in invest_prod_current ]
    #println(" Vector de l'heuristique : ", invest_heuristic)
    prod_vars = [benders_master.invest_prod[n] for n in 1:data.network.N if data.scenario[1].has_production[n] == 1]

    status = MOI.submit(
        benders_master.model, 
        MOI.HeuristicSolution(cb_data), 
        prod_vars, 
        invest_heuristic
    )
    #println("I submitted a heuristic solution, and the status was: ", status)
end

#########################################################################################
# User options
#########################################################################################
N = 3
graph_density = 70
seed = 0
time_graph = @elapsed network = create_network(N, graph_density, seed, plotGraph = false, drawGraph = true)


seed >= 0 ? Random.seed!(seed) : nothing

scenarios = 1
time_steps = 1
demand_range = 1:10
prod_cost_range = 10:40
unsupplied_cost = 1000
epsilon_flow = 1.0
grad_prod = 1.0
invest_cost_range = 500:1000
invest_prod_range = 500:1000
time_data = @elapsed data = investment_problem_data_generator(scenarios, network, time_steps, demand_range, 
prod_cost_range, unsupplied_cost, epsilon_flow, grad_prod, invest_cost_range, invest_prod_range)

# Evaluation du probleme stochastique
time_stoch_prob_creation = @elapsed stoch_prob = create_invest_optim_problem(data)

time_solve_stoch_prob = @elapsed solve(stoch_prob, false)
val_stoch_prob = objective_value(stoch_prob.model)
println("Problem value : ", val_stoch_prob)
println("Optimal solution")
println(value.(stoch_prob.invest_prod))
for e in data.network.edges
    println(e, "  ", value(stoch_prob.invest_flow[e]))
end
println()
println("*************************")

# Benders problems creation
time_master_prob_creation = @elapsed benders_master = create_master_benders_problem(data)

println()
println("Subproblems creation")
println()

subproblems =  Vector{BendersSubroblem}(undef, scenarios)
for s in 1:data.S
    subproblems[s] = create_benders_subproblem(data, s)
    #println(subproblems[s].model)
end

MOI.set(
    benders_master.model,
    MOI.LazyConstraintCallback(),
    benders_lazy_constraint_callback,
)

MOI.set(
    benders_master.model,
    MOI.UserCutCallback(),
    benders_user_cut_callback,
)

MOI.set(
    benders_master.model,
    MOI.HeuristicCallback(),
    my_heuristic_callback,
)


# Two phases resolution
# 1. LP relaxation
unset_integer.(benders_master.invest_prod)

benders_sequential(benders_master, subproblems, data)

println(value.(benders_master.invest_prod))


####################################################
##################### STOP HERE ####################
####################################################
exit()

# 2. MIP problem
set_integer.(benders_master.invest_prod)

iter_num = 0
user_iter = 0
heuristic_iter = 0
best_bound = 1e12
unset_silent(benders_master.model)
optimize!(benders_master.model)


#@time benders_sequential(benders_master, subproblems, data)

println("Lazy calls : ", iter_num)
println("User callback calls :", user_iter)
println("Heuristic calls : ", heuristic_iter)

