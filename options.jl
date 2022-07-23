#################################################
# Options struct and constructor
#################################################
mutable struct GlobalOptions
    N
    graph_density
    seed
    scenarios
    time_steps
    demand_range
    prod_cost_range
    unsupplied_cost
    flow_cost_range
    grad_prod
    invest_flow_range
    invest_prod_range
    flow_init_max

    unsupplied_tolerance
    max_unsupplied

    algorithm

    cut_aggregation
    step_size
    stab_center_tol
    init_mean_value_solution

    heuristic_frequency
    heuristic_strategy

    bilevel_mode
    primal_big_M
    dual_big_M

    time_limit
end

function define_default_options()
    GlobalOptions(
        1, 100, 0, 1, 1, 
        1:1, 1:1, 1000, 1:1,
        1.0, 1:1, 1:1, 0.0, 
        1e-3, 3, "benders", 
        "monocut", 0.5, 0.3, 1,
        "Opt", "Rand","SOS1",
        1e3, 1e3, 60
    )
end

function read_option_file(path_to_file)

    options = define_default_options()

    f = open(path_to_file, "r")
    for line in readlines(f)  
        parse_option_line(options, line)
    end
    close(f)

    return options
end

function parse_option_line(options, line)
    if line == ""
        return
    elseif line[1] == '#'
        return
    else
        keyword = rstrip(split(line, "=")[1])
        value = lstrip(rstrip(split(line, "=")[2]))
        
        # Instance options
        if keyword == "N"
            options.N = parse(Int64, value)
        elseif keyword == "graph_density"
            options.graph_density = parse(Float64, value)
        elseif keyword == "seed"
            options.seed = parse(Int64, value)
        elseif keyword == "scenarios"
            options.scenarios = parse(Int64, value)
        elseif keyword == "time_steps"
            options.time_steps = parse(Int64, value)

        elseif keyword == "demand_range"
            min_val = parse(Int64, split(value, ":")[1])
            max_val = parse(Int64, split(value, ":")[2])
            options.demand_range = min_val:max_val
        elseif keyword == "prod_cost_range"
            min_val = parse(Int64, split(value, ":")[1])
            max_val = parse(Int64, split(value, ":")[2])
            options.prod_cost_range = min_val:max_val
        elseif keyword == "unsupplied_cost"
            options.unsupplied_cost = parse(Float64, value)
        elseif keyword == "flow_cost_range"
            min_val = parse(Int64, split(value, ":")[1])
            max_val = parse(Int64, split(value, ":")[2])
            options.flow_cost_range = min_val:max_val
        elseif keyword == "grad_prod"
            options.grad_prod = parse(Float64, value)
        elseif keyword == "invest_flow_range"
            min_val = parse(Int64, split(value, ":")[1])
            max_val = parse(Int64, split(value, ":")[2])
            options.invest_flow_range = min_val:max_val
        elseif keyword == "invest_prod_range"
            min_val = parse(Int64, split(value, ":")[1])
            max_val = parse(Int64, split(value, ":")[2])
            options.invest_prod_range = min_val:max_val
        elseif keyword == "flow_init_max"
            options.flow_init_max = parse(Float64, value)
        
        elseif keyword == "unsupplied_tolerance"
            options.unsupplied_tolerance = parse(Float64, value)
        elseif keyword == "max_unsupplied"
            options.max_unsupplied = parse(Float64, value)

        # Algorithm options
        elseif keyword == "algorithm"
            options.algorithm = value
        
        elseif keyword == "cut_aggregation"
            options.cut_aggregation = value
        elseif keyword == "step_size"
            options.step_size = parse(Float64, value)
        elseif keyword == "stab_center_tol"
            options.stab_center_tol = parse(Float64, value)
        elseif keyword == "init_mean_value_solution"
            options.init_mean_value_solution = parse(Bool, value)

        elseif keyword == "heuristic_frequency"
            options.heuristic_frequency = value
        elseif keyword == "heuristic_strategy"
            options.heuristic_strategy = value

        elseif keyword == "bilevel_mode"
            options.bilevel_mode = value
        elseif keyword == "primal_big_M"
            options.primal_big_M = parse(Float64, value)
        elseif keyword == "dual_big_M"
            options.dual_big_M = parse(Float64, value)

        elseif keyword == "time_limit"
            options.time_limit = parse(Float64, value)
        else
            println("Unknown option ", keyword)
            exit()
        end
    end
end

function print_options(options)
    println("###################")
    println("# OPTIONS")
    println("###################")
    for n in fieldnames(typeof(options))
        @printf("%-25s", n)
        println(getfield(options,n))
    end
    println()
end



mutable struct Algorithm
    name
    
    # Bilevel parameters
    bilevel_mode # BilevelJuMP.SOS1Mode(), BilevelJuMP.IndicatorMode(), BilevelJuMP.FortunyAmatMcCarlMode()
    
    # Benders algorithms parameters
    cut_aggregation # monocut or multicut
    step_size # Step size in in-out stabilization
    stab_center_tol # Tolerance to accept a new better solution in (0, 1.0)
                    # computed as UB - tol*Gap
    init_mean_value_solution
    
    # Heuristic Parameters
    heuristic_frequency # All or Opt (check all invest solution or only optimal solution of Benders)
    heuristic_strategy # Rand or Min, Rand: Solution given by solver, Min: best solution with auxiliary prob solving

    time_limit
end

function create_algo(options)
    if options.algorithm == "benders"
        return create_benders(options.cut_aggregation, options.step_size, 
            options.stab_center_tol, options.init_mean_value_solution, 
            options.time_limit)
    elseif options.algorithm == "bilevel"
        return create_bilevel(options.bilevel_mode, options.time_limit; 
            primal_ub=options.primal_big_M, dual_ub=options.dual_big_M)
    elseif options.algorithm == "heuristic"
        return create_heuristic(options.cut_aggregation, options.step_size, 
            options.stab_center_tol, options.init_mean_value_solution, 
            options.heuristic_frequency, options.heuristic_strategy, 
            options.time_limit)
    else
        println("Unknown algorithm ", options.algorithm)
        exit()
    end
end

function create_bilevel(mode, time_limit; primal_ub=1e3, dual_ub=1e3)
    if mode == "SOS1"
        return Algorithm("bilevel", BilevelJuMP.SOS1Mode(), 
            any, any, any, any, any, any, time_limit)
    elseif mode == "Indicators"
        return Algorithm("bilevel", BilevelJuMP.IndicatorMode(), 
            any, any, any, any, any, any, time_limit)
    elseif mode == "Big-M"
        return Algorithm("bilevel", BilevelJuMP.FortunyAmatMcCarlMode(
            primal_big_M = primal_ub, dual_big_M=dual_ub), 
            any, any, any, any, any, any, time_limit)
    else
        println("Unknown bilevel mode ", mode)
        exit()
    end
end

function create_benders(cut_aggregation, step_size, stab_center_tol, 
    init_mean_value_solution, time_limit)
    
    return Algorithm("benders", any, cut_aggregation, step_size, 
        stab_center_tol, init_mean_value_solution, any, any, time_limit)
end

function create_heuristic(cut_aggregation, step_size, stab_center_tol, 
    init_mean_value_solution, frequency, strategy, time_limit)
    
    return Algorithm("heuristic", any, cut_aggregation, step_size,
        stab_center_tol, init_mean_value_solution, frequency, strategy, time_limit)
end