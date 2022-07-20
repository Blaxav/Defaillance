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
    epsilon_flow
    grad_prod
    invest_cost_range
    invest_prod_range
    flow_init_max
    algorithm
end

function define_default_options()
    GlobalOptions(
        1, 100, 0, 1, 1, 
        1:1, 1:1, 1000, 0.1,
        1.0, 1:1, 1:1, 0.0, "stochastic"
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
        elseif keyword == "epsilon_flow"
            options.epsilon_flow = parse(Float64, value)
        elseif keyword == "grad_prod"
            options.grad_prod = parse(Float64, value)
        elseif keyword == "invest_cost_range"
            min_val = parse(Int64, split(value, ":")[1])
            max_val = parse(Int64, split(value, ":")[2])
            options.invest_cost_range = min_val:max_val
        elseif keyword == "invest_prod_range"
            min_val = parse(Int64, split(value, ":")[1])
            max_val = parse(Int64, split(value, ":")[2])
            options.invest_prod_range = min_val:max_val
        elseif keyword == "flow_init_max"
            options.flow_init_max = parse(Float64, value)
        elseif keyword == "algorithm"
            options.algorithm = value
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