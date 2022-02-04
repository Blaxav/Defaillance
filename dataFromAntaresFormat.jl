#########################################################################################
# Import management
#########################################################################################
if "--install-packages" in ARGS
    import Pkg
    Pkg.add("Random")
    Pkg.add("Plots")
end

if "--compile" in ARGS
    include("compile_plotpackage.jl")
    println("Package Plot.jl compiled. Relaunch without '--compile' and with '--sysimage sys_plots.so'")
    exit
end

try
    using Random
    using Plots
catch
    println("ERROR : Some recquired packages not found.")
    println("        Run with option --install-packages to install them")
    exit()
end


"""
function generate_graph_from_antares
Brief: Build the graph given in an antares study
Args:
    path: path to an antares study
Returns:
    an instance of Network, node_to_int a dict linking nodes names to their int id
"""
function generate_graph_from_antares(path)
    path = path * "/input2/links/"

    nodes = []
    edges = []
    positions = []
    x = 0
    y = 0

    dirs = filter(x -> isdir(joinpath(path, x)), readdir(path))
    for dir in dirs
        push!(nodes, dir)
        # Getting coordinates
        open(path * "/../areas/" * dir * "/ui.ini") do f
            for line in readlines(f)
                if first(split(line, " ")) == "x"
                    global x = parse(Int64, last(split(line, " ")))
                elseif first(split(line, " ")) == "y"
                    global y = parse(Int64, last(split(line, " ")))
                end
            end
            push!(positions, (x,y))
        end
    end
    node_to_int = Dict(zip(nodes, 1:length(nodes)))
    
    for dir in dirs
        for lk in readdir(path * dir)
            if lk != "properties.ini"
                zone = lowercase(lk[1:length(lk) - 4])
                push!(edges, Edge(node_to_int[dir], node_to_int[zone]))
            end
        end
    end

    return Network(length(nodes), positions, edges, length(edges)), node_to_int
end


"""
function get_init_link_capacities_from_antares
Brief: Finds the link capacities in an antares study
Args:
    path: path to an antares study
    network: graph of the study
    node_to_int: a dict linking nodes names to their int id
Returns:
    Dict(Edge -> capa)
"""
function get_init_link_capacities_from_antares(path, network, node_to_int)
    
    path = path * "/input2/links/"
    dirs = filter(x -> isdir(joinpath(path, x)), readdir(path))

    capacities = Dict{Edge, Float64}()
    for dir in dirs
        # Finding all link files
        for lk in readdir(path * dir)
            if lk != "properties.ini"
                zone = lowercase(lk[1:length(lk) - 4])
                # We take the max(capa_in, capa_out) as the capacity
                open(path * dir * "/" * lk) do f
                    line = readline(f)
                    capa = maximum([parse(Float64, split(line, "\t")[i]) for i in 1:2])
                    capacities[Edge(node_to_int[dir], node_to_int[zone])] = capa
                end
            end
        end
    end

    return capacities
end

"""
function get_load_chronicles_from_antares
Brief: Finds all loads chronicles in a study antares
Args:
    path: path to an antares study
    network: graph of the study
    node_to_int: a dict linking nodes names to their int id
Returns:
    Dict(Int -> Vector{Vector{Float64}}) for each node id a vector of scenarios, 
    each containing a chronicle of load
"""
function get_load_chronicles_from_antares(path, network, node_to_int)
    path = path * "/input/load/series/"

    # For each id, we have a vector of loads chronicles (which are vectors)
    loads = Dict{Int64, Vector{Vector{Float64}}}()

    for (zone, id) in node_to_int
        loads[id] = []
        open(path * "/load_" * zone * ".txt") do f
            # Initialization
            line_init = readline(f)

            # Splitting when finding \t or spaces
            # \t zero or more time, spaces zero or more times
            # but at least one time one of them (the last +)
            for elt in split(line_init, r"[\t* *]+")
                push!(loads[id], [parse(Float64,elt)] )
            end

            for l in readlines(f)
                i = 1
                for elt in split(l, r"[\t* *]+")
                    push!(loads[id][i], parse(Float64,elt) )
                    i += 1
                end
            end
        end
    end

    return loads
end