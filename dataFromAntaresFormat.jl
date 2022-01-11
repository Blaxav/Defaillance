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

    return Network(length(nodes), positions, edges, length(edges))
end