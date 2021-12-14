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

#########################################################################################
# Network type
#########################################################################################
struct Network
    N::Int # Both the number of nodes and their indices
    positions::Vector{Tuple{Int, Int}} # Positions of the nodes on a map
    edges::Vector{Tuple{Int, Int}} # list of undirected edges
end

#########################################################################################
#########################################################################################
function help()
    println("Possible options : ")
    
    opt_dict = Dict(
        "-h, ?"                 => "display help",
        "--install-packages"    => "install recquired packages",
        "--compile"             => "Use PackageCompiler to compile Plot.jl library to largely improve perfs with the use of 'julia --sysimage sys_plots.so'",
        "-N"                    => "Number of node of the graph",
        "--seed"                => "Random seed, a negative seed leads to pure random",
        "-d"                    => "Graph density to reach",
        "-p"                    => "Proportion of node on the map : caracterized the size of the map"
    )
   
    for opt in keys(opt_dict)
        println("     ", opt, " : ", opt_dict[opt])
    end
end

function dist(p1, p2)
    return sqrt( (p1[1] - p2[1])^2 + (p1[2] - p2[2])^2)
end

"""
Function sample_positions
Biref: Samples N positions on a grid. Return a vector of tuple (Int, Int)
Args:
    N : number of points to sample
    gridRatio : ratio between X length and Y length
    proportion : Defines the surface of the grid as N/proportion
Returns:
    Vector{Tuple{Int, Int}} of integer positions
"""
function sample_positions(N; gridRatio=1.5, proportion=0.005)
    Xmax = sqrt(gridRatio*N/proportion)
    Ymax = sqrt(N/(gridRatio*proportion))

    positions = Vector{Tuple{Int, Int}}()
    distances = Array{Float64}(undef, N, N)
    for i in 1:N
        # Find unused coordinates for new point
        a = rand(1:Xmax)
        b = rand(1:Ymax)
        while (a,b) in positions
            a = rand(1:Xmax)
            b = rand(1:Ymax)
        end
        
        push!(positions, (a,b))

        distances[i,i] = 0
        # Compute distances to already sampled vertices
        for j in 1:i
            distances[i,j] = dist(positions[i], positions[j])
            distances[j,i] = dist(positions[i], positions[j])
        end
    end

    return positions, distances
end

"""
Function sample_edges_to_connexity!
Brief: Samples edges on the grid until the resulting graph is connex
    Algorithm : creates adges between all the nodes at a distance Dmax 
                and increases Dmax until the graph is connex
Args:
    N : number of nodes
    positions : position of each node on the grid
    distances : matrix of distances between each pair of vertices
    edges : Vector{Tuple{Int, Int}} of edges
Returns:
    Vector{Tuple{Int, Int}} of edges
"""
function sample_edges_to_connexity(N, positions, distances)
    
    edges = Vector{Tuple{Int, Int}}()
    
    Dmax = 0
    while is_graph_connex(N, edges) == false && Dmax < 30
        
        Dmax += 1
        for i in 1:N, j in (i+1):N
            if Dmax - 1 < distances[i,j] <= Dmax
                push!(edges, (i,j))
            end
        end
    end
    return edges
end

"""
Function is_graph_connex
Brief: Check connexity of a graph given by a set of N vertices and a list of edges
Args:
    N: Number of vertices, vertices are labelled from 1 to N
    edges: Vector{Tuple{Int, Int}} of edges
Returns:
    True if graph is connex, else False
"""
function is_graph_connex(N, edges)

    seen_vertices = Set([1])

    # List of neihbors of seen_vertices to check
    ongoing_neighbors = [1]
    current_node = 0
    while length(ongoing_neighbors) > 0
        # Take an element of ongoing_neighbors
        current_node = pop!(ongoing_neighbors)

        for (i,j) in edges
            if i == current_node && (j in seen_vertices) == false
                push!(ongoing_neighbors, j)
                push!(seen_vertices, j)
            elseif j == current_node && (i in seen_vertices) == false
                push!(ongoing_neighbors, i)
                push!(seen_vertices, i)
            end
        end

        if length(seen_vertices) == N
            return true
        end
    end

    if length(seen_vertices) == N
        return true
    end
    return false

end

"""
Function remove_edge_to_density
Brief: Try removing edges to a graph to reach a given density
Args:
    N: number of vertices
    edges: edges as pairs of vertices
    density: density to reach in percent - 0:100
    max_try: maximum number of successive tries to remove arcs resulting in a non convex graph
Returns:
    nothing
"""
function remove_edge_to_density!(N, edges, recquired_density; max_try=100)
    
    successive_fail = 0
    current_density = 100*2*length(edges) / (N*(N-1))
    while current_density > recquired_density
        
        successive_fail = 0
        

        id_to_remove = rand(1:length(edges))
        edge_to_remove = edges[id_to_remove]
        while is_graph_connex(N, deleteat!(edges, id_to_remove)) == false
            
            # If not connex, place back the edge and try another one
            push!(edges, edge_to_remove)
            successive_fail += 1

            if successive_fail == max_try
                println("Could not reach density ", recquired_density, "%, final density = ", current_density, "%.")
                return nothing
            end            

            id_to_remove = rand(1:length(edges))
            edge_to_remove = edges[id_to_remove]
        end

        # Edge removed and graph is still connex
        # Computing new graph density
        current_density = 100*2*length(edges) / (N*(N-1))
    end

    if current_density > recquired_density
        println("Could not reach density ", density, "%, final density = ", current_density, "%.")
    end
    return nothing
end


"""
function create_network
Brief: Creates a random network graph
Args:
    N: number of nodes
    graph_density: percentage of arcs (100 meaning the complete graph)
    seed: random seed for reproductibility
    gridRatio: ratio between X length and Y length for the map on which points are sampled
    proportion: proportion of the points of the map on which there will be nodes, caracterized the size of the map
    drawGraph: bool to draw the graph in an external pdf file
Returns:
    postitions, edges
"""
function create_network(N, graph_density, seed; gridRatio = 1.5, proportion=0.005, drawGraph = true, plotGraph = false)
    
    seed >= 0 ? Random.seed!(seed) : nothing

    # Sample positions of the nodes on a map
    positions, distances = sample_positions(N)
    edges = sample_edges_to_connexity(N, positions, distances)
    remove_edge_to_density!(N, edges, graph_density)

    if drawGraph == true || plotGraph == true
        graph_plot_final = plot(positions, seriestype = :scatter, showaxis = false, grid = false, ticks = false, legend = false)
        for (i,j) in edges
            plot!([positions[i][1], positions[j][1]], [positions[i][2], positions[j][2]])
        end

        if plotGraph == true
            display(graph_plot_final)
            println("Press a key to continue")
            readline()
        end

        if drawGraph == true
            savefig(graph_plot_final, "graph.pdf")
        end
    end

    return Network(N, positions, edges)
end

#########################################################################################
# Main program
#########################################################################################
function main()
    
    if "?" in ARGS || "-h" in ARGS
        help()
    end

    N = 100
    graph_density = 2.5
    seed = -1
    for i in 1:length(ARGS)
        if ARGS[i] == "-N"
            N = parse(Int, ARGS[i+1])
        elseif ARGS[i] == "-d"
            graph_density = parse(Float64, ARGS[i+1])
        elseif ARGS[i] == "-p"
            proportion = parse(Float64, ARGS[i+1])
        elseif ARGS[i] == "--seed"
            seed = parse(Int, ARGS[i+1])
        end
    end
    println("coucou")
    create_network(N, graph_density, seed, plotGraph = true, drawGraph = false)
end

if PROGRAM_FILE == "graphGenerator.jl"
    main()
end