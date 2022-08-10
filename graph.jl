#########################################################################################
# Import management
#########################################################################################
include("import_management.jl")
include("options.jl")

#########################################################################################
# Network type
#########################################################################################
struct Edge
    from::Int
    to::Int
end

struct Network
    N::Int # Both the number of nodes and their indices
    positions::Vector{Tuple{Int, Int}} # Positions of the nodes on a map
    edges::Vector{Edge} # list of undirected edges
    n_edges::Int # Number of edges
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
    edges : Vector{Edge} of edges
Returns:
    Vector{Edge} of edges
"""
function sample_edges_to_connexity(N, positions, distances)
    
    edges = Vector{Edge}()
    
    Dmax = 0
    while is_graph_connex(N, edges) == false
        
        Dmax += 1
        for i in 1:N, j in (i+1):N
            if Dmax - 1 < distances[i,j] <= Dmax
                push!(edges, Edge(i,j))
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
    edges: Vector{Edge} of edges
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

        for e in edges
            if e.from == current_node && (e.to in seen_vertices) == false
                push!(ongoing_neighbors, e.to)
                push!(seen_vertices, e.to)
            elseif e.to == current_node && (e.from in seen_vertices) == false
                push!(ongoing_neighbors, e.from)
                push!(seen_vertices, e.from)
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
    edges: edges (struct with from::Int and to::Int)
    density: density to reach in percent - 0:100
    max_try: maximum number of successive tries to remove arcs resulting in a non convex graph
Returns:
    nothing
"""
function remove_edge_to_density!(N, edges, recquired_density; max_try=100)
    
    successive_fail = 0
    current_density = 100*2*length(edges) / (N*(N-1))
    println("Initial graph density = ", current_density)
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
function create_network(N, graph_density, seed; gridRatio = 1.5, proportion = 0.01, drawGraph = true, plotGraph = false)
    
    seed >= 0 ? Random.seed!(seed) : nothing

    # Sample positions of the nodes on a map
    positions, distances = sample_positions(N; gridRatio = gridRatio, proportion = proportion)
    edges = sample_edges_to_connexity(N, positions, distances)

    # Forcing to draw the initial graph
    #graph_plot_final = plot(positions, seriestype = :scatter, showaxis = false, grid = false, ticks = false, legend = false)
    #for e in edges
    #    plot!([positions[e.from][1], positions[e.to][1]], [positions[e.from][2], positions[e.to][2]])
    #end
    #savefig(graph_plot_final, "graph_init_N30.pdf")

    #for graph_density in [15,10,8,6]
    #    remove_edge_to_density!(N, edges, graph_density)

        # Forcing to draw the initial graph
    #    graph_plot_final = plot(positions, seriestype = :scatter, showaxis = false, grid = false, ticks = false, legend = false)
    #    for e in edges
    #        plot!([positions[e.from][1], positions[e.to][1]], [positions[e.from][2], positions[e.to][2]])
    #    end
    #    savefig(graph_plot_final, "graph_final_N30_dens" * string(graph_density) * ".pdf")
    #end
    #exit()

    remove_edge_to_density!(N, edges, graph_density)

    if drawGraph == true || plotGraph == true
        graph_plot_final = plot(positions, seriestype = :scatter, showaxis = false, grid = false, ticks = false, legend = false)
        for e in edges
            plot!([positions[e.from][1], positions[e.to][1]], [positions[e.from][2], positions[e.to][2]])
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

    return Network(N, positions, edges, length(edges))
end

function create_network(options; drawGraph = true, plotGraph = false)
    create_network(
        options.N, options.graph_density, options.seed; 
        gridRatio = 1.5, proportion=0.005, 
        drawGraph = drawGraph, plotGraph = plotGraph
    )
end

#########################################################################################
# Main program
#########################################################################################
function main()
    
    if "?" in ARGS || "-h" in ARGS
        help()
    end

    N = 15
    graph_density = 30
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
    
    net = create_network(N, graph_density, seed, plotGraph = true, drawGraph = false)
    
end

if PROGRAM_FILE == "graphGenerator.jl"
    main()
end