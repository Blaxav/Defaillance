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
    dirs = filter(x -> isdir(joinpath(path, x)), readdir(path))

    for dir in dirs
        println(dir)
    end
end