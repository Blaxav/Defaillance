#########################################################################################
# Import management
#########################################################################################
if "--install-packages" in ARGS
    import Pkg
    Pkg.add("Random")
    Pkg.add("Plots")
    Pkg.add("JuMP")
    Pkg.add("Gurobi")
end

if "--compile" in ARGS
    include("compile_plotpackage.jl")
    println("Package Plot.jl compiled. Relaunch without '--compile' and with '--sysimage sys_plots.so'")
    exit
end

try
    using Random
    using Plots
    using Printf
    using JuMP, Gurobi
catch
    println("ERROR : Some recquired packages not found.")
    println("        Run with option --install-packages to install them")
    exit()
end
