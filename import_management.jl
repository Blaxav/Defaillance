#########################################################################################
# Import management
#########################################################################################
if "--install-packages" in ARGS
    import Pkg
    Pkg.add("Random")
    Pkg.add("Plots")
    Pkg.add("JuMP")
    Pkg.add("Gurobi")
    Pkg.add("Dualization")
    Pkg.add("BilevelJuMP")
end

if "--compile" in ARGS
    include("compile_plotpackage.jl")
    println("Packages compiled. Relaunch without '--compile' and with '--sysimage sys_plots.so'")
    exit
end

try
    using Random
    using Plots
    using Printf
    using JuMP, Gurobi
    using BilevelJuMP
catch
    println("ERROR : Some recquired packages not found.")
    println("        Run with option --install-packages to install them")
    exit()
end
