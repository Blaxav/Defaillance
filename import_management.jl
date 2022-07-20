#########################################################################################
# Import management
#########################################################################################
if "--install-packages" in ARGS
    import Pkg
    Pkg.add("Random")
    Pkg.add("Plots")
    Pkg.add("JuMP")
    #Pkg.add("Gurobi")
    Pkg.add("CPLEX")
    Pkg.add("Dualization")
    Pkg.add("BilevelJuMP")
end

if "--compile" in ARGS
    include("compile_external_packages.jl")
    println("Packages compiled. Relaunch without '--compile' and with '--sysimage lib.dll or lib.so'")
    exit
end

try
    using Random
    using Plots
    using Printf
    using JuMP#, Gurobi
    using CPLEX
    using BilevelJuMP
catch
    println("ERROR : Some recquired packages not found.")
    println("        Run with option --install-packages to install them")
    exit()
end
