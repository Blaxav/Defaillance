# Installation of PackageCompiler
if "--install-package" in ARGS
    import Pkg
    Pkg.add("PackageCompiler")
end

try
    # Creates an image of Plot.jl (causing long loading and launching times) in sys_plots.so
    using PackageCompiler
    create_sysimage([:Plots, :JuMP, :Gurobi, :CPLEX, :BilevelJuMP, :Dualization, :Random], sysimage_path="compiled_libs/unsupplied_libs.dll", precompile_execution_file="unsupplied_network_management.jl")
catch
    println("ERROR : Some recquired packages not found.")
    println("        Run with option --install-packages to install them")
    exit()
end