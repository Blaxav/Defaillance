# Installation of PackageCompiler
if "--install-package" in ARGS
    import Pkg
    Pkg.add("PackageCompiler")
end

try
    # Creates an image of Plot.jl (causing long loading and launching times) in sys_plots.so
    using PackageCompiler
    #create_sysimage(:Plots, sysimage_path="compiled_lib/graph_libs.dll", precompile_execution_file="graphGenerator.jl")
    create_sysimage([:Plots, :JuMP, :CPLEX, :BilevelJuMP], sysimage_path="compiled_libs/probGen_libs.dll", precompile_execution_file="problemGenerator.jl")
catch
    println("ERROR : Some recquired packages not found.")
    println("        Run with option --install-packages to install them")
    exit()
end