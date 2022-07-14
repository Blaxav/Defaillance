using JuMP, CPLEX

model = Model(CPLEX.Optimizer)


# invest variables
@variable(model, 0 <= x)
@variable(model, 0 <= y)

@constraint(model, ctr1, y <= -x + 5)
@constraint(model, ctr2, y <= -10*x + 32)

@objective(model, Min, x + y)

optimize!(model)

println("Optimal value = ", objective_value(model))
println("x = ", value(x))
println("y = ", value(y))