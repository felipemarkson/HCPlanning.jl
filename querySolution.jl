import Gurobi
import JuMP
include("src/HCPlanning.jl")
HCP = HCPlanning
Query = HCP.SolutionQuery

path = "solutions/24_bus_ev1.jump_sol"

sol = Query.read_jumpsol_file(path)
model = HCP.build_model("data/24bus/", Gurobi.Optimizer; mode=HCP.ev)
HCP.set_mult_objective(model, Gurobi)
for var in sol.binaries
    Query.set_value!(model, var)
end
x = JuMP.all_variables(model)
JuMP.set_silent(model)
HCP.run_optimizer!(model, x)

JuMP.value.(model[:cᴵₜ])
