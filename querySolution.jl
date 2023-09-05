import Gurobi
import JuMP
include("src/HCPlanning.jl")
HCP = HCPlanning
Query = HCP.SolutionQuery

mode = HCP.ev


max_j = 10
if mode == HCP.gen
    global max_j = 2
end

ci_pv = zeros(max_j)
for j = 1:max_j
    path = "solutions/24_bus_$mode$j.jump_sol"

    sol = Query.read_jumpsol_file(path)
    model = HCP.build_model("data/24bus/", Gurobi.Optimizer; mode=HCP.ev)
    HCP.set_mult_objective(model, Gurobi)
    for var in sol.binaries
        Query.set_value!(model, var)
    end
    x = JuMP.all_variables(model)
    JuMP.set_silent(model)
    HCP.run_optimizer!(model, x)

    i = 7.1 / 100
    ci = JuMP.value.(model[:cᴵₜ])
    ci_pv[j] = sum(ci[t] / (1 + i)^t for t = 1:3)
end

for j=1:max_j in 
    println("$mode $j $(ci_pv[j]/1e6)")
end

