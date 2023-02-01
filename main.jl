begin 
import HCPlanning
import Gurobi
import HiGHS
using Dates
HCP = HCPlanning
import JuMP
MOI = JuMP.MOI

GRB_DEFAULT = 0
GRB_FEASIBLE = 1
GRB_OPTIMAL = 2
GRB_BESTBOUND = 3


function solve_small!(model_big, path2small)
    model = HCP.build_model(path2small)
    JuMP.set_optimizer(model, Gurobi.Optimizer, add_bridges=false)
    # JuMP.set_silent(model)
    JuMP.set_optimizer_attribute(model, "MIPGap", 1/100)
    JuMP.set_optimizer_attribute(model, "MIPFocus", 1)
    x_small = JuMP.all_variables(model)
    HCP.set_cptv_obj(model)
    _ = HCP.run_optimizer!(model, x_small)

    for sym in [:xˡₛᵣₖₜ, :xˢˢₛₜ, :xᴺᵀₛₖₜ, :xᵖₛₖₜ,:yˡₛᵣₖₜ, :yᵖₛₖₜ, :yᵗʳₛₖₜ]
        for key in eachindex(model[sym])
            value = JuMP.value(model[sym][key])
            variable_big = model_big[sym][key]        
            JuMP.set_start_value(variable_big, value)
        end
    end
end


path2main = "data/138bus_4stages/"
path2small = "data/138bus_1stage/"
model = HCP.build_model(path2main)
JuMP.set_optimizer(model, Gurobi.Optimizer)
JuMP.set_optimizer_attribute(model, "MIPGap", 1e-6)
JuMP.set_optimizer_attribute(model, "Presolve", 2)
JuMP.set_optimizer_attribute(model, "Heuristics", 10/100)
# JuMP.set_optimizer_attribute(model, "NodefileStart", 1.0)
# JuMP.set_optimizer_attribute(model, "MIPFocus", GRB_OPTIMAL)
# JuMP.set_silent(model)
x = JuMP.all_variables(model)

# Optimize to cost
HCP.set_cptv_obj(model)
solve_small!(model, path2small)

sol_cost = HCP.run_optimizer!(model, x)

ci_cost = HCP.get_ci(model)
@show ci_cost

ctpv_cost = HCP.get_ctpv(model)
@show ctpv_cost

sol_cost, hc_cost = HCP.calc_hc(model, x, sol_cost)
end
# @show hc_cost

# println("----")
# # Optimize to HC
# HCP.set_hc_obj(model)
# sol_hc = HCP.run_optimizer!(model, x, sol_cost)

# hc_hc = HCP.get_hc(model)
# @show hc_hc
# ci_hc = HCP.get_ci(model)
# @show ci_hc

# sol_hc, ctpv_hc = HCP.calc_ctpv(model, x, sol_hc)
# @show ctpv_hc
# # Change the objective function 
# ε = 1000

# while true
#     println("----")
#     HCP.set_hc_obj(model, 1/(100*ctpv_hc))
#     HCP.add_cptv_limit(model, ctpv_hc - ε)
#     global sol_hc = HCP.run_optimizer!(model, x, sol_hc)

#     global hc_hc = HCP.get_hc(model)
#     @show hc_hc
#     global ci_hc = HCP.get_ci(model)
#     @show ci_hc
    
#     sol_hc_l, ctpv_hc_l = HCP.calc_ctpv(model, x, sol_hc)
#     global sol_hc = sol_hc_l
#     global ctpv_hc = ctpv_hc_l
#     @show ctpv_hc
    
# end
# end