begin 
import HCPlanning
import Gurobi
import HiGHS
using Dates
HCP = HCPlanning
import JuMP
MOI = JuMP.MOI


path2main = "data/24bus_1stage"
model = HCP.build_model(path2main)
JuMP.set_optimizer(model, Gurobi.Optimizer, add_bridges=false)
JuMP.set_optimizer_attribute(model, "MIPGap", 1e-6)
JuMP.set_silent(model)
x = JuMP.all_variables(model)

# Optimize to cost
HCP.set_cptv_obj(model)
sol_cost = HCP.run_optimizer!(model, x)

ci_cost = HCP.get_ci(model)
@show ci_cost

ctpv_cost = HCP.get_ctpv(model)
@show ctpv_cost

sol_cost, hc_cost = HCP.calc_hc(model, x, sol_cost)
@show hc_cost

println("----")
# Optimize to HC
HCP.set_hc_obj(model)
sol_hc = HCP.run_optimizer!(model, x, sol_cost)

hc_hc = HCP.get_hc(model)
@show hc_hc
ci_hc = HCP.get_ci(model)
@show ci_hc

sol_hc, ctpv_hc = HCP.calc_ctpv(model, x, sol_hc)
@show ctpv_hc
# Change the objective function 
ε = 1000

while true
    println("----")
    HCP.set_hc_obj(model, 1/(100*ctpv_hc))
    HCP.add_cptv_limit(model, ctpv_hc - ε)
    global sol_hc = HCP.run_optimizer!(model, x, sol_hc)

    global hc_hc = HCP.get_hc(model)
    @show hc_hc
    global ci_hc = HCP.get_ci(model)
    @show ci_hc
    
    sol_hc_l, ctpv_hc_l = HCP.calc_ctpv(model, x, sol_hc)
    global sol_hc = sol_hc_l
    global ctpv_hc = ctpv_hc_l
    @show ctpv_hc
    
end
end
# # Enter in the loop
# while true
#     # Set initial solution
#     try
#         global x_sol = run_optimizer!(model, x, x_sol)
#         check(model)
#     catch e
#         if isa(y, ErrorException)
#             break
#         end
#     end

#     ci = JuMP.value(model[:ci])
#     hc = JuMP.value(model[:hc])
#     ctpv = JuMP.value(model[:cᵀᴾⱽ]) / 1e6
#     println(now())
#     @show ci
#     @show hc
#     @show ctpv
#     println("-------------")
#     sol = JuMP.value.(x)
#     HCP.change_ci_constraint(model, ci - ε)
#     JuMP.set_start_value.(x, sol)


# end
