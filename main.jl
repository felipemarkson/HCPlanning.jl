begin 
import HCPlanning
import Gurobi
using Dates
using LoggingExtras
HCP = HCPlanning
import JuMP
MOI = JuMP.MOI

GRB_DEFAULT = 0
GRB_FEASIBLE = 1
GRB_OPTIMAL = 2
GRB_BESTBOUND = 3


date_format = "yyyy-mm-dd HH:MM:SS"
logger = FileLogger("info.log")
function logg(msg)
    with_logger(logger) do
        end_msg = "[ $(now()) ] " * msg
        @info(end_msg)
    end
end

function solve_small!(model_big, path2small)
    model = HCP.build_model(path2small, Gurobi.Optimizer)
    JuMP.set_silent(model)
    JuMP.set_optimizer_attribute(model, "Presolve", 2)
    JuMP.set_optimizer_attribute(model, "MIPGap", 1e-6)
    # JuMP.set_optimizer_attribute(model, "MIPFocus", 1)
    x_small = JuMP.all_variables(model)
    HCP.set_hc_obj(model)
    _ = HCP.run_optimizer!(model, x_small)

    for sym in [:xˡₛᵣₖₜ, :xˢˢₛₜ, :xᴺᵀₛₖₜ, :xᵖₛₖₜ,:yˡₛᵣₖₜ, :yᵖₛₖₜ, :yᵗʳₛₖₜ]
        for key in eachindex(model[sym])
            value = JuMP.value(model[sym][key])
            variable_big = model_big[sym][key]        
            JuMP.set_start_value(variable_big, value)
        end
    end
end

function config_solver!(model)
    JuMP.set_optimizer_attribute(model, "MIPGap", 1e-6)
    JuMP.set_optimizer_attribute(model, "Presolve", 2)
end

logg("Start!")
path2main = "data/24bus"
path2small = "data/24bus_1stage/"
model = HCP.build_model(path2main, Gurobi.Optimizer)
HCP.set_hc_obj(model)
# JuMP.set_silent(model)
x = JuMP.all_variables(model)
config_solver!(model)

# Initial solution
solve_small!(model, path2small)
logg("Small Solved!")

# Optimize to HC
sol_hc = HCP.run_optimizer!(model, x)

hc_hc = HCP.get_hc(model)
logg("Pareto 0| HC: $hc_hc")

sol_hc, ctpv_hc = HCP.calc_ctpv(model, x, sol_hc)
logg("Pareto 0 | CTPV: $ctpv_hc")

# # Change the objective function 
logg("Starting the Pareto!")
ε = 1000

while true
    HCP.set_hc_obj(model)
    HCP.add_cptv_limit(model, ctpv_hc - ε)
    global sol_hc = HCP.run_optimizer!(model, x, sol_hc)

    global hc_hc = HCP.get_hc(model)
    logg("Pareto | HC: $hc_hc")

    sol_hc_l, ctpv_hc_l = HCP.calc_ctpv(model, x, sol_hc)
    logg("Pareto | CTPV: $ctpv_hc_l")
    global sol_hc = sol_hc_l
    global ctpv_hc = ctpv_hc_l
    
end
end