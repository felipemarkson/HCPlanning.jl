begin 
import HCPlanning
import Gurobi
using Dates
using LoggingExtras
HCP = HCPlanning
import JuMP
MOI = JuMP.MOI


path2main = "data/24bus/" #"data/138bus_4stages/"
path2small = "data/24bus_1stage/"
logger = FileLogger("info_mult_24_both_start_again.log")



function logg(msg)
    with_logger(logger) do
        end_msg = "[ $(now()) ] " * msg
        @info(end_msg)
    end
end

function solve_small!(model_big, path)
    model = HCP.build_model(path, Gurobi.Optimizer)
    JuMP.set_silent(model)
    JuMP.set_optimizer_attribute(model, "Presolve", 2)
    JuMP.set_optimizer_attribute(model, "IntegralityFocus", 3)
    JuMP.set_optimizer_attribute(model, "MIPGap", 1e-6)
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
    JuMP.set_optimizer_attribute(model, "IntegralityFocus", 3)
end

logg("Start!")
model = HCP.build_model(path2main, Gurobi.Optimizer)
HCP.set_mult_objective(model, Gurobi)
x = JuMP.all_variables(model)
config_solver!(model)

# Initial solution
solve_small!(model, path2small)
logg("Small Solved!")

sol_hc = nothing
ctpv_hc = 3e8 # INITIAL (must be a bigger value)

# Change the objective function 
logg("Starting the Pareto!")
ε = 1e3

while true
    HCP.add_cptv_limit(model, ctpv_hc - ε)
    if isnothing(sol_hc)
        global sol_hc = HCP.run_optimizer!(model, x)
    else
        global sol_hc = HCP.run_optimizer!(model, x, sol_hc)
    end

    hc = HCP.get_hc(model)
    logg("Pareto | HC: $hc")

    global ctpv_hc = HCP.get_ctpv(model)
    logg("Pareto | CTPV: $ctpv_hc")

end
end