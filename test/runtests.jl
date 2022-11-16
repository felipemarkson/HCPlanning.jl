import HCPlanning
import Gurobi
HCP = HCPlanning
import JuMP
MOI = JuMP.MOI

function runtests()
    path2main = "../data/24bus"
    model = HCP.build_model(path2main)
    HCP.change_ci_constraint(model, 6.1834282649378115e6 - 100)
    JuMP.set_optimizer(model, Gurobi.Optimizer, add_bridges=false)
    JuMP.set_optimizer_attribute(model,"MIPGap", 0)
    JuMP.optimize!(model)
    @show JuMP.value(model[:hc])
    @show JuMP.value(model[:ci])
end

runtests()