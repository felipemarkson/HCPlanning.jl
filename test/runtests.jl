import HCPlanning
import Gurobi
HCP = HCPlanning
import JuMP
MOI = JuMP.MOI

function runtests()
    path2main = "../data/24bus_1stage"
    model = HCP.build_model(path2main)
    JuMP.set_optimizer(model, Gurobi.Optimizer, add_bridges=false)
    JuMP.optimize!(model)
    @show JuMP.value(model[:hc])
    @show JuMP.value(model[:ci])
end

runtests()