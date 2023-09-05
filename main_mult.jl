begin
    import HCPlanning
    import Gurobi
    using Dates
    using LoggingExtras
    HCP = HCPlanning
    import JuMP
    MOI = JuMP.MOI


    sys = 138
    ctpv_hc = 100e8 # INITIAL (must be a bigger value)
    mode = HCP.both

    path2main = nothing
    path2small = nothing
    sol_name = nothing
    logger = nothing

    if sys == 24
        path2main = "data/24bus/" #"data/138bus_4stages/"
        path2small = "data/24bus_1stage/"
        sol_name = "24_bus_$mode"
        logger = FileLogger("info_mult_24_$mode.log")
    elseif sys == 54
        path2main = "data/54bus_1stage/"
        path2small = "data/54bus_1stage/"
        sol_name = "54_bus_$mode"
        logger = FileLogger("info_mult_54_$mode.log")
    elseif sys == 138
        path2main = "data/138bus_1stage/"
        path2small = "data/138bus_1stage/"
        sol_name = "138_bus_$mode"
        logger = FileLogger("info_mult_138_$mode.log")
    else
        throw(InvalidStateException("Not implemented!", sys))
    end



    #######################################################################################
    function logg(msg)
        with_logger(logger) do
            end_msg = "[ $(now()) ] " * msg
            @info(end_msg)
        end
    end

    function config_solver!(model)
        JuMP.set_optimizer_attribute(model, "MIPGap", 1e-6)
        JuMP.set_optimizer_attribute(model, "Presolve", 2)
        JuMP.set_optimizer_attribute(model, "IntegralityFocus", 1)
        JuMP.set_optimizer_attribute(model, "NumericFocus", 3)
    end

    function solve_small!(model_big, path)
        model = HCP.build_model(path, Gurobi.Optimizer; mode=mode)
        #JuMP.set_silent(model)
        config_solver!(model)
        x_small = JuMP.all_variables(model)
        HCP.set_hc_obj(model)
        _ = HCP.run_optimizer!(model, x_small)

        for sym in [:xˡₛᵣₖₜ, :xˢˢₛₜ, :xᴺᵀₛₖₜ, :xᵖₛₖₜ, :yˡₛᵣₖₜ, :yᵖₛₖₜ, :yᵗʳₛₖₜ]
            for key in eachindex(model[sym])
                value = JuMP.value(model[sym][key])
                variable_big = model_big[sym][key]
                JuMP.set_start_value(variable_big, value)
            end
        end
    end

    logg("Start!")
    model = HCP.build_model(path2main, Gurobi.Optimizer; mode=mode)
    x = JuMP.all_variables(model)

    HCP.set_cptv_obj(model)
    config_solver!(model)
    _ = HCP.run_optimizer!(model, x)
    ctpv_hc_min = HCP.get_ctpv(model) 
    logg("MIN CTPV: $ctpv_hc_min")
    
    HCP.set_mult_objective(model, Gurobi)
    x = JuMP.all_variables(model)
    # Change the objective function 
    logg("Starting the Pareto!")
    hc = nothing
    sol_hc= nothing
    ε = 450.0
    count = 1
    while true
        if abs(ctpv_hc - ctpv_hc_min) < ε
            exit()
        end

        HCP.add_cptv_limit(model, ctpv_hc - ε)
        if isnothing(sol_hc)
            global sol_hc = HCP.run_optimizer!(model, x)
        else
            HCP.add_hc_limit(model, hc)
            global sol_hc = HCP.run_optimizer!(model, x, sol_hc)
        end

        global hc = HCP.get_hc(model)
        logg("Pareto | HC: $hc")

        global ctpv_hc = HCP.get_ctpv(model)
        logg("Pareto | CTPV: $ctpv_hc")
        HCP.out_integer_solutions(model, sol_name * "$count")
        global count = count + 1
    end
end