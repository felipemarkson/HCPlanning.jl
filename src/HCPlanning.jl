module HCPlanning
using JuMP
import MunozDelgado2014
MD14 = MunozDelgado2014

function build_model(path2main)
    model = MD14.build_model(path2main)
    include(path2main * "/main.jl")
    αᴴᶜ = [-1, 1]
    Sᴴᶜ = [1, 2]

    Sᴿᵂ = [1, 2]
    αᴿᵂ = [0, 1]

    JuMP.@variable(model, hc_fˡₛᵣₖₜᵦₕ[l=L, s=Ωᴺ, Ωˡₛ[l][s], Kˡ[l], T, B, Sᴿᵂ, Sᴴᶜ] ≥ 0)
    JuMP.@variable(model, hc_gᴰₛₖₜᵦₕ[p=D, Ωᵖ[p], Kᵖ[p], T, B, Sᴿᵂ, Sᴴᶜ] ≥ 0)
    JuMP.@variable(model, hc_gᵗʳₛₖₜᵦₕ[tr=TR, Ωˢˢ, Kᵗʳ[tr], T, B, Sᴿᵂ, Sᴴᶜ] ≥ 0)
    JuMP.@variable(model, V_ ≤ hc_vₛₜᵦₕ[Ωᴺ, T, B, Sᴿᵂ, Sᴴᶜ] ≤ V̅, start = Vbase)
    JuMP.@variable(model, gᴴᶜ[T] ≥ 0)
    gᴿᵂₛₖₜᵦ = model[:gᴿᵂₛₖₜᵦ]

    JuMP.@expression(model, g_p_warp_hc[p=P, s=Ωᴺ, k=Kᵖ[p], t=T, b=B, h=Sᴿᵂ, u=Sᴴᶜ],
        if s ∈ Ωᵖ[p]
            if p in RW
                αᴿᵂ[h] * gᴿᵂₛₖₜᵦ[p, s, k, t, b]
            elseif p in D
                hc_gᴰₛₖₜᵦₕ[p, s, k, t, b, h, u]
            else
                error("Invalid")
            end
        else
            0.0
        end
    )
    JuMP.@expression(model, g_tr_warp_hc[tr=TR, s=Ωᴺ, k=Kᵗʳ[tr], t=T, b=B, h=Sᴿᵂ, u=Sᴴᶜ],
        if s ∈ Ωˢˢ
            hc_gᵗʳₛₖₜᵦₕ[tr, s, k, t, b, h, u]
        else
            0.0
        end
    )

    JuMP.@expression(model, g_hc_warp[s=Ωᴺ, t=T], # Only if there is load
        if s ∈ Ωᴸᴺₜ[t]
            gᴴᶜ[t]
        else
            0.0
        end
    )


    for s = Ωˢˢ, t = T, b = B, h = Sᴿᵂ, u = Sᴴᶜ
        JuMP.fix(hc_vₛₜᵦₕ[s, t, b, h, u], Vˢˢ; force=true)
    end

    JuMP.@constraint(model, eq8_hc[l=L, r=Ωᴺ, s=Ωˡₛ[l][r], k=Kˡ[l], t=T, b=B, h=Sᴿᵂ, u=Sᴴᶜ],
        hc_fˡₛᵣₖₜᵦₕ[l, s, r, k, t, b, h, u] ≤ model[:yˡₛᵣₖₜ][l, s, r, k, t] * F̅ˡₖ[l][k]
    )

    JuMP.@constraint(model, eq9_hc[tr=TR, s=Ωˢˢ, k=Kᵗʳ[tr], t=T, b=B, h=Sᴿᵂ, u=Sᴴᶜ],
        hc_gᵗʳₛₖₜᵦₕ[tr, s, k, t, b, h, u] ≤ model[:yᵗʳₛₖₜ][tr, s, k, t] * G̅ᵗʳₖ[tr][k]
    )

    JuMP.@constraint(model, eq11_hc[p=D, s=Ωᵖ[p], k=Kᵖ[p], t=T, b=B, h=Sᴿᵂ, u=Sᴴᶜ],
        hc_gᴰₛₖₜᵦₕ[p, s, k, t, b, h, u] ≤ model[:yᵖₛₖₜ][p, s, k, t] * G̅ᴰₚₖ[p][k]
    )

    JuMP.@constraint(model, eq14_hc[s=Ωᴺ, t=T, b=B, h=Sᴿᵂ, u=Sᴴᶜ],
        sum(sum(sum(
            hc_fˡₛᵣₖₜᵦₕ[l, s, r, k, t, b, h, u] - hc_fˡₛᵣₖₜᵦₕ[l, r, s, k, t, b, h, u]
            for r in Ωˡₛ[l][s])
                for k in Kˡ[l])
            for l in L)
        ==
        sum(sum(
            g_tr_warp_hc[tr, s, k, t, b, h, u]
            for k in Kᵗʳ[tr])
            for tr in TR)
        +
        sum(sum(
            g_p_warp_hc[p, s, k, t, b, h, u]
            for k in Kᵖ[p])
            for p in P)
        -
        μᵦ[b] * Dₛₜ[s, t] + αᴴᶜ[u] * g_hc_warp[s, t]
    )

    #Eq 15 and 16
    JuMP.@constraint(model, eq16_1_hc[l=L, r=Ωᴺ, s=Ωˡₛ[l][r], k=Kˡ[l], t=T, b=B, h=Sᴿᵂ, u=Sᴴᶜ],
        -Zˡₖ[l][k] * ℓₛᵣ[s, r] * hc_fˡₛᵣₖₜᵦₕ[l, s, r, k, t, b, h, u] / Vbase + (hc_vₛₜᵦₕ[s, t, b, h, u] - hc_vₛₜᵦₕ[r, t, b, h, u]) ≤ H * (1 - model[:yˡₛᵣₖₜ][l, s, r, k, t])
    )
    JuMP.@constraint(model, eq16_2_hc[l=L, r=Ωᴺ, s=Ωˡₛ[l][r], k=Kˡ[l], t=T, b=B, h=Sᴿᵂ, u=Sᴴᶜ],
        Zˡₖ[l][k] * ℓₛᵣ[s, r] * hc_fˡₛᵣₖₜᵦₕ[l, s, r, k, t, b, h, u] / Vbase - (hc_vₛₜᵦₕ[s, t, b, h, u] - hc_vₛₜᵦₕ[r, t, b, h, u]) ≤ H * (1 - model[:yˡₛᵣₖₜ][l, s, r, k, t])
    )

    for i in 1:length(T[1:end-1])
        t = T[i]
        t_next = T[i+1]
        JuMP.@constraint(model, gᴴᶜ[t] ≤ gᴴᶜ[t_next])
    end


    JuMP.@expression(model, hc, sum(gᴴᶜ[t] for t ∈ T))
    # JuMP.@constraint(model, hc_constraint, hc ≥ 0)
    JuMP.@expression(model, ci, sum(model[:cᴵₜ][t] * ((1 + i)^-t) / i for t in T))

    return model
end

function set_bi_obj(model, weight_ctpv, weight_hc)
    JuMP.@objective(model, Max, weight_hc * model[:hc] - weight_ctpv * model[:cᵀᴾⱽ])
end

function set_cptv_obj(model)
    set_bi_obj(model, 1, 0)
end

function set_hc_obj(model, weight_ctpv=0)
    set_bi_obj(model, weight_ctpv, 1)
end

function remove(model, symbol)
    JuMP.delete(model, model[symbol])
    JuMP.unregister(model, symbol)
end

function check(model)
    if JuMP.termination_status(model) != MOI.OPTIMAL
        error("The model was not solved correctly.")
    end
end

function run_optimizer!(model, x)
    JuMP.optimize!(model)
    check(model)
    return JuMP.value.(x)
end
function run_optimizer!(model, x, x_init)
    JuMP.set_start_value.(x, x_init)
    return run_optimizer!(model, x)
end

function get_hc(model)
    return JuMP.value(model[:hc])
end
function get_ctpv(model)
    return JuMP.value(model[:cᵀᴾⱽ])
end
function get_ci(model)
    return JuMP.value(model[:ci])
end

function fixing_binaries!(model)
    variables = []
    values = []
    for sym in [:xˡₛᵣₖₜ, :xˢˢₛₜ, :xᴺᵀₛₖₜ, :xᵖₛₖₜ,:yˡₛᵣₖₜ, :yᵖₛₖₜ, :yᵗʳₛₖₜ]
        for key in eachindex(model[sym])
            variable = model[sym][key]
            value = JuMP.value(variable)
            push!(variables, variable)
            push!(values, value)    
        end
    end
    JuMP.unset_binary.(variables)   
    JuMP.fix.(variables, values)
end

function unfixing_binaries!(model)
    for sym in [:xˡₛᵣₖₜ, :xˢˢₛₜ, :xᴺᵀₛₖₜ, :xᵖₛₖₜ,:yˡₛᵣₖₜ, :yᵖₛₖₜ, :yᵗʳₛₖₜ]
        for key in eachindex(model[sym])
            variable = model[sym][key]     
            JuMP.unfix(variable)
            JuMP.set_binary(variable)
        end
    end
end

function calc_ctpv(model, x, x_init)
    ci_fix = get_ci(model)
    hc_fix = get_hc(model)
    # fixing_binaries!(model)
    # JuMP.@constraint(model, ci_constraint, model[:ci] <= ci_fix)
    JuMP.@constraint(model, hc_constraint, model[:hc] >= hc_fix)
    set_cptv_obj(model)
    sol = run_optimizer!(model, x, x_init)
    ctpv = get_ctpv(model)
    remove(model, :hc_constraint)
    # remove(model, :ci_constraint)
    # unfixing_binaries!(model)
    return sol, ctpv
end

function calc_hc(model, x, x_init)
    ci_fix = get_ci(model)
    ctpv_fix = get_ctpv(model)
    fixing_binaries!(model)
    # JuMP.@constraint(model, ci_constraint, model[:ci] <= ci_fix)
    JuMP.@constraint(model, ctpv_constraint, model[:cᵀᴾⱽ] <= ctpv_fix)
    set_hc_obj(model)
    sol = run_optimizer!(model, x, x_init)
    hc = get_hc(model)
    remove(model, :ctpv_constraint)
    # remove(model, :ci_constraint)
    unfixing_binaries!(model)
    return sol, hc
end


function add_cptv_limit(model, rhs)
    try
        remove(model, :ctpv_constraint)
    catch e
        nothing
    end
    JuMP.@constraint(model, ctpv_constraint, model[:cᵀᴾⱽ] ≤ rhs)
end

end # module
