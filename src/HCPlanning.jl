module HCPlanning
using JuMP
import MunozDelgado2014
MD14 = MunozDelgado2014

function build_model(path2main)
    model = MD14.build_model(path2main)
    include(path2main * "/main.jl")
    JuMP.@variable(model, V_ ≤ vᴴᶜₛₜᵦ[Ωᴺ, T, B] ≤ V̅, start = Vbase)
    JuMP.@variable(model, gᴴᶜᵖₛₖₜᵦ[p=P, Ωᵖ[p], Kᵖ[p], T, B] ≥ 0)
    JuMP.@variable(model, 0 ≤ fᴴᶜˡₛᵣₖₜᵦ[l=L, s=Ωᴺ, Ωˡₛ[l][s], Kˡ[l], T, B])
    JuMP.@variable(model, 0 ≤ gᴴᶜᵗʳₛₖₜᵦ[tr=TR, Ωˢˢ, Kᵗʳ[tr], T, B])
    JuMP.@variable(model, 0 ≤ gᴴᶜₛₜ[Ωᴺ, T])

    JuMP.@expression(model, g_p_warp_hc[p=P, s=Ωᴺ, k=Kᵖ[p], t=T, b=B],
        if s ∈ Ωᵖ[p]
            gᴴᶜᵖₛₖₜᵦ[p, s, k, t, b]
        else
            0.0
        end
    )
    JuMP.@expression(model, g_tr_warp_hc[tr=TR, s=Ωᴺ, k=Kᵗʳ[tr], t=T, b=B],
        if s ∈ Ωˢˢ
            gᴴᶜᵗʳₛₖₜᵦ[tr, s, k, t, b]
        else
            0.0
        end
    )

    for s = Ωˢˢ, t = T, b = B
        JuMP.fix(vᴴᶜₛₜᵦ[s, t, b], Vˢˢ; force=true)
    end

    JuMP.@constraint(model, eq8_hc[l=L, r=Ωᴺ, s=Ωˡₛ[l][r], k=Kˡ[l], t=T, b=B],
        fᴴᶜˡₛᵣₖₜᵦ[l, s, r, k, t, b] ≤ model[:yˡₛᵣₖₜ][l, s, r, k, t] * F̅ˡₖ[l][k]
    )
    JuMP.@constraint(model, eq9_hc[tr=TR, s=Ωˢˢ, k=Kᵗʳ[tr], t=T, b=B],
        gᴴᶜᵗʳₛₖₜᵦ[tr, s, k, t, b] ≤ model[:yᵗʳₛₖₜ][tr, s, k, t] * G̅ᵗʳₖ[tr][k]
    )
    JuMP.@constraint(model, eq11_hc[s=Ωᵖ["C"], k=Kᵖ["C"], t=T, b=B],
        gᴴᶜᵖₛₖₜᵦ["C", s, k, t, b] ≤ model[:yᵖₛₖₜ]["C", s, k, t] * G̅ᵖₖ["C"][k]
    )
    JuMP.@constraint(model, eq12_hc[s=Ωᵖ["W"], k=Kᵖ["W"], t=T, b=B],
        gᴴᶜᵖₛₖₜᵦ["W", s, k, t, b] ≤ model[:yᵖₛₖₜ]["W", s, k, t] * minimum([G̅ᵖₖ["W"][k], Ĝᵂₛₖₜᵦ[s][k][t][b]])
    )

    JuMP.@constraint(model, eq14_hc[s=Ωᴺ, t=T, b=B], # Eq14 needs the follow fixes
        sum(sum(sum(
            fᴴᶜˡₛᵣₖₜᵦ[l, s, r, k, t, b] - fᴴᶜˡₛᵣₖₜᵦ[l, r, s, k, t, b]
            for r in Ωˡₛ[l][s])
                for k in Kˡ[l])
            for l in L)
        ==
        sum(sum(
            g_tr_warp_hc[tr, s, k, t, b]
            for k in Kᵗʳ[tr])
            for tr in TR)
        +
        sum(sum(
            g_p_warp_hc[p, s, k, t, b]
            for k in Kᵖ[p])
            for p in P)
        -
        μᵦ[b] * Dₛₜ[s, t] + gᴴᶜₛₜ[s, t]
    )

    for l = L, r = Ωᴺ, s = Ωˡₛ[l][r], k = Kˡ[l], t = T, b = B
        JuMP.@constraint(model,
            -Zˡₖ[l][k] * ℓₛᵣ[s, r] * fᴴᶜˡₛᵣₖₜᵦ[l, s, r, k, t, b] / Vbase + (vᴴᶜₛₜᵦ[s, t, b] - vᴴᶜₛₜᵦ[r, t, b]) ≤ H * (1 - model[:yˡₛᵣₖₜ][l, s, r, k, t])
        )
        JuMP.@constraint(model,
            Zˡₖ[l][k] * ℓₛᵣ[s, r] * fᴴᶜˡₛᵣₖₜᵦ[l, s, r, k, t, b] / Vbase - (vᴴᶜₛₜᵦ[s, t, b] - vᴴᶜₛₜᵦ[r, t, b]) ≤ H * (1 - model[:yˡₛᵣₖₜ][l, s, r, k, t])
        )
    end

    JuMP.@expression(model, hc, sum(sum(gᴴᶜₛₜ[s, t] for s ∈ Ωᴺ) for t ∈ T))
    # JuMP.@constraint(model, hc_constraint, hc ≥ 0)
    JuMP.@expression(model, ci, sum(model[:cᴵₜ][t] * ((1 + i)^-t) / i for t in T))
    JuMP.@constraint(model, ci_constraint, model[:ci] ≤ floatmax(Float64))
    JuMP.@objective(model, Max, hc - model[:ci] / 1e7)
    return model
end

function change_ctpv_constraint(model, rhs)
    JuMP.delete(model, model[:ctpv_constraint])
    JuMP.unregister(model, :ctpv_constraint)
    JuMP.@constraint(model, ctpv_constraint, hc ≤ rhs)
end

end # module
