module SolutionQuery
using JuMP

struct BinaryVar
    name::Symbol
    index::Vector{Union{String,Int}}
    value::Real
end

struct Solution
    binaries::Vector{BinaryVar}
    hc::Real
    ctpv::Real
end

function split_token(line)
    if isempty(split(line, ":")[2])
        return [line, ""]
    else
        return split(line, ":")
    end
end

function convert_index(x)
    value = tryparse(Int, string(x))
    if isnothing(value)
        return string(x)
    else
        return value
    end
end


function get_binaries(lines)::Vector{BinaryVar}
    vars::Vector{BinaryVar} = []
    name = ""
    symbols::Vector{Union{Symbol,Int}} = []
    for (token, value) in lines
        if occursin(":", token)
            global name = split(token, ":")[1]
            global symbols = []
            continue
        else
            global symbols = split(token, ",") |> xs -> map(convert_index, xs)
        end
        var = BinaryVar(Symbol(name), symbols, parse(Int, value))
        push!(vars, var)
    end
    return vars
end

function read_jumpsol_file(path)::Solution
    lines = []
    open(path) do io
        lines = readlines(io)
    end
    lines = map(strip, lines) |> filter(!isempty)

    line_hc = popfirst!(lines)
    value_hc = parse(Float64, split(line_hc, ":")[2])

    line_ctpv = popfirst!(lines)
    value_ctpv = parse(Float64, split(line_ctpv, ":")[2])

    lines = map(x -> replace(x, ";" => ""), lines) |> xs -> map(split_token, xs)
    return Solution(get_binaries(lines), value_hc, value_ctpv)

end


function set_value!(model, var::BinaryVar)
    try
        JuMP.fix(model[var.name][var.index...,], var.value; force=true)
    catch e
        if isa(e, KeyError)
            JuMP.fix(model[var.name][CartesianIndex(var.index...)], var.value; force=true)
        else
            throw(e)
        end
    end
    
end


end # SolutionParser