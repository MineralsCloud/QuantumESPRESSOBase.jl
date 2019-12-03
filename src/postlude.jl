using Compat: isnothing, eachrow
using Fortran90Namelists.JuliaToFortran: to_fortran
using Parameters: @with_kw

using .Namelists: Namelist, to_dict, dropdefault
using .Namelists.PWscf
using .Namelists.PHonon
using .Cards
using .Cards.PWscf
using .Cards.PHonon
using .Inputs: QuantumESPRESSOInput, namelists, cards
using .Inputs.PWscf

export PWExecutable
export asfieldname, titleof, to_qe

# The command line args are from https://github.com/QEF/q-e/blob/033aee7/Modules/command_line_options.f90#L84-L154.
@with_kw struct PWExecutable
    path::String = "pw.x"
    inp::String
    nimage::Int = 0
    npool::Int = 0
    ntg::Int = 0
    nyfft::Int = 0
    nband::Int = 0
    ndiag::Int = 0
end

"""
    asfieldname()



# Arguments

# Examples

```jldoctest
julia>
```
"""
asfieldname(::Type{T}) where {T<:InputEntry} = error("Undefined for entry $(nameof(T))!")
asfieldname(::Type{<:ControlNamelist}) = :control
asfieldname(::Type{<:SystemNamelist}) = :system
asfieldname(::Type{<:ElectronsNamelist}) = :electrons
asfieldname(::Type{<:IonsNamelist}) = :ions
asfieldname(::Type{<:CellNamelist}) = :cell
asfieldname(::Type{<:PhNamelist}) = :inputph
asfieldname(::Type{<:Q2rNamelist}) = :input
asfieldname(::Type{<:MatdynNamelist}) = :input
asfieldname(::Type{<:DynmatNamelist}) = :input
asfieldname(::Type{<:AtomicSpeciesCard}) = :atomic_species
asfieldname(::Type{<:AtomicPositionsCard}) = :atomic_positions
asfieldname(::Type{<:KPointsCard}) = :k_points
asfieldname(::Type{<:CellParametersCard}) = :cell_parameters

"""
    titleof(::Type{<:InputEntry})

Return the title of the input entry in Quantum ESPRESSO.

The definition `titleof(x) = titleof(typeof(x))` is provided for convenience so that
instances can be passed instead of types.

# Examples

```jldoctest
julia> titleof(SystemNamelist())
"SYSTEM"

julia> titleof(SystemNamelist)
"SYSTEM"
```
"""
titleof(x::InputEntry) = titleof(typeof(x))
titleof(::Type{T}) where {T<:InputEntry} = error("Undefined for entry $(nameof(T))!")
titleof(::Type{<:ControlNamelist}) = "CONTROL"
titleof(::Type{<:SystemNamelist}) = "SYSTEM"
titleof(::Type{<:ElectronsNamelist}) = "ELECTRONS"
titleof(::Type{<:IonsNamelist}) = "IONS"
titleof(::Type{<:CellNamelist}) = "CELL"
titleof(::Type{<:PhNamelist}) = "INPUTPH"
titleof(::Type{<:Q2rNamelist}) = "INPUT"
titleof(::Type{<:MatdynNamelist}) = "INPUT"
titleof(::Type{<:DynmatNamelist}) = "INPUT"
titleof(::Type{<:AtomicSpeciesCard}) = "ATOMIC_SPECIES"
titleof(::Type{<:AtomicPositionsCard}) = "ATOMIC_POSITIONS"
titleof(::Type{<:KPointsCard}) = "K_POINTS"
titleof(::Type{<:CellParametersCard}) = "CELL_PARAMETERS"

"""
    to_qe(x, indent::AbstractString = "    ", sep::AbstractString = " ")

Return a string representing the object, valid form Quantum ESPRESSO's input.
"""
function to_qe(
    dict::AbstractDict;
    indent::AbstractString = "    ",
    sep::AbstractString = " ",
)::String
    content = ""
    f = string ∘ to_fortran
    for (key, value) in dict
        if value isa Vector
            for (i, x) in enumerate(value)
                isnothing(x) && continue
                content *= indent * join(["$key($i)", "=", "$(f(x))\n"], sep)
            end
        else
            content *= indent * join(["$key", "=", "$(f(value))\n"], sep)
        end
    end
    return content
end
function to_qe(
    nml::Namelist;
    indent::AbstractString = "    ",
    sep::AbstractString = " ",
    verbose::Bool = false,
)
    namelist_name = (titleof ∘ typeof)(nml)
    f = verbose ? to_dict : dropdefault
    inner_content = to_qe(f(nml); indent = indent, sep = sep)
    return """
    &$namelist_name
    $inner_content/
    """
end
function to_qe(data::AtomicSpecies; sep::AbstractString = " ")::String
    return join([getfield(data, i) for i in 1:nfields(data)], sep)
end
function to_qe(
    card::AtomicSpeciesCard;
    indent::AbstractString = "    ",
    sep::AbstractString = " ",
)
    return """
    ATOMIC_SPECIES
    $(join([indent * to_qe(x, sep = sep) for x in card.data], "\n"))
    """
end
function to_qe(data::AtomicPosition; sep::AbstractString = " ", verbose::Bool = false)::String
    verbose && return join([data.atom; data.pos; data.if_pos], sep)
    return join([data.atom; data.pos], sep)
end
function to_qe(
    card::AtomicPositionsCard;
    indent::AbstractString = "    ",
    sep::AbstractString = " ",
)
    return """
    ATOMIC_POSITIONS$sep{ $(card.option) }
    $(join([indent * to_qe(x, sep = sep) for x in card.data], "\n"))
    """
end
function to_qe(
    card::CellParametersCard;
    indent::AbstractString = "    ",
    sep::AbstractString = " ",
)
    return """
    CELL_PARAMETERS$sep{ $(card.option) }
    $(join([indent * join(row, sep) for row in eachrow(card.data)], "\n"))
    """
end
function to_qe(data::MonkhorstPackGrid; sep::AbstractString = " ")::String
    return join([data.grid; data.offsets], sep)
end
function to_qe(data::GammaPoint)
    return ""
end
function to_qe(data::SpecialKPoint; sep::AbstractString = " ")::String
    return join([data.coordinates; data.weight], sep)
end
function to_qe(
    card::KPointsCard;
    indent::AbstractString = "    ",
    sep::AbstractString = " ",
)::String
    content = "K_POINTS$sep{ $(card.option) }\n"
    if card.option in ("gamma", "automatic")
        content *= indent * to_qe(card.data) * "\n"
    else  # option in ("tpiba", "crystal", "tpiba_b", "crystal_b", "tpiba_c", "crystal_c")
        content *= "$(length(card.data))\n"
        for x in card.data
            content *= indent * to_qe(x, sep = sep) * "\n"
        end
    end
    return content
end
function to_qe(
    card::QPointsSpecsCard;
    indent::AbstractString = "    ",
    sep::AbstractString = " ",
)::String
    content = "$(length(card.data))\n"
    for p in card.data
        content *= indent * join([p.coordinates; p.weight], sep) * "\n"
    end
    return content
end
function to_qe(
    input::QuantumESPRESSOInput;
    indent::AbstractString = "    ",
    sep::AbstractString = " ",
    verbose::Bool = false,
)::String
    content = ""
    for namelist in namelists(input)
        content *= to_qe(namelist, indent = indent, sep = sep, verbose = verbose)
    end
    for card in cards(input)
        content *= to_qe(card, indent = indent, sep = sep)
    end
    return content
end
