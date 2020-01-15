module QuantumESPRESSOBase

using Compat: isnothing

export asfieldname, titleof, to_qe, cell_volume, bravais_lattice

"Represent any component of a `QuantumESPRESSOInput`."
abstract type InputEntry end

"""
    asfieldname(::Type{<:InputEntry})
    asfieldname(::InputEntry)

Return the field name of a `Namelist` or a `Card` in a `QuantumESPRESSOInput`.

# Examples

```jldoctest
julia> using QuantumESPRESSOBase; using QuantumESPRESSOBase.Namelists.PWscf: SystemNamelist

julia> asfieldname(SystemNamelist) == asfieldname(SystemNamelist()) == :system
true
```
"""
asfieldname(x::InputEntry) = asfieldname(typeof(x))

"""
    titleof(::Type{<:InputEntry})
    titleof(::InputEntry)

Return the title of the input entry in Quantum ESPRESSO.

The definition `titleof(x) = titleof(typeof(x))` is provided for convenience so that
instances can be passed instead of types.

# Examples

```jldoctest
julia> using QuantumESPRESSOBase; using QuantumESPRESSOBase.Namelists.PWscf: SystemNamelist

julia> titleof(SystemNamelist()) == titleof(SystemNamelist) == "SYSTEM"
true
```
"""
titleof(x::InputEntry) = titleof(typeof(x))

to_fortran(v::Int) = string(v)
function to_fortran(v::Float32; scientific::Bool = false)
    str = string(v)
    scientific && return replace(str, r"f"i => "e")
    return str
end
function to_fortran(v::Float64; scientific::Bool = false)
    str = string(v)
    scientific && return replace(str, r"e"i => "d")
    return string(v)
end
function to_fortran(v::Bool)
    v ? ".true." : ".false."
end
function to_fortran(v::AbstractString)
    return "'$v'"
end

"""
    to_qe(x; indent = ' '^4, delim = ' ')

Return a string representing the object, valid form Quantum ESPRESSO's input.
"""
function to_qe(dict::AbstractDict; indent = ' '^4, delim = ' ')::String
    content = ""
    f = string ∘ to_fortran
    for (key, value) in dict
        if value isa Vector
            for (i, x) in enumerate(value)
                isnothing(x) && continue
                content *= indent * join(["$key($i)", "=", "$(f(x))\n"], delim)
            end
        else
            content *= indent * join(["$key", "=", "$(f(value))\n"], delim)
        end
    end
    return content
end

function cell_volume end

include("bravais_lattice.jl")
include("Setters.jl")
include("Namelists/Namelists.jl")
include("Cards/Cards.jl")
include("Inputs/Inputs.jl")
include("CLI.jl")

end # module
