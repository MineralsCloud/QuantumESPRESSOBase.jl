"""
# module Namelists



# Examples

```jldoctest
julia>
```
"""
module Namelists

# ============================================================================ #
#                               Import and export                              #
# ============================================================================ #
using FilePaths: AbstractPath, extension, exists
using Fortran90Namelists.JuliaToFortran: to_fortran
import JSON
using MLStyle: @match
using DataStructures: OrderedDict

using QuantumESPRESSOBase: InputEntry

export Namelist, to_dict, dropdefault
# ============================================================================ #


abstract type Namelist <: InputEntry end

"""
    to_dict(nml; keeporder = true)

Convert a `Namelist` to a dictionary.

# Arguments
- `nml::Namelist`: the namelist to be converted.
- `keeporder::Bool = true`: whether or not use the default order of parameters in QE's docs.
"""
function to_dict(nml::Namelist; keeporder::Bool = true)
    dict = (keeporder ? OrderedDict{Symbol,Any}() : Dict{Symbol,Any}())
    for n in propertynames(nml)
        dict[n] = getproperty(nml, n)
    end
    return dict
end

function dropdefault(nml::Namelist)
    default = typeof(nml)()
    result = filter!(item -> item.second != getfield(default, item.first), to_dict(nml))
    isempty(result) && @info "Every entry in the namelist is the default value!"
    return result
end

function Base.dump(path::AbstractPath, nml::Namelist)
    exists(path) || touch(path)
    entries = Dict(key => to_fortran(value) for (key, value) in to_dict(nml))
    iswritable(path) || error("File $(path) not writable!")
    open(path, "r+") do io
        @match extension(path) begin
            "json" => JSON.print(io, entries)
            "yaml" || "yml" => @warn "Currently not supported!"
            _ => error("Unknown extension type given!")
        end
    end
end # function Base.dump


# ============================================================================ #
#                                    Include                                   #
# ============================================================================ #
include("PWscf.jl")
include("CP.jl")
include("PHonon.jl")
# ============================================================================ #

end
