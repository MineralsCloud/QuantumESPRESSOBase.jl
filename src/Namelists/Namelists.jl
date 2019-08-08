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
using Parameters: type2dict

using QuantumESPRESSOBase: InputEntry

export Namelist, to_dict, dropdefault
# ============================================================================ #


abstract type Namelist <: InputEntry end

function to_dict(nml::Namelist)::Dict{Symbol,Any}
    return type2dict(nml)
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
        if extension(path) == "json"
            JSON.print(io, entries)
        elseif extension(path) == "yaml" || extension(path) == "yml"
            @warn "Currently not supported!"
        else
            error("Unknown extension type given!")
        end
    end
end # function Base.dump


# ============================================================================ #
#                                    Include                                   #
# ============================================================================ #
include("PW.jl")
include("PH.jl")
# ============================================================================ #

end
