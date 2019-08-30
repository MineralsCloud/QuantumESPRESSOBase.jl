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
    to_dict(nml; defaultorder = true)

Convert a `Namelist` to a dictionary.

# Arguments
- `nml::Namelist`: the namelist to be converted.
- `defaultorder::Bool = true`: whether or not use the default order of parameters in QE's docs.
"""
function to_dict(nml::Namelist; defaultorder::Bool = true)
    dict = (defaultorder ? OrderedDict{Symbol,Any}() : Dict{Symbol,Any}())
    for n in propertynames(nml)
        dict[n] = getproperty(nml, n)
    end
    return dict
end

"""
    dropdefault(nml::Namelist)

Return an `AbstractDict` of non-default values of a `Namelist`.
"""
function dropdefault(nml::Namelist)
    default = typeof(nml)()  # Create a `Namelist` with all default values
    # Compare `default` with `nml`, discard the same values
    result = filter!(item -> item.second != getfield(default, item.first), to_dict(nml))
    isempty(result) && @info "Every entry in the namelist is the default value!"
    return result
end

"""
    Base.dump(path, nml::Namelist)

Serialize a `Namelist` to `path`. Currently, only JSON and YAML formats are supported.
"""
function Base.dump(path::AbstractPath, nml::Namelist)
    exists(path) || touch(path)  # If the file does not exist, create one
    entries = Dict(key => to_fortran(value) for (key, value) in to_dict(nml))
    iswritable(path) || error("File $(path) is not writable!")
    open(path, "r+") do io
        @match extension(path) begin  # If the extension of the file is:
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
