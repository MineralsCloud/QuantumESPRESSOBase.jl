"""
# module Namelists



# Examples

```jldoctest
julia>
```
"""
module Namelists

using OrderedCollections: OrderedDict

using QuantumESPRESSOBase: InputEntry, titleof, to_qe

import QuantumESPRESSOBase

export to_dict, dropdefault

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

include("PWscf.jl")
include("CP.jl")
include("PHonon.jl")

function QuantumESPRESSOBase.to_qe(
    nml::Namelist;
    indent = ' '^4,
    delim = ' ',
    newline = '\n',
)
    namelist_name = titleof(nml)
    content = to_qe(dropdefault(nml); indent = indent, delim = delim)
    return "&$namelist_name" * newline * content * '/'
end

end
