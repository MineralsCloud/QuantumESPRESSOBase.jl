"""
# module Inputs



# Examples

```julia
julia>
```
"""
module Inputs

using AbInitioSoftwareBase.Inputs: Input, InputEntry, Namelist, Card, Setter, groupname
using OrderedCollections: OrderedDict

export getoption,
    optionpool,
    groupname,
    required_namelists,
    optional_namelists,
    required_cards,
    optional_cards,
    allnamelists,
    allcards

"""
    dropdefault(nml::Namelist)

Return an [`OrderedDict`](https://juliacollections.github.io/OrderedCollections.jl/latest/ordered_containers.html#OrderedDicts-and-OrderedSets-1) of non-default values of a `Namelist`.
"""
function dropdefault(nml::Namelist)
    default = typeof(nml)()  # Create a `Namelist` with all default values
    # Compare `default` with `nml`, discard the same values
    result = Iterators.filter(item -> item.second != getfield(default, item.first), nml)
    # for (drivingarg, passivearg) in _coupledargs(typeof(nml))
    # rule
    # end
    if isempty(result)
        @info "Every entry in the namelist is the default value!"
    end
    return OrderedDict{Symbol,Any}(result)
end

"""
    getoption(card::Card)

Return a `String` representing the option of a `Card`.

!!! warning
    Do not use `card.option` to access the option since it may not exist.
"""
getoption(card::Card) = hasproperty(card, :option) ? getproperty(card, :option) : ""

"""
    optionpool(T::Type{<:Card})

Return the allowed options for `Card` `T`.

# Examples
```julia
julia> optionpool(AtomicPositionsCard)
("alat", "bohr", "angstrom", "crystal", "crystal_sg")

julia> optionpool(CellParametersCard)
("alat", "bohr", "angstrom")

julia> optionpool(SpecialKPointsCard)
("tpiba", "crystal", "tpiba_b", "crystal_b", "tpiba_c", "crystal_c")
```
"""
function optionpool end

"Represent input files of executables (such as `pw.x` and `cp.x`)."
abstract type QuantumESPRESSOInput <: Input end

function allnamelists end

function allcards end

function required_namelists end

function optional_namelists end

function required_cards end

function optional_cards end

struct VerbositySetter <: Setter
    v::String
    function VerbositySetter(v)
        @assert v in ("high", "low")
        return new(v)
    end
end

include("crystallography.jl")
include("PWscf/PWscf.jl")
# include("CP/CP.jl")
include("PHonon/PHonon.jl")

end
