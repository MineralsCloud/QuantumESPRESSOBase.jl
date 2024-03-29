using AbInitioSoftwareBase: Input, InputEntry, Namelist, Card, Setter, groupname
using OrderedCollections: OrderedDict

export QuantumESPRESSOInput, hasoption, getoption, optionpool, groupname

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
getoption(card::Card) = hasoption(card) ? card.option : nothing

hasoption(::Type{T}) where {T} = hasfield(T, :option)
hasoption(card::Card) = hasproperty(card, :option)

"""
    optionpool(card::Card)
    optionpool(T::Type{<:Card})

Return the allowed options for a `Card` or a `Card` type.

# Examples
```jldoctest
julia> optionpool(AtomicPositionsCard)
("alat", "bohr", "angstrom", "crystal", "crystal_sg")

julia> optionpool(CellParametersCard)
("alat", "bohr", "angstrom")

julia> optionpool(SpecialPointsCard)
("tpiba", "crystal", "tpiba_b", "crystal_b", "tpiba_c", "crystal_c")
```
"""
function optionpool end

"Represent input files of executables (such as `pw.x` and `cp.x`)."
abstract type QuantumESPRESSOInput <: Input end

struct VerbositySetter <: Setter
    v::String
    function VerbositySetter(v)
        @assert v in ("high", "low")
        return new(v)
    end
end
