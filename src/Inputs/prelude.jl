export AbstractInput, namelists, cards

abstract type AbstractInput end

# A helper function to implement `namelists` and `cards`. It should not be exported.
_filterfields(obj, ::Type{T}) where {T} =
    filter(x -> isa(x, T), map(x -> getfield(obj, x), fieldnames(typeof(obj))) |> collect)

"""
    namelists(input::AbstractInput)

Return a vector of `Namelist`s of a `AbstractInput`'s subtypes.
"""
namelists(input::AbstractInput) = _filterfields(input, Namelist)

"""
    cards(input::AbstractInput)

Return a vector of `Card`s of a `AbstractInput`'s subtypes.
"""
cards(input::AbstractInput) = _filterfields(input, Card)
