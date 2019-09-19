export AbstractInput, namelists, cards

abstract type AbstractInput end

# A helper function to implement `namelists` and `cards`. It should not be exported.
_filterfields(obj, ::Type{T}) where {T} =
    filter(x -> isa(x, T), map(x -> getfield(obj, x), fieldnames(typeof(obj))) |> collect)

"""
    namelists(input::PWscfInput)

Return a vector of `Namelist`s of a `PWscfInput`.
"""
namelists(input::PWscfInput) = _filterfields(input, Namelist)

"""
    cards(input::PWscfInput)

Return a vector of `Card`s of a `PWscfInput`.
"""
cards(input::PWscfInput) = _filterfields(input, Card)
