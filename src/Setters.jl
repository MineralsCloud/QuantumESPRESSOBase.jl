module Setters

using Unitful: AbstractQuantity

import Setfield

export VerbositySetter, FiniteTemperatureSetter, CellParametersSetter

export makelens, preset_values

abstract type BatchSetter end

struct VerbositySetter{T} <: BatchSetter
    function VerbositySetter{T}() where {T}
        T ∈ (:high, :low) ||
        throw(ArgumentError("the type parameter must be either `:high` or `:low`!"))
        return new()
    end # function VerbositySetter{T}
end
VerbositySetter(x) = VerbositySetter{x}()

struct FiniteTemperatureSetter{N} <: BatchSetter
    function FiniteTemperatureSetter{N}() where {N}
        N isa Union{Real,AbstractQuantity} ||
        throw(ArgumentError("the type parameter must be a real or a quantity!"))
        return new()
    end # function FiniteTemperatureSetter
end
FiniteTemperatureSetter(x) = FiniteTemperatureSetter{x}()

"""
    CellParametersSetter <: BatchSetter

Generate automatically a `CellParametersCard` for a `PWInput` or `CPInput` if its `cell_parameters` field is `nothing`.

Sometimes the `ibrav` field of a `PWInput` is not `0`, with its `cell_parameters` field to be empty.
But there are cases we want to write its `CellParametersCard` explicitly. This function will take a `PWInput` described
above and generate a new `PWInput` with its `ibrav = 0` and `cell_parameters` not empty.
"""
struct CellParametersSetter <: BatchSetter end

function makelens end

function preset_values end

Setfield.set(template, T::BatchSetter) =
    set(template, makelens(template, T), preset_values(template, T))

end
