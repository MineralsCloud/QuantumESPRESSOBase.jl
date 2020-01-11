module Setters

using Setfield: set
using Unitful: AbstractQuantity

import Setfield

export VerbositySetter,
    FiniteTemperatureSetter, CellParametersSetter, AlatPressSetter, LensMaker
export make, preset_values

abstract type BatchSetter end

struct VerbositySetter{T} <: BatchSetter
    function VerbositySetter{T}() where {T}
        T ∈ (:high, :low) ||
        throw(ArgumentError("the type parameter must be either `:high` or `:low`!"))
        return new()
    end
end
VerbositySetter(x) = VerbositySetter{x}()

struct FiniteTemperatureSetter{N} <: BatchSetter
    function FiniteTemperatureSetter{N}() where {N}
        N isa Union{Real,AbstractQuantity} ||
        throw(ArgumentError("the type parameter must be a real or a quantity!"))
        return new()
    end
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

struct AlatPressSetter <: BatchSetter end

struct LensMaker{S<:BatchSetter,T} end
LensMaker{S}(::T) where {S<:BatchSetter,T} = LensMaker{S,T}()
LensMaker(S::BatchSetter) = T -> LensMaker{S,T}()
LensMaker(::S, ::T) where {S<:BatchSetter,T} = LensMaker{S,T}()

make(l::LensMaker) = error("`make` is not defined for `$l`!")
make(maker::LensMaker{<:BatchSetter,T}, makers::LensMaker{<:BatchSetter,T}...) where {T} =
    make(maker) ∘ make(makers...)

function preset_values end

Setfield.set(template, T::BatchSetter) =
    set(template, make(LensMaker(T, template)), preset_values(T, template))

end
