module Setters

using Unitful: AbstractQuantity

export VerbositySetter, FiniteTemperatureSetter, CellParametersSetter

export batchset

abstract type BatchSetter end

struct VerbositySetter{T} <: BatchSetter
    function VerbositySetter{T}() where {T}
        T ∈ (:high, :low) || throw(ArgumentError("the type parameter must be either `:high` or `:low`!"))
        return new()
    end # function VerbositySetter{T}
end
VerbositySetter(x) = VerbositySetter{x}()

struct FiniteTemperatureSetter{N} <: BatchSetter
    function FiniteTemperatureSetter{N}() where {N}
        N isa Union{Real,AbstractQuantity} || throw(ArgumentError("the type parameter must be a real or a quantity!"))
        return new()
    end # function FiniteTemperatureSetter
end
FiniteTemperatureSetter(x) = FiniteTemperatureSetter{x}()

struct CellParametersSetter <: BatchSetter end

function batchset end

end
