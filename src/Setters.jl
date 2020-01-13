module Setters

using Setfield: set
using Unitful: AbstractQuantity

import Setfield

export VerbositySetter,
    FiniteTemperatureSetter,
    CellParametersSetter,
    AlatPressSetter,
    CalculationSetter,
    LensMaker
export make, preset_values

abstract type Setter end

struct VerbositySetter{T} <: Setter
    function VerbositySetter{T}() where {T}
        T ∈ (:high, :low) ||
        throw(ArgumentError("the type parameter must be either `:high` or `:low`!"))
        return new()
    end
end
VerbositySetter(x) = VerbositySetter{x}()

struct FiniteTemperatureSetter{N} <: Setter
    function FiniteTemperatureSetter{N}() where {N}
        N isa Union{Real,AbstractQuantity} ||
        throw(ArgumentError("the type parameter must be a real or a quantity!"))
        return new()
    end
end
FiniteTemperatureSetter(x) = FiniteTemperatureSetter{x}()

"""
    CellParametersSetter <: Setter

Generate automatically a `CellParametersCard` for a `PWInput` or `CPInput` if its `cell_parameters` field is `nothing`.

Sometimes the `ibrav` field of a `PWInput` is not `0`, with its `cell_parameters` field to be empty.
But there are cases we want to write its `CellParametersCard` explicitly. This function will take a `PWInput` described
above and generate a new `PWInput` with its `ibrav = 0` and `cell_parameters` not empty.
"""
struct CellParametersSetter <: Setter end

struct AlatPressSetter <: Setter end

struct CalculationSetter{T} <: Setter
    function CalculationSetter{T}() where {T}
        @assert T ∈ (
            :scf,
            :cp,
            :nscf,
            :bands,
            :relax,
            :md,
            :vc_relax,
            :vc_md,
            :vc_cp,
            :cp_wf,
            :vc_cp_wf,
        )
        return new()
    end
end
CalculationSetter(x) = CalculationSetter{x}()

struct LensMaker{S<:Setter,T} end
LensMaker{S}(::T) where {S<:Setter,T} = LensMaker{S,T}()
LensMaker(S::Setter) = T -> LensMaker{S,T}()
LensMaker(::S, ::T) where {S<:Setter,T} = LensMaker{S,T}()

make(l::LensMaker) = error("`make` is not defined for `$l`!")
make(maker::LensMaker{<:Setter,T}, makers::LensMaker{<:Setter,T}...) where {T} =
    make(maker) ∘ make(makers...)

function preset_values end

function upgrade end

Setfield.set(template, T::Setter) =
    set(template, make(LensMaker(T, template)), preset_values(T, template))

end
