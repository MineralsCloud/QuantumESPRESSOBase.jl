module Setters

export VerbositySetter

export batchset

abstract type BatchSetter end
struct VerbositySetter{T} <: BatchSetter
    function VerbositySetter{T}() where {T}
        T ∈ (:high, :low) ||
        throw(ArgumentError("the type parameter must be either `:high` or `:low`!"))
        return new()
    end # function VerbositySetter{T}
end
VerbositySetter(x) = VerbositySetter{x}()

function batchset end

end
