using Compat: isnothing
using LinearAlgebra: Diagonal, det, cross

export BravaisLattice
export asfieldname, titleof, to_qe, cell_volume, reciprocalof, supercell

"""
    InputEntry

Represent any component of a `QuantumESPRESSOInput`.

Hierachy of `InputEntry`:
```
QuantumESPRESSOBase.InputEntry
├─ QuantumESPRESSOBase.Cards.Card
│  ├─ CP.AbstractCellParametersCard
│  │  ├─ CP.CellParametersCard
│  │  └─ CP.RefCellParametersCard
│  ├─ CP.AtomicForcesCard
│  ├─ CP.AtomicPositionsCard
│  ├─ CP.AtomicSpeciesCard
│  ├─ CP.AtomicVelocitiesCard
│  ├─ CP.KPointsCard
│  ├─ PHonon.AbstractCellParametersCard
│  │  └─ PHonon.CellParametersCard
│  ├─ PHonon.AtomicForcesCard
│  ├─ PHonon.AtomicPositionsCard
│  ├─ PHonon.AtomicSpeciesCard
│  ├─ PHonon.KPointsCard
│  ├─ PWscf.AbstractCellParametersCard
│  │  └─ PWscf.CellParametersCard
│  ├─ PWscf.AtomicForcesCard
│  ├─ PWscf.AtomicPositionsCard
│  ├─ PWscf.AtomicSpeciesCard
│  └─ PWscf.KPointsCard
└─ QuantumESPRESSOBase.Namelists.Namelist
   ├─ CP.CellNamelist
   ├─ CP.ControlNamelist
   ├─ CP.ElectronsNamelist
   ├─ CP.IonsNamelist
   ├─ CP.PressAiNamelist
   ├─ CP.SystemNamelist
   ├─ CP.WannierNamelist
   ├─ PHonon.DynmatNamelist
   ├─ PHonon.MatdynNamelist
   ├─ PHonon.PhNamelist
   ├─ PHonon.Q2rNamelist
   ├─ PWscf.BandsNamelist
   ├─ PWscf.CellNamelist
   ├─ PWscf.ControlNamelist
   ├─ PWscf.DosNamelist
   ├─ PWscf.ElectronsNamelist
   ├─ PWscf.IonsNamelist
   └─ PWscf.SystemNamelist
```
"""
abstract type InputEntry end

"""
    asfieldname(::Type{<:InputEntry})
    asfieldname(::InputEntry)

Return the field name of a `Namelist` or a `Card` in a `QuantumESPRESSOInput`.

# Examples

```jldoctest
julia> using QuantumESPRESSOBase; using QuantumESPRESSOBase.Namelists.PWscf: SystemNamelist

julia> asfieldname(SystemNamelist) == asfieldname(SystemNamelist()) == :system
true
```
"""
asfieldname(x::InputEntry) = asfieldname(typeof(x))

"""
    titleof(::Type{<:InputEntry})
    titleof(::InputEntry)

Return the title of the input entry in Quantum ESPRESSO.

The definition `titleof(x) = titleof(typeof(x))` is provided for convenience so that
instances can be passed instead of types.

# Examples

```jldoctest
julia> using QuantumESPRESSOBase; using QuantumESPRESSOBase.Namelists.PWscf: SystemNamelist

julia> titleof(SystemNamelist()) == titleof(SystemNamelist) == "SYSTEM"
true
```
"""
titleof(x::InputEntry) = titleof(typeof(x))

to_fortran(v::Int) = string(v)
function to_fortran(v::Float32; scientific::Bool = false)
    str = string(v)
    scientific && return replace(str, r"f"i => "e")
    return str
end
function to_fortran(v::Float64; scientific::Bool = false)
    str = string(v)
    scientific && return replace(str, r"e"i => "d")
    return string(v)
end
function to_fortran(v::Bool)
    v ? ".true." : ".false."
end
function to_fortran(v::AbstractString)
    return "'$v'"
end

"""
    to_qe(x; indent = ' '^4, delim = ' ')

Return a `String` representing the object, which is valid for Quantum ESPRESSO's input.
"""
function to_qe(dict::AbstractDict; indent = ' '^4, delim = ' ')::String
    content = ""
    f = string ∘ to_fortran
    for (key, value) in dict
        if value isa Vector
            for (i, x) in enumerate(value)
                isnothing(x) && continue
                content *= indent * join(["$key($i)", "=", "$(f(x))\n"], delim)
            end
        else
            content *= indent * join(["$key", "=", "$(f(value))\n"], delim)
        end
    end
    return content
end

function cell_volume end

struct BravaisLattice{I}
    celldm::Vector{Union{Nothing,Float64}}
    function BravaisLattice{I}(celldm) where {I}
        @assert I ∈ union(1:1:14, (-3, -5, -9, 91, -12, -13))  # It can't be `0`!
        @assert(
            if I == 14
                length(celldm) == 6
            elseif I ∈ (5, -5, 12, 13)
                4 <= length(celldm) <= 6
            elseif I ∈ (4, 6, 7, 8, 9, -9, 91, 10, 11)  # `91` is new from QE 6.4
                3 <= length(celldm) <= 6
            elseif I == -13  # `-13` is new from QE 6.4
                5 <= length(celldm) <= 6
            else
                1 <= length(celldm) <= 6
            end,
            "`celldm` must have length between 1 to 6! See `ibrav`'s doc!"
        )
        return new(celldm)
    end
end

"""
    (::BravaisLattice{I})()

Return a ``3 × 3`` matrix representing the Bravais lattice from `BravaisLattice{I}`.
"""
(bravais::BravaisLattice{1})() = bravais.celldm[1] * [
    1 0 0
    0 1 0
    0 0 1
]
(bravais::BravaisLattice{2})() = bravais.celldm[1] / 2 * [
    -1 0 1
    0 1 1
    -1 1 0
]
(bravais::BravaisLattice{3})() = bravais.celldm[1] / 2 * [
    1 1 1
    -1 1 1
    -1 -1 1
]
(bravais::BravaisLattice{-3})() = bravais.celldm[1] / 2 * [
    -1 1 1
    1 -1 1
    1 1 -1
]
(bravais::BravaisLattice{4})() =
    bravais.celldm[1] * [
        1 0 0
        -1 / 2 √3 / 2 0
        0 0 bravais.celldm[3]
    ]
function (bravais::BravaisLattice{5})()
    c = bravais.celldm[4]
    tx = sqrt((1 - c) / 2)
    ty = sqrt((1 - c) / 6)
    tz = sqrt((1 + 2c) / 3)
    return bravais.celldm[1] * [
        tx -ty tz
        0 2ty tz
        -tx -ty tz
    ]
end
function (bravais::BravaisLattice{-5})()
    ap = bravais.celldm[1] / √3
    c = bravais.celldm[4]
    ty = sqrt((1 - c) / 6)
    tz = sqrt((1 + 2c) / 3)
    u = tz - 2 * √2 * ty
    v = tz + √2 * ty
    return ap * [
        u v v
        v u v
        v v u
    ]
end
(bravais::BravaisLattice{6})() = bravais.celldm[1] * [
    1 0 0
    0 1 0
    0 0 bravais.celldm[3]
]
(bravais::BravaisLattice{7})() =
    bravais.celldm[1] / 2 * [
        1 -1 bravais.celldm[3] / bravais.celldm[1]
        1 1 bravais.celldm[3] / bravais.celldm[1]
        -1 -1 bravais.celldm[3] / bravais.celldm[1]
    ]
(bravais::BravaisLattice{8})() =
    bravais.celldm[1] * [
        1 0 0
        0 bravais.celldm[2] 0
        0 0 bravais.celldm[3]
    ]
(bravais::BravaisLattice{9})() =
    bravais.celldm[1] * [
        1 / 2 bravais.celldm[2] / 2 0
        -1 / 2 bravais.celldm[2] / 2 0
        0 0 bravais.celldm[3]
    ]
(bravais::BravaisLattice{-9})() =
    bravais.celldm[1] * [
        1 / 2 -bravais.celldm[2] / 2 0
        1 / 2 bravais.celldm[2] / 2 0
        0 0 bravais.celldm[3]
    ]
(bravais::BravaisLattice{91})() =
    bravais.celldm[1] * [
        1 0 0
        0 bravais.celldm[2] / 2 -bravais.celldm[3] / 2
        0 bravais.celldm[2] / 2 bravais.celldm[3] / 2
    ]  # New from QE 6.4
(bravais::BravaisLattice{10})() =
    bravais.celldm[1] * [
        1 / 2 0 bravais.celldm[3] / 2
        1 / 2 bravais.celldm[2] / 2 0
        0 bravais.celldm[2] / 2 bravais.celldm[3] / 2
    ]
(bravais::BravaisLattice{11})() =
    bravais.celldm[1] * [
        1 / 2 0 bravais.celldm[3] / 2
        1 / 2 bravais.celldm[2] / 2 0
        0 bravais.celldm[2] / 2 bravais.celldm[3] / 2
    ]
(bravais::BravaisLattice{12})() =
    bravais.celldm[1] * [
        1 0 0
        bravais.celldm[2] * bravais.celldm[4] bravais.celldm[2] *
                                              sqrt(1 - bravais.celldm[4]^2) 0
        0 0 bravais.celldm[3]
    ]
(bravais::BravaisLattice{-12})() =
    bravais.celldm[1] * [
        1 0 0
        0 bravais.celldm[2] 0
        bravais.celldm[3] * bravais.celldm[5] 0 bravais.celldm[3] *
                                                sqrt(1 - bravais.celldm[5]^2)
    ]
(bravais::BravaisLattice{13})() =
    bravais.celldm[1] * [
        1 / 2 0 -bravais.celldm[3] / 2
        bravais.celldm[2] * bravais.celldm[4] bravais.celldm[2] *
                                              sqrt(1 - bravais.celldm[4]^2) 0
        1 / 2 0 bravais.celldm[3] / 2
    ]
(bravais::BravaisLattice{-13})() =
    bravais.celldm[1] * [
        1 / 2 -bravais.celldm[2] / 2 0
        1 / 2 bravais.celldm[2] / 2 0
        bravais.celldm[3] * bravais.celldm[5] 0 bravais.celldm[3] *
                                                sqrt(1 - bravais.celldm[5]^2)
    ]
(bravais::BravaisLattice{14})() =
    bravais.celldm[1] * [
        1 0 0
        bravais.celldm[2] * bravais.celldm[6] bravais.celldm[2] *
                                              sqrt(1 - bravais.celldm[6]^2) 0
        bravais.celldm[3] * bravais.celldm[5] bravais.celldm[3] * (
            bravais.celldm[4] - bravais.celldm[5] * bravais.celldm[6]
        ) / sqrt(1 - bravais.celldm[6]^2) bravais.celldm[3] * sqrt(
            1 + 2 * bravais.celldm[4] * bravais.celldm[5] * bravais.celldm[6] -
            bravais.celldm[4]^2 - bravais.celldm[5]^2 - bravais.celldm[6]^2,
        ) / sqrt(1 - bravais.celldm[6]^2)
    ]

function reciprocalof(mat::AbstractMatrix)
    @assert size(mat) == (3, 3)
    volume = det(mat)
    a1, a2, a3 = mat[1, :], mat[2, :], mat[3, :]
    return 2π / volume * [cross(a2, a3) cross(a3, a1) cross(a1, a2)]
end # function reciprocalof
"""
    reciprocalof(bravais::BravaisLattice)

Return a ``3 × 3`` matrix representing the reciprocal lattice from a `BravaisLattice`.
"""
function reciprocalof(bravais::BravaisLattice)
    return reciprocalof(bravais())
end # function reciprocalof

"""
    supercell(cell::AbstractMatrix, expansion::AbstractMatrix{<:Integer})

Allow the supercell to be a tilted extension of `cell`.
"""
function supercell(cell::AbstractMatrix, expansion::AbstractMatrix{<:Integer})
    @assert(det(expansion) != 0, "matrix `expansion` cannot be a singular integer matrix!")
    return expansion * cell
end # function supercell
"""
    supercell(cell::AbstractMatrix, expansion::AbstractVector{<:Integer})

Return a supercell based on `cell` and expansion coefficients.
"""
function supercell(cell::AbstractMatrix, expansion::AbstractVector{<:Integer})
    @assert length(expansion) == 3
    return supercell(cell, Diagonal(expansion))
end # function supercell
