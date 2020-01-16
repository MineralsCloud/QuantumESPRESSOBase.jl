using Compat: isnothing
using LinearAlgebra: Diagonal, det, cross

export asfieldname,
    titleof, to_qe, cell_volume, direct_lattice, reciprocal_lattice, supercell

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

"""
    direct_lattice(ibrav::Integer, celldm::AbstractVector)

Return a ``3 × 3`` matrix representing the Bravais lattice (in real space) from `ibrav` and `celldm`.
"""
direct_lattice(ibrav::Integer, celldm::AbstractVector) = _direct_lattice(Val(ibrav), celldm)
# These are helper methods and should not be exported.
_direct_lattice(::Val{1}, celldm::AbstractVector) = celldm[1] * [
    1 0 0
    0 1 0
    0 0 1
]
_direct_lattice(::Val{2}, celldm::AbstractVector) = celldm[1] / 2 * [
    -1 0 1
    0 1 1
    -1 1 0
]
_direct_lattice(::Val{3}, celldm::AbstractVector) = celldm[1] / 2 * [
    1 1 1
    -1 1 1
    -1 -1 1
]
_direct_lattice(::Val{-3}, celldm::AbstractVector) = celldm[1] / 2 * [
    -1 1 1
    1 -1 1
    1 1 -1
]
_direct_lattice(::Val{4}, celldm::AbstractVector) =
    celldm[1] * [
        1 0 0
        -1 / 2 √3 / 2 0
        0 0 celldm[3]
    ]
function _direct_lattice(::Val{5}, celldm::AbstractVector)
    c = celldm[3]
    tx = sqrt((1 - c) / 2)
    ty = sqrt((1 - c) / 6)
    tz = sqrt((1 + 2c) / 3)
    return celldm[1] * [
        tx -ty tz
        0 2ty tz
        -tx -ty tz
    ]
end
function _direct_lattice(::Val{-5}, celldm::AbstractVector)
    ap = celldm[1] / √3
    c = celldm[3]
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
_direct_lattice(::Val{6}, celldm::AbstractVector) = celldm[1] * [
    1 0 0
    0 1 0
    0 0 celldm[3]
]
_direct_lattice(::Val{7}, celldm::AbstractVector) =
    celldm[1] / 2 * [
        1 -1 celldm[3] / celldm[1]
        1 1 celldm[3] / celldm[1]
        -1 -1 celldm[3] / celldm[1]
    ]
_direct_lattice(::Val{8}, celldm::AbstractVector) =
    celldm[1] * [
        1 0 0
        0 celldm[2] 0
        0 0 celldm[3]
    ]
_direct_lattice(::Val{9}, celldm::AbstractVector) =
    celldm[1] * [
        1 / 2 celldm[2] / 2 0
        -1 / 2 celldm[2] / 2 0
        0 0 celldm[3]
    ]
_direct_lattice(::Val{-9}, celldm::AbstractVector) =
    celldm[1] * [
        1 / 2 -celldm[2] / 2 0
        1 / 2 celldm[2] / 2 0
        0 0 celldm[3]
    ]
_direct_lattice(::Val{10}, celldm::AbstractVector) =
    celldm[1] * [
        1 / 2 0 celldm[3] / 2
        1 / 2 celldm[2] / 2 0
        0 celldm[2] / 2 celldm[3] / 2
    ]
_direct_lattice(::Val{11}, celldm::AbstractVector) =
    celldm[1] * [
        1 / 2 0 celldm[3] / 2
        1 / 2 celldm[2] / 2 0
        0 celldm[2] / 2 celldm[3] / 2
    ]
_direct_lattice(::Val{12}, celldm::AbstractVector) =
    celldm[1] * [
        1 0 0
        celldm[2] * celldm[4] celldm[2] * sqrt(1 - celldm[4]^2) 0
        0 0 celldm[3]
    ]
_direct_lattice(::Val{-12}, celldm::AbstractVector) =
    celldm[1] * [
        1 0 0
        0 celldm[2] 0
        celldm[3] * celldm[5] 0 celldm[3] * sqrt(1 - celldm[5]^2)
    ]
_direct_lattice(::Val{13}, celldm::AbstractVector) =
    celldm[1] * [
        1 / 2 0 -celldm[3] / 2
        celldm[2] * celldm[4] celldm[2] * sqrt(1 - celldm[4]^2) 0
        1 / 2 0 celldm[3] / 2
    ]
_direct_lattice(::Val{-13}, celldm::AbstractVector) =
    celldm[1] * [
        1 / 2 -celldm[2] / 2 0
        1 / 2 celldm[2] / 2 0
        celldm[3] * celldm[5] 0 celldm[3] * sqrt(1 - celldm[5]^2)
    ]
_direct_lattice(::Val{14}, celldm::AbstractVector) =
    celldm[1] * [
        1 0 0
        celldm[2] * celldm[6] celldm[2] * sqrt(1 - celldm[6]^2) 0
        celldm[3] * celldm[5] celldm[3] * (celldm[4] - celldm[5] * celldm[6]) /
                              sqrt(1 - celldm[6]^2) celldm[3] * sqrt(
            1 + 2 * celldm[4] * celldm[5] * celldm[6] - celldm[4]^2 - celldm[5]^2 -
            celldm[6]^2,
        ) / sqrt(1 - celldm[6]^2)
    ]

function reciprocal_lattice(mat::AbstractMatrix)
    @assert size(mat) == (3, 3)
    volume = det(mat)
    a1, a2, a3 = mat[1, :], mat[2, :], mat[3, :]
    return 2π / volume * [cross(a2, a3) cross(a3, a1) cross(a1, a2)]
end # function reciprocal_lattice
"""
    reciprocal_lattice(ibrav::Integer, celldm::AbstractVector)

Return a ``3 × 3`` matrix representing the reciprocal lattice from `ibrav` and `celldm`.
"""
function reciprocal_lattice(ibrav::Integer, celldm::AbstractVector)
    return reciprocal_lattice(direct_lattice(ibrav, celldm))
end # function reciprocal_lattice

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
