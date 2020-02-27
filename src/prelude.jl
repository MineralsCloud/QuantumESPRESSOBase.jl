using LinearAlgebra: Diagonal, det, cross

using Compat: isnothing
using Crystallography:
    BravaisLattice,
    Triclinic,
    Monoclinic,
    Orthorhombic,
    Tetragonal,
    Cubic,
    Trigonal,
    Hexagonal,
    Centering,
    BaseCentered,
    Primitive,
    BodyCentered,
    FaceCentered,
    RhombohedralCentered,
    BaseCentered,
    Lattice,
    CellParameters

import Crystallography

export asfieldname, titleof, to_qe, reciprocalof, supercell

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

Crystallography.BravaisLattice(ibrav::Integer) = BravaisLattice(Val(ibrav))
Crystallography.BravaisLattice(::Val{1}) = BravaisLattice(Cubic(), Primitive())
Crystallography.BravaisLattice(::Val{2}) = BravaisLattice(Cubic(), FaceCentered())
Crystallography.BravaisLattice(::Val{3}) = BravaisLattice(Cubic(), BodyCentered())
Crystallography.BravaisLattice(::Val{4}) = BravaisLattice(Hexagonal(), Primitive())
Crystallography.BravaisLattice(::Val{5}) =
    BravaisLattice(Hexagonal(), RhombohedralCentered())
Crystallography.BravaisLattice(::Val{-5}) =
    BravaisLattice(Hexagonal(), RhombohedralCentered())
Crystallography.BravaisLattice(::Val{6}) = BravaisLattice(Tetragonal(), Primitive())
Crystallography.BravaisLattice(::Val{7}) = BravaisLattice(Tetragonal(), BodyCentered())
Crystallography.BravaisLattice(::Val{8}) = BravaisLattice(Orthorhombic(), Primitive())
Crystallography.BravaisLattice(::Val{9}) = BravaisLattice(Orthorhombic(), BaseCentered(:B))
Crystallography.BravaisLattice(::Val{-9}) = BravaisLattice(Orthorhombic(), BaseCentered(:B))
Crystallography.BravaisLattice(::Val{91}) = BravaisLattice(Orthorhombic(), BaseCentered(:A))  # New in QE 6.5
Crystallography.BravaisLattice(::Val{10}) = BravaisLattice(Orthorhombic(), FaceCentered())
Crystallography.BravaisLattice(::Val{11}) = BravaisLattice(Orthorhombic(), BodyCentered())
Crystallography.BravaisLattice(::Val{12}) = BravaisLattice(Monoclinic(), Primitive())
Crystallography.BravaisLattice(::Val{-12}) = BravaisLattice(Monoclinic(), Primitive())
Crystallography.BravaisLattice(::Val{13}) = BravaisLattice(Monoclinic(), BaseCentered(:C))
Crystallography.BravaisLattice(::Val{-13}) = BravaisLattice(Monoclinic(), BaseCentered(:B))  # New in QE 6.5
Crystallography.BravaisLattice(::Val{14}) = BravaisLattice(Triclinic(), Primitive())

struct AxesSetting{N} end
AxesSetting(N::Int) = AxesSetting{N}()

Lattice(::BravaisLattice{Cubic,Primitive}, p::CellParameters) = Lattice(p[1] * [
    1 0 0
    0 1 0
    0 0 1
])
Lattice(::BravaisLattice{Cubic,FaceCentered}, p::CellParameters) = Lattice(p[1] / 2 * [
    -1 0 1
    0 1 1
    -1 1 0
])
function Lattice(
    ::BravaisLattice{Cubic,BodyCentered},
    p::CellParameters,
    ::AxesSetting{1},
)
    return Lattice(p[1] / 2 * [
        1 1 1
        -1 1 1
        -1 -1 1
    ])
end # function Lattice
function Lattice(
    ::BravaisLattice{Cubic,BodyCentered},
    p::CellParameters,
    ::AxesSetting{2},
)
    return Lattice(p[1] / 2 * [
        -1 1 1
        1 -1 1
        1 1 -1
    ])
end # function Lattice
Lattice(::BravaisLattice{Hexagonal{3},Primitive}, p::CellParameters) =
    Lattice(p[1] * [
        1 0 0
        -1 / 2 √3 / 2 0
        0 0 p[3] / p[1]
    ])
function Lattice(
    ::BravaisLattice{Hexagonal{3},RhombohedralCentered},
    p::CellParameters,
    ::AxesSetting{1},
)
    r = cos(p[4])
    tx = sqrt((1 - r) / 2)
    ty = sqrt((1 - r) / 6)
    tz = sqrt((1 + 2r) / 3)
    return Lattice(p[1] * [
        tx -ty tz
        0 2ty tz
        -tx -ty tz
    ])
end
function Lattice(
    ::BravaisLattice{Hexagonal{3},RhombohedralCentered},
    p::CellParameters,
    ::AxesSetting{2},
)
    ap = p[1] / √3
    γ = acos(p[4])
    ty = sqrt((1 - γ) / 6)
    tz = sqrt((1 + 2γ) / 3)
    u = tz - 2 * √2 * ty
    v = tz + √2 * ty
    return Lattice(ap * [
        u v v
        v u v
        v v u
    ])
end
Lattice(::BravaisLattice{Tetragonal,Primitive}, p::CellParameters) =
    Lattice(p[1] * [
        1 0 0
        0 1 0
        0 0 p[3] / p[1]
    ])
function Lattice(::BravaisLattice{Tetragonal,BodyCentered}, p::CellParameters)
    r = p[3] / p[1]
    return Lattice(p[1] / 2 * [
        1 -1 r
        1 1 r
        -1 -1 r
    ])
end # function Lattice
Lattice(::BravaisLattice{Orthorhombic,Primitive}, p::CellParameters) = Lattice([
    p[1] 0 0
    0 p[2] 0
    0 0 p[3]
])
Lattice(::BravaisLattice{Orthorhombic,BaseCentered{:B}}, p::CellParameters, ::AxesSetting{1}) = Lattice([
    p[1] / 2 p[2] / 2 0
   -p[1] / 2 p[2] / 2 0
    0        0        p[3]
])
Lattice(::BravaisLattice{Orthorhombic,BaseCentered{:B}}, p::CellParameters, ::AxesSetting{2}) = Lattice([
    p[1] / 2 -p[2] / 2 0
    p[1] / 2 p[2] / 2  0
    0        0         p[3]
])
Lattice(::BravaisLattice{Orthorhombic,BaseCentered{:A}}, p::CellParameters) = Lattice([
    p[1] 0 0
    0 p[2] / 2 -p[3] / 2
    0 p[2] / 2 p[3] / 2
])  # New in QE 6.4
function Lattice(::BravaisLattice{Orthorhombic,FaceCentered}, p::CellParameters)
    a, b, c = p.x
    return Lattice([
        a 0 c
        a b 0
        0 b c
    ] / 2)
end # function Lattice
function Lattice(::BravaisLattice{Orthorhombic,BodyCentered}, p::CellParameters)
    a, b, c = p.x
    return Lattice([
        a b c
        -a b c
        -a -b c
    ] / 2)
end
function Lattice(::BravaisLattice{Monoclinic,Primitive}, p::CellParameters, ::AxesSetting{1})
    a, b, c = p.x
    return Lattice([
        a 0 0
        b * cos(p[6]) b * sin(p[6]) 0
        0 0 c
    ])
end
function Lattice(::BravaisLattice{Monoclinic,Primitive}, p::CellParameters, ::AxesSetting{2})
    a, b, c = p.x
    return Lattice([
        a 0 0
        0 b 0
        c * cos(p[5]) 0 c * sin(p[5])
    ])
end
function Lattice(::BravaisLattice{Monoclinic,BaseCentered{:C}}, p::CellParameters)
    a, b, c = p.x
    return Lattice([
        a / 2 0 -c / 2
        b * cos(p[6]) b * sin(p[6]) 0
        a / 2 0 c / 2
    ])
end
function Lattice(::BravaisLattice{Monoclinic,BaseCentered{:B}}, p::CellParameters)
    a, b, c = p.x
    return Lattice([
        a / 2 b / 2 0
        -a / 2 b / 2 0
        c * cos(p[5]) 0 c * sin(p[5])
    ])
end
function Lattice(::BravaisLattice{Monoclinic,BaseCentered{:B}}, p::CellParameters)
    a, b, c = p.x
    α, β, γ = p.y
    zz = c * sqrt(1 + 2cos(α) * cos(β) * cos(γ) - cos(α)^2 - cos(β)^2 - cos(γ)^2) / sin(γ)
    return Lattice([
        a 0 0
        b * cos(γ) b * sin(γ) 0
        c * cos(β) c * (cos(α) - cos(β) * cos(γ)) / sin(γ) zz
    ])
end # function Lattice

function reciprocalof(mat::AbstractMatrix)
    @assert size(mat) == (3, 3)
    volume = abs(det(mat))
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
