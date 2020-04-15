module QuantumESPRESSOBase

using Crystallography

import Crystallography

Crystallography.BravaisLattice(ibrav::Integer) = _BravaisLattice(Val(ibrav))
# These are helper methods and should not be exported!
_BravaisLattice(::Val{N}) where {N} = error("Bravais lattice undefined for `ibrav = $N`!")
_BravaisLattice(::Val{1}) = PrimitiveCubic()
_BravaisLattice(::Val{2}) = FaceCenteredCubic()
_BravaisLattice(::Val{3}) = BodyCenteredCubic()
_BravaisLattice(::Val{4}) = PrimitiveHexagonal()
_BravaisLattice(::Val{5}) = RCenteredHexagonal()
_BravaisLattice(::Val{-5}) = RCenteredHexagonal()
_BravaisLattice(::Val{6}) = PrimitiveTetragonal()
_BravaisLattice(::Val{7}) = BodyCenteredTetragonal()
_BravaisLattice(::Val{8}) = PrimitiveOrthorhombic()
_BravaisLattice(::Val{9}) = BCenteredOrthorhombic()
_BravaisLattice(::Val{-9}) = BCenteredOrthorhombic()
_BravaisLattice(::Val{91}) = ACenteredOrthorhombic()  # New in QE 6.5
_BravaisLattice(::Val{10}) = FaceCenteredOrthorhombic()
_BravaisLattice(::Val{11}) = BodyCenteredOrthorhombic()
_BravaisLattice(::Val{12}) = PrimitiveMonoclinic()
_BravaisLattice(::Val{-12}) = PrimitiveMonoclinic()
_BravaisLattice(::Val{13}) = CCenteredMonoclinic()
_BravaisLattice(::Val{-13}) = BCenteredMonoclinic()  # New in QE 6.5
_BravaisLattice(::Val{14}) = PrimitiveTriclinic()

struct AxesSetting{N} end
AxesSetting(N::Int) = AxesSetting{N}()

Crystallography.Lattice(::PrimitiveCubic, p::CellParameters) = Lattice(p[1] * [
    1 0 0
    0 1 0
    0 0 1
], true)
Crystallography.Lattice(::FaceCenteredCubic, p::CellParameters) = Lattice(p[1] / 2 * [
    -1 0 1
    0 1 1
    -1 1 0
], true)
function Crystallography.Lattice(::BodyCenteredCubic, p::CellParameters, ::AxesSetting{1})
    return Lattice(p[1] / 2 * [
        1 1 1
        -1 1 1
        -1 -1 1
    ], true)
end # function Lattice
function Crystallography.Lattice(::BodyCenteredCubic, p::CellParameters, ::AxesSetting{2})
    return Lattice(p[1] / 2 * [
        -1 1 1
        1 -1 1
        1 1 -1
    ], true)
end # function Lattice
Crystallography.Lattice(::PrimitiveHexagonal, p::CellParameters) = Lattice(p[1] * [
    1 0 0
    -1 / 2 √3 / 2 0
    0 0 p[3] / p[1]
], true)
function Crystallography.Lattice(::RCenteredHexagonal, p::CellParameters, ::AxesSetting{1})
    r = cos(p[4])
    tx = sqrt((1 - r) / 2)
    ty = sqrt((1 - r) / 6)
    tz = sqrt((1 + 2r) / 3)
    return Lattice(p[1] * [
        tx -ty tz
        0 2ty tz
        -tx -ty tz
    ], true)
end
function Crystallography.Lattice(::RCenteredHexagonal, p::CellParameters, ::AxesSetting{2})
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
    ], true)
end
Crystallography.Lattice(::PrimitiveTetragonal, p::CellParameters) = Lattice(p[1] * [
    1 0 0
    0 1 0
    0 0 p[3] / p[1]
], true)
function Crystallography.Lattice(::BodyCenteredTetragonal, p::CellParameters)
    r = p[3] / p[1]
    return Lattice(p[1] / 2 * [
        1 -1 r
        1 1 r
        -1 -1 r
    ], true)
end # function Lattice
Crystallography.Lattice(::PrimitiveOrthorhombic, p::CellParameters) = Lattice(
    [
        p[1] 0 0
        0 p[2] 0
        0 0 p[3]
    ],
    true,
)
Crystallography.Lattice(::BCenteredOrthorhombic, p::CellParameters, ::AxesSetting{1}) =
    Lattice(
        [
            p[1] / 2 p[2] / 2 0
            -p[1] / 2 p[2] / 2 0
            0 0 p[3]
        ],
        true,
    )
Crystallography.Lattice(::BCenteredOrthorhombic, p::CellParameters, ::AxesSetting{2}) =
    Lattice(
        [
            p[1] / 2 -p[2] / 2 0
            p[1] / 2 p[2] / 2 0
            0 0 p[3]
        ],
        true,
    )
Crystallography.Lattice(::Tuple{Orthorhombic,BaseCentering{:A}}, p::CellParameters) =
    Lattice(
        [
            p[1] 0 0
            0 p[2] / 2 -p[3] / 2
            0 p[2] / 2 p[3] / 2
        ],
        true,
    )  # New in QE 6.4
function Crystallography.Lattice(::FaceCenteredOrthorhombic, p::CellParameters)
    a, b, c = p.x
    return Lattice([
        a 0 c
        a b 0
        0 b c
    ] / 2, true)
end # function Lattice
function Crystallography.Lattice(::BodyCenteredOrthorhombic, p::CellParameters)
    a, b, c = p.x
    return Lattice([
        a b c
        -a b c
        -a -b c
    ] / 2, true)
end
function Crystallography.Lattice(::PrimitiveMonoclinic, p::CellParameters, ::AxesSetting{1})
    a, b, c = p.x
    return Lattice(
        [
            a 0 0
            b * cos(p[6]) b * sin(p[6]) 0
            0 0 c
        ],
        true,
    )
end
function Crystallography.Lattice(::PrimitiveMonoclinic, p::CellParameters, ::AxesSetting{2})
    a, b, c = p.x
    return Lattice(
        [
            a 0 0
            0 b 0
            c * cos(p[5]) 0 c * sin(p[5])
        ],
        true,
    )
end
function Crystallography.Lattice(::CCenteredMonoclinic, p::CellParameters)
    a, b, c = p.x
    return Lattice(
        [
            a / 2 0 -c / 2
            b * cos(p[6]) b * sin(p[6]) 0
            a / 2 0 c / 2
        ],
        true,
    )
end
function Crystallography.Lattice(::BCenteredMonoclinic, p::CellParameters)
    a, b, c = p.x
    return Lattice(
        [
            a / 2 b / 2 0
            -a / 2 b / 2 0
            c * cos(p[5]) 0 c * sin(p[5])
        ],
        true,
    )
end
function Crystallography.Lattice(::PrimitiveTriclinic, p::CellParameters)
    a, b, c = p.x
    α, β, γ = p.y
    zz =
        c * sqrt(1 + 2 * cos(α) * cos(β) * cos(γ) - cos(α)^2 - cos(β)^2 - cos(γ)^2) / sin(γ)
    return Lattice(
        [
            a 0 0
            b * cos(γ) b * sin(γ) 0
            c * cos(β) c * (cos(α) - cos(β) * cos(γ)) / sin(γ) zz
        ],
        true,
    )
end # function Lattice

include("Setters.jl")
include("Inputs/Inputs.jl")
include("CLI.jl")

end # module
