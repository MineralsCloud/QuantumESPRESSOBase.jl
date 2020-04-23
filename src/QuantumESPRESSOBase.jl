module QuantumESPRESSOBase

using Crystallography

import Crystallography

Crystallography.Bravais(ibrav::Integer) = _Bravais(Val(ibrav))
# These are helper methods and should not be exported!
_Bravais(::Val{N}) where {N} = error("Bravais lattice undefined for `ibrav = $N`!")
_Bravais(::Val{1}) = PrimitiveCubic()
_Bravais(::Val{2}) = FaceCenteredCubic()
_Bravais(::Val{3}) = BodyCenteredCubic()
_Bravais(::Val{4}) = PrimitiveHexagonal()
_Bravais(::Val{5}) = RCenteredHexagonal()
_Bravais(::Val{-5}) = RCenteredHexagonal()
_Bravais(::Val{6}) = PrimitiveTetragonal()
_Bravais(::Val{7}) = BodyCenteredTetragonal()
_Bravais(::Val{8}) = PrimitiveOrthorhombic()
_Bravais(::Val{9}) = BCenteredOrthorhombic()
_Bravais(::Val{-9}) = BCenteredOrthorhombic()
_Bravais(::Val{91}) = ACenteredOrthorhombic()  # New in QE 6.5
_Bravais(::Val{10}) = FaceCenteredOrthorhombic()
_Bravais(::Val{11}) = BodyCenteredOrthorhombic()
_Bravais(::Val{12}) = PrimitiveMonoclinic()
_Bravais(::Val{-12}) = PrimitiveMonoclinic()
_Bravais(::Val{13}) = CCenteredMonoclinic()
_Bravais(::Val{-13}) = BCenteredMonoclinic()  # New in QE 6.5
_Bravais(::Val{14}) = PrimitiveTriclinic()

struct AxesSetting{N} end
AxesSetting(N::Int) = AxesSetting{N}()

Crystallography.Lattice(::PrimitiveCubic, p) = Lattice(p[1] * [
    1 0 0
    0 1 0
    0 0 1
], true)
Crystallography.Lattice(::FaceCenteredCubic, p) = Lattice(p[1] / 2 * [
    -1 0 1
    0 1 1
    -1 1 0
], true)
function Crystallography.Lattice(::BodyCenteredCubic, p, ::AxesSetting{1})
    return Lattice(p[1] / 2 * [
        1 1 1
        -1 1 1
        -1 -1 1
    ], true)
end # function Lattice
function Crystallography.Lattice(::BodyCenteredCubic, p, ::AxesSetting{2})
    return Lattice(p[1] / 2 * [
        -1 1 1
        1 -1 1
        1 1 -1
    ], true)
end # function Lattice
Crystallography.Lattice(::PrimitiveHexagonal, p) = Lattice(p[1] * [
    1 0 0
    -1 / 2 √3 / 2 0
    0 0 p[3] / p[1]
], true)
function Crystallography.Lattice(::RCenteredHexagonal, p, ::AxesSetting{1})
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
function Crystallography.Lattice(::RCenteredHexagonal, p, ::AxesSetting{2})
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
Crystallography.Lattice(::PrimitiveTetragonal, p) = Lattice(p[1] * [
    1 0 0
    0 1 0
    0 0 p[3] / p[1]
], true)
function Crystallography.Lattice(::BodyCenteredTetragonal, p)
    r = p[3] / p[1]
    return Lattice(p[1] / 2 * [
        1 -1 r
        1 1 r
        -1 -1 r
    ], true)
end # function Lattice
Crystallography.Lattice(::PrimitiveOrthorhombic, p) = Lattice(
    [
        p[1] 0 0
        0 p[2] 0
        0 0 p[3]
    ],
    true,
)
Crystallography.Lattice(::BCenteredOrthorhombic, p, ::AxesSetting{1}) =
    Lattice(
        [
            p[1] / 2 p[2] / 2 0
            -p[1] / 2 p[2] / 2 0
            0 0 p[3]
        ],
        true,
    )
Crystallography.Lattice(::BCenteredOrthorhombic, p, ::AxesSetting{2}) =
    Lattice(
        [
            p[1] / 2 -p[2] / 2 0
            p[1] / 2 p[2] / 2 0
            0 0 p[3]
        ],
        true,
    )
Crystallography.Lattice(::Tuple{Orthorhombic,BaseCentering{:A}}, p) =
    Lattice(
        [
            p[1] 0 0
            0 p[2] / 2 -p[3] / 2
            0 p[2] / 2 p[3] / 2
        ],
        true,
    )  # New in QE 6.4
function Crystallography.Lattice(::FaceCenteredOrthorhombic, p)
    a, b, c = p.x
    return Lattice([
        a 0 c
        a b 0
        0 b c
    ] / 2, true)
end # function Lattice
function Crystallography.Lattice(::BodyCenteredOrthorhombic, p)
    a, b, c = p.x
    return Lattice([
        a b c
        -a b c
        -a -b c
    ] / 2, true)
end
function Crystallography.Lattice(::PrimitiveMonoclinic, p, ::AxesSetting{1})
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
function Crystallography.Lattice(::PrimitiveMonoclinic, p, ::AxesSetting{2})
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
function Crystallography.Lattice(::CCenteredMonoclinic, p)
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
function Crystallography.Lattice(::BCenteredMonoclinic, p)
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
function Crystallography.Lattice(::PrimitiveTriclinic, p)
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
