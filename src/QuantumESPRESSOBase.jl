module QuantumESPRESSOBase

using Crystallography:
    Lattice,
    PrimitiveCubic,
    FaceCenteredCubic,
    BodyCenteredCubic,
    PrimitiveHexagonal,
    RCenteredHexagonal,
    PrimitiveTetragonal,
    BodyCenteredTetragonal,
    PrimitiveOrthorhombic,
    BCenteredOrthorhombic,
    ACenteredOrthorhombic,
    FaceCenteredOrthorhombic,
    BodyCenteredOrthorhombic,
    PrimitiveMonoclinic,
    CCenteredMonoclinic,
    BCenteredMonoclinic,
    PrimitiveTriclinic

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

"""
    Crystallography.Lattice(::Bravais, p[, obverse::Bool])

Create a Bravais lattice from the exact lattice type and cell parameters `p` (not `celldm`!).

The first elements of `p` are `a`, `b`, `c`; the last 3 are `α`, `β`, `γ` (in radians).
"""
Crystallography.Lattice(::PrimitiveCubic, p, args...) = Lattice(p[1] * [
    1 0 0
    0 1 0
    0 0 1
])
Crystallography.Lattice(::FaceCenteredCubic, p, args...) = Lattice(p[1] / 2 * [
    -1 0 1
    0 1 1
    -1 1 0
])
function Crystallography.Lattice(::BodyCenteredCubic, p, obverse::Bool = true)
    if obverse
        return Lattice(p[1] / 2 * [
            1 1 1
            -1 1 1
            -1 -1 1
        ])
    else
        return Lattice(p[1] / 2 * [
            -1 1 1
            1 -1 1
            1 1 -1
        ])
    end
end
Crystallography.Lattice(::PrimitiveHexagonal, p, args...) = Lattice(p[1] * [
    1 0 0
    -1 / 2 √3 / 2 0
    0 0 p[3] / p[1]
])
function Crystallography.Lattice(::RCenteredHexagonal, p, obverse::Bool = true)
    if obverse
        r = cos(p[4])
        tx = sqrt((1 - r) / 2)
        ty = sqrt((1 - r) / 6)
        tz = sqrt((1 + 2r) / 3)
        return Lattice(p[1] * [
            tx -ty tz
            0 2ty tz
            -tx -ty tz
        ])
    else
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
end
Crystallography.Lattice(::PrimitiveTetragonal, p, args...) = Lattice(p[1] * [
    1 0 0
    0 1 0
    0 0 p[3] / p[1]
])
function Crystallography.Lattice(::BodyCenteredTetragonal, p, args...)
    r = p[3] / p[1]
    return Lattice(p[1] / 2 * [
        1 -1 r
        1 1 r
        -1 -1 r
    ])
end
Crystallography.Lattice(::PrimitiveOrthorhombic, p, args...) = Lattice([
    p[1] 0 0
    0 p[2] 0
    0 0 p[3]
])
function Crystallography.Lattice(::BCenteredOrthorhombic, p, obverse::Bool = true)
    a, b, c = Iterators.take(p, 3)
    if obverse
        return Lattice([
            a / 2 b / 2 0
            -a / 2 b / 2 0
            0 0 c
        ])
    else
        return Lattice([
            a / 2 -b / 2 0
            a / 2 b / 2 0
            0 0 c
        ])
    end
end
Crystallography.Lattice(::ACenteredOrthorhombic, p, args...) = Lattice([
    p[1] 0 0
    0 p[2] / 2 -p[3] / 2
    0 p[2] / 2 p[3] / 2
])  # New in QE 6.4
function Crystallography.Lattice(::FaceCenteredOrthorhombic, p, args...)
    a, b, c = Iterators.take(p, 3)
    return Lattice([
        a 0 c
        a b 0
        0 b c
    ] / 2)
end
function Crystallography.Lattice(::BodyCenteredOrthorhombic, p, args...)
    a, b, c = Iterators.take(p, 3)  # Every `p` that is an iterable can be used
    return Lattice([
        a b c
        -a b c
        -a -b c
    ] / 2)
end
function Crystallography.Lattice(::PrimitiveMonoclinic, p, obverse::Bool = true)
    a, b, c, _, β, γ = p
    if obverse
        return Lattice([
            a 0 0
            b * cos(γ) b * sin(γ) 0
            0 0 c
        ])
    else
        return Lattice([
            a 0 0
            0 b 0
            c * cos(β) 0 c * sin(β)
        ])
    end
end
function Crystallography.Lattice(::CCenteredMonoclinic, p, args...)
    a, b, c = Iterators.take(p, 3)
    return Lattice([
        a / 2 0 -c / 2
        b * cos(p[6]) b * sin(p[6]) 0
        a / 2 0 c / 2
    ])
end
function Crystallography.Lattice(::BCenteredMonoclinic, p, args...)
    a, b, c = Iterators.take(p, 3)
    return Lattice([
        a / 2 b / 2 0
        -a / 2 b / 2 0
        c * cos(p[5]) 0 c * sin(p[5])
    ])
end
function Crystallography.Lattice(::PrimitiveTriclinic, p, args...)
    a, b, c, α, β, γ = p  # Every `p` that is an iterable can be used
    δ = c * sqrt(1 + 2 * cos(α) * cos(β) * cos(γ) - cos(α)^2 - cos(β)^2 - cos(γ)^2) / sin(γ)
    return Lattice([
        a 0 0
        b * cos(γ) b * sin(γ) 0
        c * cos(β) c * (cos(α) - cos(β) * cos(γ)) / sin(γ) δ
    ])
end

include("Setters.jl")
include("Inputs/Inputs.jl")
include("CLI.jl")

end # module
