module QuantumESPRESSOBase

using Crystallography:
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

import Crystallography: Bravais, Lattice

function Bravais(ibrav::Integer)
    if ibrav == 1
        return PrimitiveCubic()
    elseif ibrav == 2
        return FaceCenteredCubic()
    elseif ibrav == 3
        return BodyCenteredCubic()
    elseif ibrav == 4
        return PrimitiveHexagonal()
    elseif ibrav == 5
        return RCenteredHexagonal()
    elseif ibrav == -5
        return RCenteredHexagonal()
    elseif ibrav == 6
        return PrimitiveTetragonal()
    elseif ibrav == 7
        return BodyCenteredTetragonal()
    elseif ibrav == 8
        return PrimitiveOrthorhombic()
    elseif ibrav == 9
        return BCenteredOrthorhombic()
    elseif ibrav == -9
        return BCenteredOrthorhombic()
    elseif ibrav == 91
        return ACenteredOrthorhombic()  # New in QE 6.5
    elseif ibrav == 10
        return FaceCenteredOrthorhombic()
    elseif ibrav == 11
        return BodyCenteredOrthorhombic()
    elseif ibrav == 12
        return PrimitiveMonoclinic()
    elseif ibrav == -12
        return PrimitiveMonoclinic()
    elseif ibrav == 13
        return CCenteredMonoclinic()
    elseif ibrav == -13
        return BCenteredMonoclinic()  # New in QE 6.5
    elseif ibrav == 14
        return PrimitiveTriclinic()
    else
        error("Bravais lattice undefined for `ibrav = $ibrav`!")
    end
end

"""
    Lattice(::Bravais, p[, obverse::Bool])

Create a Bravais lattice from the exact lattice type and cell parameters `p` (not `celldm`!).

The first elements of `p` are `a`, `b`, `c`; the last 3 are `α`, `β`, `γ` (in radians).
"""
Lattice(::PrimitiveCubic, p, args...) = Lattice(p[1] * [
    1 0 0
    0 1 0
    0 0 1
])
Lattice(::FaceCenteredCubic, p, args...) = Lattice(p[1] / 2 * [
    -1 0 1
    0 1 1
    -1 1 0
])
function Lattice(::BodyCenteredCubic, p, obverse::Bool = true)
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
Lattice(::PrimitiveHexagonal, p, args...) = Lattice(p[1] * [
    1 0 0
    -1 / 2 √3 / 2 0
    0 0 p[3] / p[1]
])
function Lattice(::RCenteredHexagonal, p, obverse::Bool = true)
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
Lattice(::PrimitiveTetragonal, p, args...) = Lattice(p[1] * [
    1 0 0
    0 1 0
    0 0 p[3] / p[1]
])
function Lattice(::BodyCenteredTetragonal, p, args...)
    r = p[3] / p[1]
    return Lattice(p[1] / 2 * [
        1 -1 r
        1 1 r
        -1 -1 r
    ])
end
Lattice(::PrimitiveOrthorhombic, p, args...) = Lattice([
    p[1] 0 0
    0 p[2] 0
    0 0 p[3]
])
function Lattice(::BCenteredOrthorhombic, p, obverse::Bool = true)
    a, b, c = p[1:3]
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
Lattice(::ACenteredOrthorhombic, p, args...) = Lattice([
    p[1] 0 0
    0 p[2] / 2 -p[3] / 2
    0 p[2] / 2 p[3] / 2
])  # New in QE 6.4
function Lattice(::FaceCenteredOrthorhombic, p, args...)
    a, b, c = p[1:3]
    return Lattice([
        a 0 c
        a b 0
        0 b c
    ] / 2)
end
function Lattice(::BodyCenteredOrthorhombic, p, args...)
    a, b, c = p[1:3]
    return Lattice([
        a b c
        -a b c
        -a -b c
    ] / 2)
end
function Lattice(::PrimitiveMonoclinic, p, obverse::Bool = true)
    a, b, c, _, β, γ = p[1:6]
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
function Lattice(::CCenteredMonoclinic, p, args...)
    a, b, c = p[1:3]
    return Lattice([
        a / 2 0 -c / 2
        b * cos(p[6]) b * sin(p[6]) 0
        a / 2 0 c / 2
    ])
end
function Lattice(::BCenteredMonoclinic, p, args...)
    a, b, c = p[1:3]
    return Lattice([
        a / 2 b / 2 0
        -a / 2 b / 2 0
        c * cos(p[5]) 0 c * sin(p[5])
    ])
end
function Lattice(::PrimitiveTriclinic, p, args...)
    a, b, c, α, β, γ = p  # Every `p` that is an iterable can be used
    δ = c * sqrt(1 + 2 * cos(α) * cos(β) * cos(γ) - cos(α)^2 - cos(β)^2 - cos(γ)^2) / sin(γ)
    return Lattice([
        a 0 0
        b * cos(γ) b * sin(γ) 0
        c * cos(β) c * (cos(α) - cos(β) * cos(γ)) / sin(γ) δ
    ])
end

include("Inputs/Inputs.jl")
include("CLI.jl")

end # module
