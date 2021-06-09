using CrystallographyBase:
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

import CrystallographyBase: Bravais, Lattice

function Bravais(ibrav::Integer)
    if ibrav == 1
        return PrimitiveCubic(true)
    elseif ibrav == 2
        return FaceCenteredCubic(true)
    elseif ibrav == 3
        return BodyCenteredCubic(true)
    elseif ibrav == -3
        return BodyCenteredCubic(false)
    elseif ibrav == 4
        return PrimitiveHexagonal(true)
    elseif ibrav == 5
        return RCenteredHexagonal(true)
    elseif ibrav == -5
        return RCenteredHexagonal(false)
    elseif ibrav == 6
        return PrimitiveTetragonal(true)
    elseif ibrav == 7
        return BodyCenteredTetragonal(true)
    elseif ibrav == 8
        return PrimitiveOrthorhombic(true)
    elseif ibrav == 9
        return BCenteredOrthorhombic(true)
    elseif ibrav == -9
        return BCenteredOrthorhombic(false)
    elseif ibrav == 91
        return ACenteredOrthorhombic(true)  # New in QE 6.5
    elseif ibrav == 10
        return FaceCenteredOrthorhombic(true)
    elseif ibrav == 11
        return BodyCenteredOrthorhombic(true)
    elseif ibrav == 12
        return PrimitiveMonoclinic(true)
    elseif ibrav == -12
        return PrimitiveMonoclinic(false)
    elseif ibrav == 13
        return CCenteredMonoclinic(true)
    elseif ibrav == -13
        return BCenteredMonoclinic(true)  # New in QE 6.5
    elseif ibrav == 14
        return PrimitiveTriclinic(true)
    else
        throw(ArgumentError("Bravais lattice undefined for `ibrav = $ibrav`!"))
    end
end

"""
    Lattice(::Bravais, p)

Create a Bravais lattice from the exact lattice type and cell parameters `p` (not `celldm`!).

The first elements of `p` are `a`, `b`, `c`; the last 3 are `α`, `β`, `γ` (in radians).
"""
Lattice(::PrimitiveCubic, p) = Lattice(p[1] * [
    1 0 0
    0 1 0
    0 0 1
])
Lattice(::FaceCenteredCubic, p) = Lattice(p[1] / 2 * [
    -1 0 1
    0 1 1
    -1 1 0
])
function Lattice(bravais::BodyCenteredCubic, p)
    if bravais.obverse
        return Lattice(p[1] / 2 * [
            1 1 1
            -1 1 1
            -1 -1 1
        ])
    else  # -3
        return Lattice(p[1] / 2 * [
            -1 1 1
            1 -1 1
            1 1 -1
        ])
    end
end
Lattice(::PrimitiveHexagonal, p) = Lattice(p[1] * [
    1 0 0
    -1/2 √3/2 0
    0 0 p[3]
])
function Lattice(bravais::RCenteredHexagonal, p)
    cosγ = p[4]
    ty = sqrt((1 - cosγ) / 6)
    tz = sqrt((1 + 2cosγ) / 3)
    if bravais.obverse
        tx = sqrt((1 - cosγ) / 2)
        return Lattice(p[1] * [
            tx -ty tz
            0 2ty tz
            -tx -ty tz
        ])
    else  # -5
        a′ = p[1] / √3
        u = tz - 2 * √2 * ty
        v = tz + √2 * ty
        return Lattice(a′ * [
            u v v
            v u v
            v v u
        ])
    end
end
Lattice(::PrimitiveTetragonal, p) = Lattice(p[1] * [
    1 0 0
    0 1 0
    0 0 p[3]
])
function Lattice(::BodyCenteredTetragonal, celldm)
    r = celldm[3]
    return Lattice(celldm[1] / 2 * [
        1 -1 r
        1 1 r
        -1 -1 r
    ])
end
Lattice(::PrimitiveOrthorhombic, p) = Lattice(p[1] * [
    1 0 0
    0 p[2] 0
    0 0 p[3]
])
function Lattice(bravais::BCenteredOrthorhombic, p)
    a, b, c = p[1], p[1] * p[2], p[1] * p[3]
    if bravais.obverse
        return Lattice([
            a/2 b/2 0
            -a/2 b/2 0
            0 0 c
        ])
    else
        return Lattice([
            a/2 -b/2 0
            a/2 b/2 0
            0 0 c
        ])
    end
end
function Lattice(::ACenteredOrthorhombic, p)
    a, r1, r2 = p[1:3]
    return Lattice(a * [
        1 0 0
        0 r1/2 -r2/2
        0 r1/2 r2/2
    ])
end  # New in QE 6.4
function Lattice(::FaceCenteredOrthorhombic, p)
    a, b, c = p[1], p[1] * p[2], p[1] * p[3]
    return Lattice([
        a 0 c
        a b 0
        0 b c
    ] / 2)
end
function Lattice(::BodyCenteredOrthorhombic, p)
    a, b, c = p[1], p[1] * p[2], p[1] * p[3]
    return Lattice([
        a b c
        -a b c
        -a -b c
    ] / 2)
end
function Lattice(bravais::PrimitiveMonoclinic, p)
    if bravais.obverse
        a, r1, r2, cosγ = p[1:4]
        return Lattice(a * [
            1 0 0
            r1*cosγ r1*sin(acos(cosγ)) 0
            0 0 r2
        ])
    else
        a, r1, r2, _, cosβ = p[1:5]
        return Lattice(a * [
            1 0 0
            0 r1 0
            r2*cosβ 0 r2*sin(acos(cosβ))
        ])
    end
end
function Lattice(::CCenteredMonoclinic, p)
    a, r1, r2, cosγ = p[1:4]
    return Lattice(a * [
        1/2 0 -r2/2
        r1*cosγ r1*sin(acos(cosγ)) 0
        1/2 0 r2/2
    ])
end
function Lattice(::BCenteredMonoclinic, p)
    a, r1, r2, _, cosβ = p[1:3]
    return Lattice(a * [
        1/2 r1/2 0
        -1/2 r1/2 0
        r2*cosβ 0 r2*sin(acos(cosβ))
    ])
end
function Lattice(::PrimitiveTriclinic, p)
    a, r1, r2, cosα, cosβ, cosγ = p[1:6]  # Every `p` that is an iterable can be used
    sinγ = sin(acos(cosγ))
    δ = r2 * sqrt(1 + 2 * cosα * cosβ * cosγ - cosα^2 - cosβ^2 - cosγ^2) / sinγ
    return Lattice(a * [
        1 0 0
        r1*cosγ r1*sinγ 0
        r2*cosβ r2*(cosα-cosβ*cosγ)/sinγ δ
    ])
end
