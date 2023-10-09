using StaticArrays: MVector, SVector

@enum Ibrav begin
    PrimitiveCubic = 1
    FaceCenteredCubic = 2
    BodyCenteredCubic = 3
    BodyCenteredCubic2 = -3
    PrimitiveHexagonal = 4
    RCenteredHexagonal = 5
    RCenteredHexagonal2 = -5
    PrimitiveTetragonal = 6
    BodyCenteredTetragonal = 7
    PrimitiveOrthorhombic = 8
    BCenteredOrthorhombic = 9
    BCenteredOrthorhombic2 = -9
    ACenteredOrthorhombic = 91  # New in QE 6.5=91
    FaceCenteredOrthorhombic = 10
    BodyCenteredOrthorhombic = 11
    PrimitiveMonoclinic = 12
    PrimitiveMonoclinic2 = -12
    CCenteredMonoclinic = 13
    BCenteredMonoclinic2 = -13  # New in QE 6.5=-13
    PrimitiveTriclinic = 14
end

latticevectors(p, ibrav::Ibrav) = latticevectors(p, Val(Int(ibrav)))
latticevectors(p, ::Val{1}) = p[1] * [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
latticevectors(p, ::Val{2}) = p[1] / 2 * [[-1, 0, 1], [0, 1, 1], [-1, 1, 0]]
latticevectors(p, ::Val{3}) = p[1] / 2 * [[1, 1, 1], [-1, 1, 1], [-1, -1, 1]]
latticevectors(p, ::Val{-3}) = p[1] / 2 * [[-1, 1, 1], [1, -1, 1], [1, 1, -1]]
latticevectors(p, ::Val{4}) = p[1] * [[1, 0, 0], [-1 / 2, √3 / 2, 0], [0, 0, p[3]]]
function latticevectors(p, ::Val{5})
    cosγ = p[4]
    ty = sqrt((1 - cosγ) / 6)
    tz = sqrt((1 + 2cosγ) / 3)
    tx = sqrt((1 - cosγ) / 2)
    return p[1] * [[tx, -ty, tz], [0, 2ty, tz], [-tx, -ty, tz]]
end
function latticevectors(p, ::Val{-5})
    cosγ = p[4]
    ty = sqrt((1 - cosγ) / 6)
    tz = sqrt((1 + 2cosγ) / 3)
    a′ = p[1] / √3
    u = tz - 2 * √2 * ty
    v = tz + √2 * ty
    return a′ * [[u, v, v], [v, u, v], [v, v, u]]
end
latticevectors(p, ::Val{6}) = p[1] * [[1, 0, 0], [0, 1, 0], [0, 0, p[3]]]
function latticevectors(p, ::Val{7})
    r = p[3]
    return p[1] / 2 * [[1, -1, r], [1, 1, r], [-1, -1, r]]
end
latticevectors(p, ::Val{8}) = p[1] * [[1, 0, 0], [0, p[2], 0], [0, 0, p[3]]]
function latticevectors(p, ::Val{9})
    a, b, c = p[1] .* (1, p[2], p[3])
    return [[a / 2, b / 2, 0], [-a / 2, b / 2, 0], [0, 0, c]]
end
function latticevectors(p, ::Val{-9})
    a, b, c = p[1] .* (1, p[2], p[3])
    return [[a / 2, -b / 2, 0], [a / 2, b / 2, 0], [0, 0, c]]
end
function latticevectors(p, ::Val{91})
    a, r1, r2 = p[1:3]
    return a * [[1, 0, 0], [0, r1 / 2, -r2 / 2], [0, r1 / 2, r2 / 2]]
end  # New in QE 6.4
function latticevectors(p, ::Val{10})
    a, b, c = p[1], p[1] * p[2], p[1] * p[3]
    return [[a, 0, c], [a, b, 0], [0, b, c]] / 2
end
function latticevectors(p, ::Val{11})
    a, b, c = p[1] .* (1, p[2], p[3])
    return [[a, b, c], [-a, b, c], [-a, -b, c]] / 2
end
function latticevectors(p, ::Val{12})
    a, r1, r2, cosγ = p[1:4]
    return a * [[1, 0, 0], [r1 * cosγ, r1 * sin(acos(cosγ)), 0], [0, 0, r2]]
end
function latticevectors(p, ::Val{-12})
    a, r1, r2, _, cosβ = p[1:5]
    return a * [[1, 0, 0], [0, r1, 0], [r2 * cosβ, 0, r2 * sin(acos(cosβ))]]
end
function latticevectors(p, ::Val{13})
    a, r1, r2, cosγ = p[1:4]
    return a *
           [[1 / 2, 0, -r2 / 2], [r1 * cosγ, r1 * sin(acos(cosγ)), 0], [1 / 2, 0, r2 / 2]]
end
function latticevectors(p, ::Val{-13})
    a, r1, r2, _, cosβ = p[1:3]
    return a *
           [[1 / 2, r1 / 2, 0], [-1 / 2, r1 / 2, 0], [r2 * cosβ, 0, r2 * sin(acos(cosβ))]]
end
function latticevectors(p, ::Val{14})
    a, r1, r2, cosα, cosβ, cosγ = p[1:6]  # Every `p` that is an iterable can be used
    sinγ = sin(acos(cosγ))
    δ = r2 * sqrt(1 + 2 * cosα * cosβ * cosγ - cosα^2 - cosβ^2 - cosγ^2) / sinγ
    return a * [
        [1, 0, 0],
        [r1 * cosγ, r1 * sinγ, 0],
        [r2 * cosβ, r2 * (cosα - cosβ * cosγ) / sinγ, δ],
    ]
end

"""
    ReciprocalPoint(x, y, z, w)

Represent a special point of the 3D Brillouin zone. Each of them has a weight `w`.
"""
struct ReciprocalPoint{T}
    coordinates::SVector{3,T}
    weight::Float64
end
ReciprocalPoint(coordinates, weight) =
    ReciprocalPoint(SVector{3,eltype(coordinates)}(coordinates), weight)
ReciprocalPoint(x, y, z, w) = ReciprocalPoint(SVector(x, y, z), w)
