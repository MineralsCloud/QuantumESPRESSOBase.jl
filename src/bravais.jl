using QuantumESPRESSOBase.Namelists.PWscf

export bravais_lattice

"""
    bravais_lattice(nml::SystemNamelist)

Return a 3x3 matrix representing the Bravais lattice from `nml`.
"""
bravais_lattice(nml::SystemNamelist) = bravais_lattice(nml.ibrav, nml.celldm)
"""
    bravais_lattice(ibrav::Integer, celldm::AbstractVector{Union{Missing, Float64}})

Return a 3x3 matrix representing the Bravais lattice from `ibrav` and `celldm`.
"""
bravais_lattice(ibrav::Integer, celldm::AbstractVector{Union{Missing, Float64}}) = bravais_lattice(Val(ibrav), celldm)
bravais_lattice(::Val{1}, celldm::AbstractVector{Union{Missing, Float64}}) = celldm[1] * [
    1 0 0
    0 1 0
    0 0 1
]
bravais_lattice(::Val{2}, celldm::AbstractVector{Union{Missing, Float64}}) = celldm[1] / 2 * [
    -1 0 1
     0 1 1
    -1 1 0
]
bravais_lattice(::Val{3}, celldm::AbstractVector{Union{Missing, Float64}}) = celldm[1] / 2 * [
    1  1 1
   -1  1 1
   -1 -1 1
]
bravais_lattice(::Val{-3}, celldm::AbstractVector{Union{Missing, Float64}}) = celldm[1] / 2 * [
   -1  1  1
    1 -1  1
    1  1 -1
]
bravais_lattice(::Val{4}, celldm::AbstractVector{Union{Missing, Float64}}) = celldm[1] * [
    1  0  0
    -1/2 sqrt(3)/2  0
    0  0 celldm[3]
]
function bravais_lattice(::Val{5}, celldm::AbstractVector{Union{Missing, Float64}})
    c = cellem[3]
    tx = sqrt((1-c)/2)
    ty = sqrt((1-c)/6)
    tz = sqrt((1+2c)/3)
    return
    celldm[1] * [
    tx  -ty  tz
    0  2ty   tz
    -tx  -ty tz
]
end
function bravais_lattice(::Val{-5}, celldm::AbstractVector{Union{Missing, Float64}})
    ap = celldm[1] / sqrt(3)
    c = cellem[3]
    ty = sqrt((1-c)/6)
    tz = sqrt((1+2c)/3)
    u = tz - 2*sqrt(2)*ty
    v = tz + sqrt(2)*ty
    retrun
    ap * [
    u  v  v
    v  u  v
    v  v  u
]
end
bravais_lattice(::Val{6}, celldm::AbstractVector{Union{Missing, Float64}}) = celldm[1] * [
    1  0  0
    0  1  0
    0  0 celldm[3]
]
bravais_lattice(::Val{7}, celldm::AbstractVector{Union{Missing, Float64}}) = celldm[1] / 2 * [
     1   -1  celldm[3]/celldm[1]
     1    1  celldm[3]/celldm[1]
    -1   -1  celldm[3]/celldm[1]
]
bravais_lattice(::Val{8}, celldm::AbstractVector{Union{Missing, Float64}}) = celldm[1] * [
    1      0          0
    0   celldm[2]     0
    0      0       celldm[3]
]
bravais_lattice(::Val{9}, celldm::AbstractVector{Union{Missing, Float64}}) = celldm[1] * [
     1/2      celldm[2]/2          0
    -1/2      celldm[2]/2          0
      0           0            celldm[3]
]
bravais_lattice(::Val{-9}, celldm::AbstractVector{Union{Missing, Float64}}) = celldm[1] * [
     1/2      -celldm[2]/2          0
     1/2       celldm[2]/2          0
      0            0            celldm[3]
]
bravais_lattice(::Val{91}, celldm::AbstractVector{Union{Missing, Float64}}) = celldm[1] * [
     1           0              0
     0     celldm[2]/2    -celldm[3]/2
     0     celldm[2]/2     celldm[3]/2
]
bravais_lattice(::Val{10}, celldm::AbstractVector{Union{Missing, Float64}}) = celldm[1] * [
     1/2           0          celldm[3]/2
     1/2      celldm[2]/2          0
      0       celldm[2]/2     celldm[3]/2
]
bravais_lattice(::Val{11}, celldm::AbstractVector{Union{Missing, Float64}}) = celldm[1] * [
     1/2           0          celldm[3]/2
     1/2      celldm[2]/2          0
      0       celldm[2]/2     celldm[3]/2
]
bravais_lattice(::Val{12}, celldm::AbstractVector{Union{Missing, Float64}}) = celldm[1] * [
              1                              0                       0
     celldm[2]*celldm[4]      celldm[2]*sqrt(1-celldm[4]^2)          0
              0                              0                   celldm[3]
]
bravais_lattice(::Val{-12}, celldm::AbstractVector{Union{Missing, Float64}}) = celldm[1] * [
              1                    0                       0
              0                 celldm[2]                  0
    celldm[3]*celldm[5]            0            celldm[3]*sqrt(1-celldm[5]^2)
]
bravais_lattice(::Val{13}, celldm::AbstractVector{Union{Missing, Float64}}) = celldm[1] * [
           1/2                            0                        -celldm[3]/2
    celldm[2]*celldm[4]       celldm[2]* sqrt(1-celldm[4]^2)             0
          1/2                             0                         celldm[3]/2
]
bravais_lattice(::Val{-13}, celldm::AbstractVector{Union{Missing, Float64}}) = celldm[1] * [
           1/2                 -celldm[2]/2                      0
           1/2                  celldm[2]/2                      0
    celldm[3]*celldm[5]              0             celldm[3]*sqrt(1-celldm[5]^2)
]
bravais_lattice(::Val{14}, celldm::AbstractVector{Union{Missing, Float64}}) = celldm[1] * [
             1                                    0                                                                                  0
    celldm[2]*celldm[6]              celldm[2]* sqrt(1-celldm[6]^2)                                                                   0
    celldm[3]*celldm[5]  celldm[3]*(celldm[4]-celldm[5]*celldm[6])/sqrt(1-celldm[6]^2)   celldm[3] * sqrt(1 + 2*celldm[4]*celldm[5]*celldm[6] - celldm[4]^2 - celldm[5]^2 - celldm[6]^2) / sqrt(1-celldm[6]^2)
]
