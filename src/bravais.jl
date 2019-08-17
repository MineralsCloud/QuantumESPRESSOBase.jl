export bravais_lattice

bravais_lattice(nml::SystemNamelist) = bravais_lattice(nml.ibrav, nml.celldm)
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
