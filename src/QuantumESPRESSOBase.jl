module QuantumESPRESSOBase

using Compat: isnothing, eachrow
using Fortran90Namelists.JuliaToFortran: to_fortran
using Parameters: @with_kw

export PWCommand
export asfieldname, titleof, to_qe

abstract type InputEntry end

"""
    asfieldname()



# Arguments

# Examples

```jldoctest
julia>
```
"""
function asfieldname end

"""
    titleof(::Type{<:InputEntry})

Return the title of the input entry in Quantum ESPRESSO.

The definition `titleof(x) = titleof(typeof(x))` is provided for convenience so that
instances can be passed instead of types.

# Examples

```jldoctest
julia> titleof(SystemNamelist())
"SYSTEM"

julia> titleof(SystemNamelist)
"SYSTEM"
```
"""
titleof(x::InputEntry) = titleof(typeof(x))

"""
    to_qe(x, indent::AbstractString = "    ", sep::AbstractString = " ")

Return a string representing the object, valid form Quantum ESPRESSO's input.
"""
function to_qe(
    dict::AbstractDict;
    indent::AbstractString = "    ",
    sep::AbstractString = " ",
)::String
    content = ""
    f = string ∘ to_fortran
    for (key, value) in dict
        if value isa Vector
            for (i, x) in enumerate(value)
                isnothing(x) && continue
                content *= indent * join(["$key($i)", "=", "$(f(x))\n"], sep)
            end
        else
            content *= indent * join(["$key", "=", "$(f(value))\n"], sep)
        end
    end
    return content
end

"""
    PWCommand(exe = "pw.x", inp, nimage = 0, npool = 0, ntg = 0, nyfft = 0, nband = 0, ndiag = 0)

Represent the executable for the PW calculation. Query each field for more information.
"""
@with_kw struct PWCommand
    # docs from https://www.quantum-espresso.org/Doc/user_guide/node18.html
    exe::String = "pw.x"
    inp::String
    """
    Processors can then be divided into different "images", each corresponding to a
    different self-consistent or linear-response calculation, loosely coupled to others.
    """
    nimage::Int = 0
    """
    Each image can be subpartitioned into "pools", each taking care of a group of k-points.
    """
    npool::Int = 0
    """
    In order to allow good parallelization of the 3D FFT when the number of processors
    exceeds the number of FFT planes, FFTs on Kohn-Sham states are redistributed to "task"
    groups so that each group can process several wavefunctions at the same time.
    Alternatively, when this is not possible, a further subdivision of FFT planes is
    performed.
    """
    ntg::Int = 0
    """
    Task group paralleization and nyfft parallelization, both introduced to improve scaling
    coexist and are in part interchangeable if TG is available it's faster that NYFFT
    because it communicates larger data chuncks less times. But sometimes it is not
    available as for instance when metagga is used or realus or for conjugate gradient.
    `-nyfft` can be used. `-ntg` and `-nyfft` are both allowed flags with the same values.
    These variables are kept separated to help understanding which operation belong to TG or
    to NYFFT. This can enable to make them different if the need arises.
    """
    nyfft::Int = ntg
    """
    Each pool is subpartitioned into "band groups", each taking care of a group of Kohn-Sham
    orbitals (also called bands, or wavefunctions). Especially useful for calculations with
    hybrid functionals.
    """
    nband::Int = 0
    """
    A further level of parallelization, independent on PW or k-point parallelization, is the
    parallelization of subspace diagonalization / iterative orthonormalization. Both
    operations required the diagonalization of arrays whose dimension is the number of
    Kohn-Sham states (or a small multiple of it). All such arrays are distributed block-like
    across the "linear-algebra group", a subgroup of the pool of processors, organized in a
    square 2D grid. As a consequence the number of processors in the linear-algebra group is
    given by ``n^2``, where n is an integer; ``n^2`` must be smaller than the number of
    processors in the PW group. The diagonalization is then performed in parallel using
    standard linear algebra operations. (This diagonalization is used by, but should not be
    confused with, the iterative Davidson algorithm). The preferred option is to use ELPA
    and ScaLAPACK; alternative built-in algorithms are anyway available.
    """
    ndiag::Int = 0
end

include("Namelists/Namelists.jl")
include("Cards/Cards.jl")
include("bravais_lattice.jl")
include("Inputs/Inputs.jl")

end # module
