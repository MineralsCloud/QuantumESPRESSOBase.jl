module CLI

export PWExec

# See https://stackoverflow.com/a/44446042/3260253
"""
    PWExec(inp; which = "pw.x", nimage = 0, npool = 0, ntg = 0, nyfft = 0, nband = 0, ndiag = 0)

Represent the executable for the PW calculation. Query each field for more information.
"""
mutable struct PWExec
    # docs from https://www.quantum-espresso.org/Doc/user_guide/node18.html
    inp::String
    which::String
    """
    Processors can then be divided into different "images", each corresponding to a
    different self-consistent or linear-response calculation, loosely coupled to others.
    """
    nimage::UInt
    """
    Each image can be subpartitioned into "pools", each taking care of a group of k-points.
    """
    npool::UInt
    """
    In order to allow good parallelization of the 3D FFT when the number of processors
    exceeds the number of FFT planes, FFTs on Kohn-Sham states are redistributed to "task"
    groups so that each group can process several wavefunctions at the same time.
    Alternatively, when this is not possible, a further subdivision of FFT planes is
    performed.
    """
    ntg::UInt
    """
    Task group paralleization and nyfft parallelization, both introduced to improve scaling
    coexist and are in part interchangeable if TG is available it's faster that NYFFT
    because it communicates larger data chuncks less times. But sometimes it is not
    available as for instance when metagga is used or realus or for conjugate gradient.
    `-nyfft` can be used. `-ntg` and `-nyfft` are both allowed flags with the same values.
    These variables are kept separated to help understanding which operation belong to TG or
    to NYFFT. This can enable to make them different if the need arises.
    """
    nyfft::UInt
    """
    Each pool is subpartitioned into "band groups", each taking care of a group of Kohn-Sham
    orbitals (also called bands, or wavefunctions). Especially useful for calculations with
    hybrid functionals.
    """
    nband::UInt
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
    ndiag::UInt
end
PWExec(
    inp;
    which = "pw.x",
    nimage = 0,
    npool = 0,
    ntg = 0,
    nyfft = 0,
    nband = 0,
    ndiag = 0,
) = PWExec(inp, which, nimage, npool, ntg, nyfft, nband, ndiag)

function Base.Cmd(exec::PWExec)
    options = String[]
    for f in Iterators.drop(fieldnames(typeof(exec)), 2)  # 3 to end
        v = getfield(exec, f)
        if !iszero(v)
            push!(options, "-$f", string(v))
        end
    end
    return Cmd([exec.which, options..., "-inp", exec.inp])
end # function Base.Cmd

end
