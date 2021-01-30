module CLI

using AbInitioSoftwareBase.CLI: Executable, Mpiexec
using Preferences: @load_preference, @set_preferences!

import AbInitioSoftwareBase.CLI: scriptify

export PWExec, PhExec, Q2rExec, MatdynExec, setbinpath, binpath, scriptify

# const REDIRECTION_OPERATORS = ("-inp", "1>", "2>")
# See https://www.quantum-espresso.org/Doc/pw_user_guide/node21.html 5.0.0.3

abstract type QuantumESPRESSOExec <: Executable end
# See https://github.com/QEF/q-e/blob/884a6f8/Modules/command_line_options.f90
struct PWExec <: QuantumESPRESSOExec
    nimage::UInt
    npool::UInt
    ntg::UInt
    nyfft::UInt
    nband::UInt
    ndiag::UInt
end
PWExec(; nimage = 0, npool = 0, ntg = 0, nyfft = 0, nband = 0, ndiag = 0) =
    PWExec(nimage, npool, ntg, nyfft, nband, ndiag)
# See https://www.quantum-espresso.org/Doc/ph_user_guide/node14.html
struct PhExec <: QuantumESPRESSOExec
    nimage::UInt
    npool::UInt
end
PhExec(; nimage = 0, npool = 0) = PhExec(nimage, npool)

struct Q2rExec <: QuantumESPRESSOExec end

struct MatdynExec <: QuantumESPRESSOExec end

function _prescriptify(  # Never export!
    x::QuantumESPRESSOExec,
    stdin,
    stdout,
    stderr,
    use_shell,
    input_not_read,
)
    args = [binpath(typeof(x))]
    if x isa PWExec
        for k in (:nimage, :npool, :ntg, :nyfft, :nband, :ndiag)
            v = getfield(x, k)
            if !iszero(v)
                push!(args, "-$k", string(v))
            end
        end
    end
    if use_shell
        @warn "using shell maybe error prone!"
        if input_not_read
            REDIRECTION_OPERATORS = ("-inp", "1>", "2>")
        else
            REDIRECTION_OPERATORS = ("<", "1>", "2>")
        end
        for (i, v) in enumerate((stdin, stdout, stderr))
            if v !== nothing
                push!(args, REDIRECTION_OPERATORS[i], "'$v'")
            end
        end
    else
        if input_not_read
            push!(args, "-inp", "$stdin")
        end
    end
    return args
end
"""
    (::PWX)(; bin = "pw.x", nimage = 0, npool = 0, ntg = 0, nyfft = 0, nband = 0, ndiag = 0, stdin = nothing, stdout = nothing, stderr = nothing)

# Arguments
- `bin`: the path to the PWscf executable, usually is `\"pw.x\"`. Better be an absolute path.
- `nimage = 0`: processors can then be divided into different "images", each corresponding to a different self-consistent or linear-response calculation, loosely coupled to others.
- `npool = 0`: each image can be subpartitioned into "pools", each taking care of a group of k-points.
- `ntg = 0`: in order to allow good parallelization of the 3D FFT when the number of processors exceeds the number of FFT planes, FFTs on Kohn-Sham states are redistributed to "task" groups so that each group can process several wavefunctions at the same time. Alternatively, when this is not possible, a further subdivision of FFT planes is performed.
- `nyfft = 0`: task group paralleization and nyfft parallelization, both introduced to improve scaling coexist and are in part interchangeable if TG is available it's faster that NYFFT because it communicates larger data chuncks less times. But sometimes it is not available as for instance when metagga is used or realus or for conjugate gradient. `-nyfft` can be used. `-ntg` and `-nyfft` are both allowed flags with the same values. These variables are kept separated to help understanding which operation belong to TG or to NYFFT. This can enable to make them different if the need arises.
- `nband = 0`: each pool is subpartitioned into "band groups", each taking care of a group of Kohn-Sham orbitals (also called bands, or wavefunctions). Especially useful for calculations with hybrid functionals.
- `ndiag = 0`: a further level of parallelization, independent on PW or k-point parallelization, is the parallelization of subspace diagonalization / iterative orthonormalization. Both operations required the diagonalization of arrays whose dimension is the number of Kohn-Sham states (or a small multiple of it). All such arrays are distributed block-like across the "linear-algebra group", a subgroup of the pool of processors, organized in a square 2D grid. As a consequence the number of processors in the linear-algebra group is given by ``n^2``, where n is an integer; ``n^2`` must be smaller than the number of processors in the PW group. The diagonalization is then performed in parallel using standard linear algebra operations. (This diagonalization is used by, but should not be confused with, the iterative Davidson algorithm). The preferred option is to use ELPA and ScaLAPACK; alternative built-in algorithms are anyway available.
- `stdin = nothing`: input
- `stdout = nothing`: output
- `stderr = nothing`: error
"""
function scriptify(
    x::QuantumESPRESSOExec;
    stdin,
    stdout = nothing,
    stderr = nothing,
    dir = dirname(stdin),  # If `stdin` path is not complete, this will save it
    use_shell = false,
    input_not_read = false,
)
    args = _prescriptify(x, stdin, stdout, stderr, use_shell, input_not_read)
    return _postscriptify(args, stdin, stdout, stderr, dir, use_shell, input_not_read)
end
# docs from https://www.quantum-espresso.org/Doc/user_guide/node18.html
function scriptify(
    mpi::Mpiexec,
    x::QuantumESPRESSOExec;
    stdin,
    stdout = nothing,
    stderr = nothing,
    dir = dirname(stdin),
    use_shell = false,
    input_not_read = false,
)
    cmd = [mpi.bin, "-n", string(mpi.np)]
    for (k, v) in mpi.args
        push!(cmd, "$k", string(v))
    end
    args = _prescriptify(x, stdin, stdout, stderr, use_shell, input_not_read)
    append!(cmd, args)
    return _postscriptify(cmd, stdin, stdout, stderr, dir, use_shell, input_not_read)
end
function _postscriptify(args, stdin, stdout, stderr, dir, use_shell, input_not_read)
    if use_shell
        mkpath(dir)
        path, _ = mktemp(dir; cleanup = false)
        open(path, "w") do io
            shell = join(args, " ")
            write(io, shell)
        end
        chmod(path, 0o755)
        return setenv(Cmd([abspath(path)]); dir = dir)
    else
        cmd = pipeline(setenv(Cmd(args); dir = dir), stdout = stdout, stderr = stderr)
        if input_not_read
            return cmd
        else
            return pipeline(cmd; stdin = stdin)
        end
    end
end

productname(::Type{PWExec}) = "pw.x"
productname(::Type{PhExec}) = "ph.x"
productname(::Type{Q2rExec}) = "q2r.x"
productname(::Type{MatdynExec}) = "matdyn.x"

function setbinpath(T::Type{<:QuantumESPRESSOExec}, path)
    @set_preferences!(productname(T) => path)
end
function binpath(T::Type{<:QuantumESPRESSOExec})
    return @load_preference(productname(T))
end

# Set default paths of the executables
setbinpath(PWExec, "pw.x")
setbinpath(PhExec, "ph.x")
setbinpath(Q2rExec, "q2r.x")
setbinpath(MatdynExec, "matdyn.x")

end
