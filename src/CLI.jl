module CLI

using AbInitioSoftwareBase.CLI: MpiExec

export MpiExec, PWCmd, PhCmd, Q2rCmd, MatdynCmd

struct PWCmd
    bin
    nimage::UInt
    npool::UInt
    ntg::UInt
    nyfft::UInt
    nband::UInt
    ndiag::UInt
end
PWCmd(; bin = "pw.x", nimage = 0, npool = 0, ntg = 0, nyfft = 0, nband = 0, ndiag = 0) =
    PWCmd(bin, nimage, npool, ntg, nyfft, nband, ndiag)

struct PhCmd
    bin
end
PhCmd(; bin = "ph.x") = PhCmd(bin)

struct Q2rCmd
    bin
end
Q2rCmd(; bin = "q2r.x") = Q2rCmd(bin)

struct MatdynCmd
    bin
end
MatdynCmd(; bin = "q2r.x") = MatdynCmd(bin)

"""
    (::PWCmd)(; bin = "pw.x", nimage = 0, npool = 0, ntg = 0, nyfft = 0, nband = 0, ndiag = 0, stdin = nothing, stdout = nothing, stderr = nothing)

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
function (x::PWCmd)(;
    stdin = nothing,
    stdout = nothing,
    stderr = nothing,
    asstring = false,
    input_redirect = false,
)
    options = String[]
    for k in (:nimage, :npool, :ntg, :nyfft, :nband, :ndiag)
        v = getfield(x, k)
        if !iszero(v)
            push!(options, "-$k", string(v))
        end
    end
    if asstring
        @warn "using string commands maybe error prone! Use with care!"
        for (f, v) in zip((:stdin, :stdout, :stderr), (stdin, stdout, stderr))
            if v !== nothing
                push!(options, redir[f], "'$v'")
            end
        end
        return join((x.bin, options...), " ")
    else
        if input_redirect
            return pipeline(
                Cmd([x.bin; options]);
                stdin = stdin,
                stdout = stdout,
                stderr = stderr,
            )
        else
            return pipeline(
                Cmd([x.bin, options..., "-inp", "$stdin"]);
                stdout = stdout,
                stderr = stderr,
            )
        end
    end
end
# docs from https://www.quantum-espresso.org/Doc/user_guide/node18.html

function (x::PhCmd)(;
    stdin = nothing,
    stdout = nothing,
    stderr = nothing,
    asstring = false,
    input_redirect = false,
)
    options = String[]
    if asstring
        @warn "using string commands maybe error prone! Use with care!"
        for (f, v) in zip((:stdin, :stdout, :stderr), (stdin, stdout, stderr))
            if v !== nothing
                push!(options, redir[f], "'$v'")
            end
        end
        return join((x.bin, options...), " ")
    else
        if input_redirect
            return pipeline(
                Cmd([x.bin; options]);
                stdin = stdin,
                stdout = stdout,
                stderr = stderr,
            )
        else
            return pipeline(
                Cmd([x.bin, options..., "-inp", "$stdin"]);
                stdout = stdout,
                stderr = stderr,
            )
        end
    end
end

function (x::Q2rCmd)(;
    stdin = nothing,
    stdout = nothing,
    stderr = nothing,
    asstring = false,
    input_redirect = false,
)
    options = String[]
    if asstring
        @warn "using string commands maybe error prone! Use with care!"
        for (f, v) in zip((:stdin, :stdout, :stderr), (stdin, stdout, stderr))
            if v !== nothing
                push!(options, redir[f], "'$v'")
            end
        end
        return join((x.bin, options...), " ")
    else
        if input_redirect
            return pipeline(
                Cmd([x.bin; options]);
                stdin = stdin,
                stdout = stdout,
                stderr = stderr,
            )
        else
            return pipeline(
                Cmd([x.bin, options..., "-inp", "$stdin"]);
                stdout = stdout,
                stderr = stderr,
            )
        end
    end
end

function (x::MatdynCmd)(;
    stdin = nothing,
    stdout = nothing,
    stderr = nothing,
    asstring = false,
    input_redirect = false,
)
    options = String[]
    if asstring
        @warn "using string commands maybe error prone! Use with care!"
        for (f, v) in zip((:stdin, :stdout, :stderr), (stdin, stdout, stderr))
            if v !== nothing
                push!(options, redir[f], "'$v'")
            end
        end
        return join((x.bin, options...), " ")
    else
        if input_redirect
            return pipeline(
                Cmd([x.bin; options]);
                stdin = stdin,
                stdout = stdout,
                stderr = stderr,
            )
        else
            return pipeline(
                Cmd([x.bin, options..., "-inp", "$stdin"]);
                stdout = stdout,
                stderr = stderr,
            )
        end
    end
end

const redir = (stdin = "-inp", stdout = "1>", stderr = "2>")
# See https://www.quantum-espresso.org/Doc/pw_user_guide/node21.html

function Base.:∘(mpi::MpiExec, pw::PWCmd)
    function (;
        stdin = nothing,
        stdout = nothing,
        stderr = nothing,
        asstring = false,
        input_redirect = false,
    )
        args = String[]
        for f in (:host, :hostfile)
            v = getfield(mpi, f)
            if !isempty(v)
                push!(args, "-$f", v)
            end
        end
        for (k, v) in mpi.args
            push!(args, "-$k", string(v))
        end
        push!(args, pw.bin)
        for f in (:nimage, :npool, :ntg, :nyfft, :nband, :ndiag)
            v = getfield(pw, f)
            if !iszero(v)
                push!(args, "-$f", string(v))
            end
        end
        if asstring
            @warn "using string commands maybe error prone! Use with care!"
            for (f, v) in zip((:stdin, :stdout, :stderr), (stdin, stdout, stderr))
                if v !== nothing
                    push!(args, redir[f], "'$v'")
                end
            end
            return join(
                (
                    mpi.bin,
                    "-n",
                    mpi.np,
                    "--mca",
                    "btl_vader_single_copy_mechanism",
                    "none",
                    args...,
                ),
                " ",
            )
        else
            if input_redirect
                return pipeline(
                    Cmd([
                        mpi.bin,
                        "-n",
                        string(mpi.np),
                        "--mca",
                        "btl_vader_single_copy_mechanism",
                        "none",
                        args...,
                    ]);
                    stdin = stdin,
                    stdout = stdout,
                    stderr = stderr,
                )
            else
                return pipeline(
                    Cmd([
                        mpi.bin,
                        "-n",
                        string(mpi.np),
                        args...,
                        "--mca",
                        "btl_vader_single_copy_mechanism",
                        "none",
                        "-inp",
                        "$stdin",
                    ]);
                    stdout = stdout,
                    stderr = stderr,
                )
            end
        end
    end
end
function Base.:∘(mpi::MpiExec, ph::PhCmd)
    function (;
        stdin = nothing,
        stdout = nothing,
        stderr = nothing,
        asstring = false,
        input_redirect = false,
    )
        args = String[]
        for f in (:host, :hostfile)
            v = getfield(mpi, f)
            if !isempty(v)
                push!(args, "-$f", v)
            end
        end
        for (k, v) in mpi.args
            push!(args, "-$k", string(v))
        end
        push!(args, ph.bin)
        if asstring
            @warn "using string commands maybe error prone! Use with care!"
            for (f, v) in zip((:stdin, :stdout, :stderr), (stdin, stdout, stderr))
                if v !== nothing
                    push!(args, redir[f], "'$v'")
                end
            end
            return join(
                (
                    mpi.bin,
                    "-n",
                    mpi.np,
                    "--mca",
                    "btl_vader_single_copy_mechanism",
                    "none",
                    args...,
                ),
                " ",
            )
        else
            if input_redirect
                return pipeline(
                    Cmd([
                        mpi.bin,
                        "-n",
                        string(mpi.np),
                        "--mca",
                        "btl_vader_single_copy_mechanism",
                        "none",
                        args...,
                    ]);
                    stdin = stdin,
                    stdout = stdout,
                    stderr = stderr,
                )
            else
                return pipeline(
                    Cmd([
                        mpi.bin,
                        "-n",
                        string(mpi.np),
                        args...,
                        "--mca",
                        "btl_vader_single_copy_mechanism",
                        "none",
                        "-inp",
                        "$stdin",
                    ]);
                    stdout = stdout,
                    stderr = stderr,
                )
            end
        end
    end
end
function Base.:∘(mpi::MpiExec, q2r::Q2rCmd)
    function (;
        stdin = nothing,
        stdout = nothing,
        stderr = nothing,
        asstring = false,
        input_redirect = false,
    )
        args = String[]
        for f in (:host, :hostfile)
            v = getfield(mpi, f)
            if !isempty(v)
                push!(args, "-$f", v)
            end
        end
        for (k, v) in mpi.args
            push!(args, "-$k", string(v))
        end
        push!(args, q2r.bin)
        if asstring
            @warn "using string commands maybe error prone! Use with care!"
            for (f, v) in zip((:stdin, :stdout, :stderr), (stdin, stdout, stderr))
                if v !== nothing
                    push!(args, redir[f], "'$v'")
                end
            end
            return join(
                (
                    mpi.bin,
                    "-n",
                    mpi.np,
                    "--mca",
                    "btl_vader_single_copy_mechanism",
                    "none",
                    args...,
                ),
                " ",
            )
        else
            if input_redirect
                return pipeline(
                    Cmd([
                        mpi.bin,
                        "-n",
                        string(mpi.np),
                        "--mca",
                        "btl_vader_single_copy_mechanism",
                        "none",
                        args...,
                    ]);
                    stdin = stdin,
                    stdout = stdout,
                    stderr = stderr,
                )
            else
                return pipeline(
                    Cmd([
                        mpi.bin,
                        "-n",
                        string(mpi.np),
                        args...,
                        "--mca",
                        "btl_vader_single_copy_mechanism",
                        "none",
                        "-inp",
                        "$stdin",
                    ]);
                    stdout = stdout,
                    stderr = stderr,
                )
            end
        end
    end
end
function Base.:∘(mpi::MpiExec, matdyn::MatdynCmd)
    function (;
        stdin = nothing,
        stdout = nothing,
        stderr = nothing,
        asstring = false,
        input_redirect = false,
    )
        args = String[]
        for f in (:host, :hostfile)
            v = getfield(mpi, f)
            if !isempty(v)
                push!(args, "-$f", v)
            end
        end
        for (k, v) in mpi.args
            push!(args, "-$k", string(v))
        end
        push!(args, matdyn.bin)
        if asstring
            @warn "using string commands maybe error prone! Use with care!"
            for (f, v) in zip((:stdin, :stdout, :stderr), (stdin, stdout, stderr))
                if v !== nothing
                    push!(args, redir[f], "'$v'")
                end
            end
            return join(
                (
                    mpi.bin,
                    "-n",
                    mpi.np,
                    "--mca",
                    "btl_vader_single_copy_mechanism",
                    "none",
                    args...,
                ),
                " ",
            )
        else
            if input_redirect
                return pipeline(
                    Cmd([
                        mpi.bin,
                        "-n",
                        string(mpi.np),
                        "--mca",
                        "btl_vader_single_copy_mechanism",
                        "none",
                        args...,
                    ]);
                    stdin = stdin,
                    stdout = stdout,
                    stderr = stderr,
                )
            else
                return pipeline(
                    Cmd([
                        mpi.bin,
                        "-n",
                        string(mpi.np),
                        args...,
                        "--mca",
                        "btl_vader_single_copy_mechanism",
                        "none",
                        "-inp",
                        "$stdin",
                    ]);
                    stdout = stdout,
                    stderr = stderr,
                )
            end
        end
    end
end

end
