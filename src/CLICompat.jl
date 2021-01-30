# See https://github.com/QEF/q-e/blob/884a6f8/Modules/command_line_options.f90
struct PWExec <: QuantumESPRESSOExec
    bin
    nimage::UInt
    npool::UInt
    ntg::UInt
    nyfft::UInt
    nband::UInt
    ndiag::UInt
end
PWExec(; bin = "pw.x", nimage = 0, npool = 0, ntg = 0, nyfft = 0, nband = 0, ndiag = 0) =
    PWExec(bin, nimage, npool, ntg, nyfft, nband, ndiag)
# See https://www.quantum-espresso.org/Doc/ph_user_guide/node14.html
struct PhExec <: QuantumESPRESSOExec
    bin
    nimage::UInt
    npool::UInt
end
PhExec(; bin = "ph.x", nimage = 0, npool = 0) = PhExec(bin, nimage, npool)

struct Q2rExec <: QuantumESPRESSOExec
    bin
end
Q2rExec(; bin = "q2r.x") = Q2rExec(bin)

struct MatdynExec <: QuantumESPRESSOExec
    bin
end
MatdynExec(; bin = "q2r.x") = MatdynExec(bin)

binpath(x) = x.bin
