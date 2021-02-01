using Preferences: @load_preference, @set_preferences!, @has_preference

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

productname(::Type{PWExec}) = "pw.x"
productname(::Type{PhExec}) = "ph.x"
productname(::Type{Q2rExec}) = "q2r.x"
productname(::Type{MatdynExec}) = "matdyn.x"

const string_nameof = string ∘ nameof

function setbinpath(T::Type{<:QuantumESPRESSOExec}, path)
    @set_preferences!(string_nameof(T) => path)
end
function binpath(T::Type{<:QuantumESPRESSOExec})
    return @load_preference(string_nameof(T))
end
binpath(x::QuantumESPRESSOExec) = binpath(typeof(x))

# Set default paths of the executables
foreach((PWExec, PhExec, Q2rExec, MatdynExec)) do T
    if !@has_preference(string_nameof(T))
        setbinpath(T, productname(T))
    end
end
