using FilePaths

using QuantumESPRESSOBase.Inputs.PWscf

function Base.write(io::IO, pw::PWInput, debug::Bool = true)
    write(io, to_qe(pw, debug = debug))
end
function Base.write(path::AbstractPath, pw::PWInput, debug::Bool = true)
    open(path, "r+") do io
        write(io, pw, debug)
    end
end
