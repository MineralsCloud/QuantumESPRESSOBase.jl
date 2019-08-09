using FilePaths

using QuantumESPRESSOBase.Inputs.PWscf

function Base.write(io::IO, pw::PWscfInput, debug::Bool = true)
    write(io, to_qe(pw, debug = debug))
end
function Base.write(path::AbstractPath, pw::PWscfInput, debug::Bool = true)
    open(path, "r+") do io
        write(io, pw, debug)
    end
end
