using FilePaths

using QuantumESPRESSOBase.Inputs.PWscf

function Base.write(io::IO, pw::PWscfInput, verbose::Bool = true)
    write(io, to_qe(pw, verbose = verbose))
end
function Base.write(path::AbstractPath, pw::PWscfInput, verbose::Bool = true)
    open(path, "r+") do io
        write(io, pw, verbose)
    end
end
