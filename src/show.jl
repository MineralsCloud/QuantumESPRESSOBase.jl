using QuantumESPRESSOBase

function Base.show(io::IO, data::InputEntry)
    print(io, to_qe(data))
end # function Base.show
