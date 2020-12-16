struct PhInput <: QuantumESPRESSOInput
    title_line::String
    inputph::PhNamelist
    q_points::Union{Nothing,QPointsCard}
end # struct PhInput
PhInput(inputph::PhNamelist, qpts::QPointsCard) = PhInput(inputph.prefix, inputph, qpts)
PhInput(inputph::PhNamelist) = PhInput(inputph.prefix, inputph, nothing)
PhInput() = PhInput(PhNamelist().prefix, PhNamelist(), nothing)

struct Q2rInput <: QuantumESPRESSOInput
    input::Q2rNamelist
end # struct Q2rInput
Q2rInput() = Q2rInput(Q2rNamelist())

struct MatdynInput <: QuantumESPRESSOInput
    input::MatdynNamelist
    q_points::Union{Nothing,QPointsCard}
end # struct MatdynInput
MatdynInput(input) = MatdynInput(input, nothing)
MatdynInput() = MatdynInput(MatdynNamelist(), nothing)

struct DynmatInput <: QuantumESPRESSOInput
    input::DynmatNamelist
end # struct DynmatInput
DynmatInput() = DynmatInput(DynmatNamelist())
