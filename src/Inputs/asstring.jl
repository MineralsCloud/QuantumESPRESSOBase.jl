import AbInitioSoftwareBase.Inputs: asstring

export asstring

"""
    asstring(input::QuantumESPRESSOInput)

Return a `String` representing a `QuantumESPRESSOInput`, valid for Quantum ESPRESSO's input.
"""
function asstring(input::QuantumESPRESSOInput)
    return join(
        map(
            asstring,
            Iterators.filter(!isnothing, getfield(input, i) for i in 1:nfields(input)),
        ),
        newline(input),
    ) * newline(input)  # Add a new line at the end of line to prevent errors
end
"""
    asstring(nml::Namelist)

Return a `String` representing a `Namelist`, valid for Quantum ESPRESSO's input.
"""
function asstring(nml::Namelist)
    content = _nmlasstring(
        dropdefault(nml);
        indent = indent(nml),
        delimiter = delimiter(nml),
        newline = newline(nml),
    )
    return join(filter(!isempty, ("&" * groupname(nml), content, '/')), newline(nml))
end
asstring(x::AbstractString) = string(x)
function _nmlasstring(dict::AbstractDict; indent = ' '^4, delimiter = ' ', newline = '\n')
    return join(
        map(keys(dict), values(dict)) do key, value
            _nmlasstring(
                key,
                value;
                indent = indent,
                delimiter = delimiter,
                newline = newline,
            )
        end,
        newline,
    )
end
function _nmlasstring(
    key,
    value::AbstractVector;
    indent = ' '^4,
    delimiter = ' ',
    newline = '\n',
)
    return join(
        (
            indent * join((string(key, '(', i, ')'), "=", fstring(x)), delimiter) for
            (i, x) in enumerate(value) if !isnothing(x)
        ),
        newline,
    )
end
function _nmlasstring(
    key,
    value::NamedTuple;
    indent = ' '^4,
    delimiter = ' ',
    newline = '\n',
)
    return join(
        map(keys(value), values(value)) do x, y
            indent * join((string(key, '%', x), "=", fstring(y)), delimiter)
        end,
        newline,
    )
end
_nmlasstring(key, value; indent = ' '^4, delimiter = ' ', newline = '\n') =
    indent * join((string(key), "=", fstring(value)), delimiter)

newline(::Union{QuantumESPRESSOInput,Namelist,Card}) = '\n'

indent(::Namelist) = ' '^4

delimiter(::Namelist) = ' '
