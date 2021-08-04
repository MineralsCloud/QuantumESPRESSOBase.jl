import AbInitioSoftwareBase.Inputs: FormatConfig, asstring

export asstring

FormatConfig(::Type{<:Union{QuantumESPRESSOInput,Namelist}}) = FormatConfig(;
    delimiter = " ",
    newline = "\n",
    indent = ' '^4,
    float = "%f",
    int = "%i",
    bool = ".%.",
)
FormatConfig(x) = FormatConfig(typeof(x))

"""
    asstring(input::QuantumESPRESSOInput)

Return a `String` representing a `QuantumESPRESSOInput`, valid for Quantum ESPRESSO's input.
"""
function asstring(input::QuantumESPRESSOInput)
    newline = FormatConfig(input).newline
    return join(map(_asstring, getfield(input, i) for i in 1:nfields(input)), newline) *
           newline  # Add a new line at the end of line to prevent errors
end
"""
    asstring(nml::Namelist)

Return a `String` representing a `Namelist`, valid for Quantum ESPRESSO's input.
"""
function asstring(nml::Namelist)
    config = FormatConfig(nml)
    dict = dropdefault(nml)
    content = join((_asstring(key, value) for (key, value) in dict), config.newline)
    return join(filter(!isempty, ("&" * groupname(nml), content, '/')), config.newline)
end
function _asstring(key, value::AbstractVector)
    config = FormatConfig(Namelist)
    indent, delimiter, newline = config.indent, config.delimiter, config.newline
    iter = (
        indent * join((string(key, '(', i, ')'), "=", fstring(x)), delimiter) for
        (i, x) in enumerate(value) if !isnothing(x)
    )
    return join(iter, newline)
end
function _asstring(key, value::NamedTuple)
    config = FormatConfig(Namelist)
    indent, delimiter, newline = config.indent, config.delimiter, config.newline
    iter = (
        indent * join((string(key, '%', x), "=", fstring(y)), delimiter) for
        (x, y) in value
    )
    return join(iter, newline)
end
function _asstring(key, value)
    config = FormatConfig(Namelist)
    indent, delimiter = config.indent, config.delimiter
    return indent * join((string(key), "=", fstring(value)), delimiter)
end
_asstring(::Nothing) = ""
_asstring(x) = asstring(x)
