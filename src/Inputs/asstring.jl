import AbInitioSoftwareBase.Inputs: FormatConfig, asstring

export asstring

formatconfig(::Type{QuantumESPRESSOInput}) = FormatConfig(;
    delimiter = " ",
    newline = "\n",
    indent = ' '^4,
    float = "%f",
    int = "%i",
    bool = ".%.",
)
formatconfig(::Type{<:Namelist}) = FormatConfig(;
    delimiter = " ",
    newline = "\n",
    indent = ' '^4,
    float = "%f",
    int = "%i",
    bool = ".%.",
)
formatconfig(x) = formatconfig(typeof(x))

"""
    asstring(input::QuantumESPRESSOInput)

Return a `String` representing a `QuantumESPRESSOInput`, valid for Quantum ESPRESSO's input.
"""
function asstring(input::QuantumESPRESSOInput)
    newline = formatconfig(input).newline
    return join(map(_asstring, getfield(input, i) for i in 1:nfields(input)), newline) *
           newline  # Add a new line at the end of line to prevent errors
end
"""
    asstring(nml::Namelist)

Return a `String` representing a `Namelist`, valid for Quantum ESPRESSO's input.
"""
function asstring(nml::Namelist)
    config = formatconfig(nml)
    content = _asstring(dropdefault(nml))
    return join(filter(!isempty, ("&" * groupname(nml), content, '/')), config.newline)
end
function _asstring(dict::AbstractDict)
    config = formatconfig(Namelist)
    return join((_asstring(key, value) for (key, value) in dict), config.newline)
end
function _asstring(key, value::AbstractVector)
    config = formatconfig(Namelist)
    indent, delimiter, newline = config.indent, config.delimiter, config.newline
    return join(
        (
            indent * join((string(key, '(', i, ')'), "=", fstring(x)), delimiter) for
            (i, x) in enumerate(value) if !isnothing(x)
        ),
        newline,
    )
end
function _asstring(key, value::NamedTuple)
    config = formatconfig(Namelist)
    indent, delimiter, newline = config.indent, config.delimiter, config.newline
    return join(
        map(keys(value), values(value)) do x, y
            indent * join((string(key, '%', x), "=", fstring(y)), delimiter)
        end,
        newline,
    )
end
function _asstring(key, value)
    config = formatconfig(Namelist)
    indent, delimiter = config.indent, config.delimiter
    return indent * join((string(key), "=", fstring(value)), delimiter)
end
_asstring(::Nothing) = ""
_asstring(x) = asstring(x)
