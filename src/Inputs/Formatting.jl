module Formatting

export delimiter, newline, indent, floatfmt, intfmt

delimiter() = ' '
delimiter(x) = delimiter()

newline() = '\n'
newline(x) = newline()

indent() = ' '^4
indent(x) = indent()

floatfmt() = "%14.9f"
floatfmt(x) = floatfmt()

intfmt() = "%5d"
intfmt(x) = intfmt()

end
