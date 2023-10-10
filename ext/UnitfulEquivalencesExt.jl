module UnitfulEquivalencesExt

using Unitful: Temperature, Frequency, Wavenumber, ustrip, @u_str
using UnitfulEquivalences: Thermal, Spectral

import QuantumESPRESSOBase.PWscf: degauss

degauss(value::Temperature) = ustrip(u"Ry", value, Thermal())
degauss(value::Frequency) = ustrip(u"Ry", value, Spectral())
degauss(value::Wavenumber) = ustrip(u"Ry", value, Spectral())

end
