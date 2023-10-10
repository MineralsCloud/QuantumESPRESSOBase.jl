module UnitfulExt

using CrystallographyBase: Lattice
using QuantumESPRESSOBase: PWInput
using Unitful: AbstractQuantity, Energy, Length, ustrip, @u_str
using UnitfulAtomic: Ry, bohr

import QuantumESPRESSOBase.PWscf: CellParametersCard, VolumeSetter, PressureSetter, degauss

degauss(value::Energy) = ustrip(u"Ry", value)

CellParametersCard(lattice::Lattice{<:Length}) =
    CellParametersCard(Lattice(map(Base.Fix1(ustrip, u"bohr"), parent(lattice))), "bohr")

(x::VolumeSetter{<:AbstractQuantity})(template::PWInput) =
    VolumeSetter(ustrip(u"bohr^3", x.vol))(template)

(x::PressureSetter{<:AbstractQuantity})(template::PWInput) =
    PressureSetter(ustrip(u"kbar", x.press))(template)

end
