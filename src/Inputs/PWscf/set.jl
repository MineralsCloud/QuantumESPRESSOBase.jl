export VerbositySetter,
    VolumeSetter, PressureSetter, CellParametersCardSetter, AtomicPositionsCardSetter

function (s::VerbositySetter)(template::PWInput)
    @set! template.control = s(template.control)
    return template
end

function (x::ElectronicTemperatureSetter)(template::PWInput)
    @set! template.system = x(template.system)
    return template
end

struct VolumeSetter{T<:Number} <: Setter
    vol::T
end
function (x::VolumeSetter{<:Real})(template::PWInput)
    factor = cbrt(x.vol / cellvolume(template))
    if isnothing(template.cell_parameters) || optionof(template.cell_parameters) == "alat"
        @set! template.system.celldm[1] *= factor
    else
        @set! template.system.celldm = zeros(6)
        @set! template.cell_parameters =
            optconvert("bohr", CellParametersCard(template.cell_parameters.data * factor))
    end
    return template
end
(x::VolumeSetter{<:AbstractQuantity})(template::PWInput) =
    VolumeSetter(ustrip(u"bohr^3", x.vol))(template)

struct PressureSetter{T<:Number} <: Setter
    press::T
end
function (x::PressureSetter{<:Real})(template::PWInput)
    @set! template.cell.press = x.press
    return template
end
(x::PressureSetter{<:AbstractQuantity})(template::PWInput) =
    PressureSetter(ustrip(u"kbar", x.press))(template)

struct CardSetter{T} <: Setter
    card::T
end

const CellParametersCardSetter = CardSetter{CellParametersCard}
const AtomicPositionsCardSetter = CardSetter{AtomicPositionsCard}

function (x::CellParametersCardSetter)(template::PWInput)
    if optionof(x.card) == "alat"
        if isnothing(template.cell_parameters) ||
           optionof(template.cell_parameters) == "alat"
            @set! template.system.celldm = [template.system.celldm[1]]
        else  # optionof(template.cell_parameters) is "bohr" or "angstrom"
            throw(LackCellInfoError("the `CellParametersCard` does not have units!"))
        end
    else  # "bohr" or "angstrom"
        @set! template.system.celldm = zeros(6)
    end
    @set! template.system.ibrav = 0
    @set! template.cell_parameters = x.card
    return template
end
function (x::AtomicPositionsCardSetter)(template::PWInput)
    @set! template.atomic_positions = x.card
    return template
end
