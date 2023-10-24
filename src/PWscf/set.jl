export VerbositySetter,
    VolumeSetter,
    PressureSetter,
    CardSetter,
    CellParametersCardSetter,
    AtomicPositionsCardSetter

function (s::VerbositySetter)(template::PWInput)
    @reset template.control = s(template.control)
    return template
end

function (x::ElectronicTemperatureSetter)(template::PWInput)
    @reset template.system = x(template.system)
    return template
end

struct VolumeSetter{T<:Number} <: Setter
    vol::T
end
function (x::VolumeSetter{<:Real})(template::PWInput)
    factor = cbrt(x.vol / cellvolume(template))
    if isnothing(template.cell_parameters) || getoption(template.cell_parameters) == :alat
        @reset template.system.celldm[1] *= factor
    else
        @reset template.system.celldm = zeros(6)
        @reset template.cell_parameters = convertoption(
            CellParametersCard(template.cell_parameters.data * factor), :bohr
        )
    end
    return template
end
(x::VolumeSetter{<:AbstractQuantity})(template::PWInput) =
    VolumeSetter(ustrip(u"bohr^3", x.vol))(template)

struct PressureSetter{T<:Number} <: Setter
    press::T
end
function (x::PressureSetter{<:Real})(template::PWInput)
    @reset template.cell.press = x.press
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
    if getoption(x.card) == :alat
        if isnothing(template.cell_parameters) ||
            getoption(template.cell_parameters) == :alat
            @reset template.system.celldm = [template.system.celldm[1]]
        else  # getoption(template.cell_parameters) is `:bohr` or `:angstrom`
            throw(InsufficientInfoError("the `CellParametersCard` does not have units!"))
        end
    else  # "bohr" or "angstrom"
        @reset template.system.celldm = zeros(6)
    end
    @reset template.system.ibrav = 0
    @reset template.cell_parameters = x.card
    return template
end
function (x::AtomicPositionsCardSetter)(template::PWInput)
    @reset template.atomic_positions = x.card
    return template
end
