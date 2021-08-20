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

struct CellParametersCardSetter <: Setter
    card::CellParametersCard
end

struct AtomicPositionsCardSetter <: Setter
    card::AtomicPositionsCard
end

function (x::CellParametersCardSetter)(template::PWInput)
    if isnothing(template.cell_parameters)
        if optionof(x.card) in ("bohr", "angstrom")
            @set! template.cell_parameters = x.card
            @set! template.system.ibrav = 0
            @set! template.system.celldm = zeros(6)
        else
            @set! template.system.celldm = [template.system.celldm[1]]
            @warn "Please note this `CellParametersCard` might not have the same `alat` as before!"
        end
    else
        if optionof(template.cell_parameters) == "alat"
            if optionof(x.card) in ("bohr", "angstrom")
                @set! template.system.celldm = [template.system.celldm[1]]
                cell_parameters =
                    CellParametersCard(x.card.data / template.system.celldm[1], "alat")
            else  # "alat"
                @warn "Please note this `CellParametersCard` might not have the same `alat` as before!"
            end
        else
            if optionof(cell_parameters) == "alat"
                error("not matched!")
            end
        end
    end
    @set! template.cell_parameters = x.card
    return template
end
function (x::AtomicPositionsCardSetter)(template::PWInput)
    @set! template.atomic_positions = x.card
    return template
end
