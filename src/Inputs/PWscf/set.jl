export VerbositySetter, VolumeSetter, PressureSetter, StructureSetter

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

struct StructureSetter{S,T} <: Setter
    cp::S
    ap::T
end
StructureSetter(cp::CellParametersCard) = StructureSetter(cp, nothing)
StructureSetter(ap::AtomicPositionsCard) = StructureSetter(nothing, ap)
function (x::StructureSetter{CellParametersCard,Nothing})(template::PWInput)
    if isnothing(template.cell_parameters)
        if optionof(x.cp) in ("bohr", "angstrom")
            @set! template.cell_parameters = x.cp
            @set! template.system.ibrav = 0
            @set! template.system.celldm = zeros(6)
        else
            @set! template.system.celldm = [template.system.celldm[1]]
            @warn "Please note this `CellParametersCard` might not have the same `alat` as before!"
        end
    else
        if optionof(template.cell_parameters) == "alat"
            if optionof(x.cp) in ("bohr", "angstrom")
                @set! template.system.celldm = [template.system.celldm[1]]
                cell_parameters =
                    CellParametersCard(x.cp.data / template.system.celldm[1], "alat")
            else  # "alat"
                @warn "Please note this `CellParametersCard` might not have the same `alat` as before!"
            end
        else
            if optionof(cell_parameters) == "alat"
                error("not matched!")
            end
        end
    end
    @set! template.cell_parameters = x.cp
    return template
end
function (x::StructureSetter{Nothing,AtomicPositionsCard})(template::PWInput)
    @set! template.atomic_positions = x.ap
    return template
end
(x::StructureSetter{CellParametersCard,AtomicPositionsCard})(template::PWInput) =
    (StructureSetter(x.ap) ∘ StructureSetter(x.cp))(template)
