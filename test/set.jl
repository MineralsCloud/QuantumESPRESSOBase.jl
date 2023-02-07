module Setters

using Test: @testset, @test
using Unitful
using UnitfulAtomic

using QuantumESPRESSOBase.Inputs.PWscf:
    ControlNamelist, SystemNamelist, VerbositySetter, ElectronicTemperatureSetter

@testset "Apply `VerbositySetter` on `ControlNamelist`" begin
    control = ControlNamelist()
    @test (  # Default values
        control.verbosity,
        control.wf_collect,
        control.tstress,
        control.tprnfor,
        control.disk_io,
    ) == ("low", true, false, false, "low")
    @testset "Set verbosity to 'high'" begin
        result = VerbositySetter("high")(control)
        @test (
            result.verbosity,
            result.wf_collect,
            result.tstress,
            result.tprnfor,
            result.disk_io,
        ) == ("high", true, true, true, "high")
    end
    @testset "Set verbosity to 'low'" begin
        result = VerbositySetter("low")(control)
        @test (
            result.verbosity,
            result.wf_collect,
            result.tstress,
            result.tprnfor,
            result.disk_io,
        ) == ("low", false, false, false, "low")
    end
end

@testset "Apply `ElectronicTemperatureSetter` on `SystemNamelist`" begin
    system = SystemNamelist()
    @test (system.occupations, system.degauss, system.smearing) ==
        ("fixed", 0.0, "gaussian")
    for value in (
        0.0019000869380733452,
        300u"K",
        3e8u"μK",
        2.5851999786e-8u"MeV",
        0.00095004346903668u"hartree",
        6.250985736u"THz",
        20851.044012u"m^-1",
        208.51044012u"cm^-1",
        4.141947e-24u"kJ",
    )
        result = ElectronicTemperatureSetter(value)(system)
        @test result.occupations == "smearing"
        @test result.smearing == "fermi-dirac"
        @test result.degauss ≈ 0.0019000869380733452
    end
end

end
