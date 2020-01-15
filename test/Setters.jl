module Setters

using Test

using Setfield: set
using Unitful
using UnitfulAtomic

using QuantumESPRESSOBase.Namelists.PWscf
using QuantumESPRESSOBase.Setters

@testset "Test `VerbositySetter` on `ControlNamelist`" begin
    control = ControlNamelist()
    @test (  # Default values
        control.verbosity,
        control.wf_collect,
        control.tstress,
        control.tprnfor,
        control.disk_io,
    ) == ("low", true, false, false, "low")
    control = set(control, VerbositySetter(:high))
    @test (
        control.verbosity,
        control.wf_collect,
        control.tstress,
        control.tprnfor,
        control.disk_io,
    ) == ("high", true, true, true, "high")
    control = set(control, VerbositySetter(:low))
    @test (
        control.verbosity,
        control.wf_collect,
        control.tstress,
        control.tprnfor,
        control.disk_io,
    ) == ("low", false, false, false, "low")
end # testset

@testset "Test `FiniteTemperatureSetter` on `SystemNamelist`" begin
    system = SystemNamelist()
    @test (system.occupations, system.degauss, system.smearing) ==
          ("fixed", 0.0, "gaussian")
    system = set(system, FiniteTemperatureSetter(0.0019000869380733452))
    @test (system.occupations, system.degauss, system.smearing) ==
          ("smearing", 0.0019000869380733452, "fermi-dirac")
    system = set(system, FiniteTemperatureSetter(3e8 * u"μK"))  # 300K
    @test (system.occupations, system.degauss, system.smearing) ==
          ("smearing", 0.0019000869380733452, "fermi-dirac")
    system = set(system, FiniteTemperatureSetter(2.5851999786e-8 * u"MeV"))  # 300K
    @test (system.occupations, system.degauss, system.smearing) ==
          ("smearing", 0.001900086938034304, "fermi-dirac")
    system = set(system, FiniteTemperatureSetter(0.00095004346903668 * u"hartree"))  # 300K
    @test (system.occupations, system.degauss, system.smearing) ==
          ("smearing", 0.00190008693807336, "fermi-dirac")
    system = set(system, FiniteTemperatureSetter(6.250985736 * u"THz"))  # 300K
    @test (system.occupations, system.degauss, system.smearing) ==
          ("smearing", 0.0019000869377698247, "fermi-dirac")
    system = set(system, FiniteTemperatureSetter(20851.044012 * u"m^-1"))  # 300K
    @test (system.occupations, system.degauss, system.smearing) ==
          ("smearing", 0.001900086937837879, "fermi-dirac")
    system = set(system, FiniteTemperatureSetter(208.51044012 * u"cm^-1"))  # 300K
    @test (system.occupations, system.degauss, system.smearing) ==
          ("smearing", 0.001900086937837879, "fermi-dirac")
    system = set(system, FiniteTemperatureSetter(4.141947e-24 * u"kJ"))  # 300K
    @test (system.occupations, system.degauss, system.smearing) ==
          ("smearing", 0.0019000869380663146, "fermi-dirac")
    system = set(system, FiniteTemperatureSetter(4.6085375610000005e-35 * u"g"))  # 300K
    @test (system.occupations, system.degauss, system.smearing) ==
          ("smearing", 0.0019000869377760491, "fermi-dirac")
    system = set(system, FiniteTemperatureSetter(4.6085375610000005e-32 * u"mg"))  # 300K
    @test (system.occupations, system.degauss, system.smearing) ==
          ("smearing", 0.001900086937776049, "fermi-dirac")
end # testset

end
