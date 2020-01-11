using Test

using Parameters: @unpack
using Setfield: set

using QuantumESPRESSOBase.Namelists.PWscf
using QuantumESPRESSOBase.Setters

@testset "Test `VerbositySetter` on `ControlNamelist`" begin
    control = ControlNamelist()
    @test (
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
