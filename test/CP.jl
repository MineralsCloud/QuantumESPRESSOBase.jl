module CP

using Test

using Accessors
using StructArrays: StructArray

using QuantumESPRESSOBase
using QuantumESPRESSOBase.CP

@testset "Constructing `AtomicVelocity`" begin
    # Data from https://gitlab.com/QEF/q-e/blob/master/CPV/examples/autopilot-example/reference/water.autopilot.out
    x = AtomicVelocity("H", [0.140374E-04, -0.333683E-04, 0.231834E-04])
    @test_throws AssertionError @set x.atom = "hydrogen"
    @test_throws InexactError @set x.velocity = [1im, 2im, 3im]
    @test_throws AssertionError x.velocity = [1, 2, 3, 4]
    @test_throws AssertionError x.velocity = [1, 2]
    x.atom = 'H'  # Setting `atom` with a `Char` still works
    @test x.atom == "H"
    @test_throws ErrorException x.velociti = [1, 2, 3] # Given a wrong field name
    @test x == AtomicVelocity("H", [0.140374E-04, -0.333683E-04, 0.231834E-04])
    y = AtomicVelocity('H')  # Incomplete initialization
    @test_throws UndefRefError y == AtomicVelocity("H")
    as = AtomicSpecies("H", 1.00794000, "H.pbe-rrkjus_psl.1.0.0.UPF")
    v = [0.140374E-04, -0.333683E-04, 0.231834E-04]
    @test AtomicVelocity(as).atom == "H"
    @test AtomicVelocity(as, v) == AtomicVelocity("H", v)
    ap = AtomicPosition('H', [0.500000000, 0.288675130, 1.974192764])
    @test AtomicVelocity(ap).atom == "H"
    @test AtomicVelocity(ap, v) == AtomicVelocity("H", v)
end # testset

@testset "Test constructing `AtomicVelocitiesCard` from `StructArray`" begin
    a = ["H", "O"]
    v = [
        [0.140374E-04, -0.333683E-04, 0.231834E-04],
        [-0.111125E-05, -0.466724E-04, 0.972745E-05],
    ]
    card = AtomicVelocitiesCard(StructArray{AtomicVelocity}((a, v)))
    @test card.data == [
        AtomicVelocity("H", [0.140374E-04, -0.333683E-04, 0.231834E-04]),
        AtomicVelocity("O", [-0.111125E-05, -0.466724E-04, 0.972745E-05]),
    ]
end # testset

@testset "Test `push_atom!`" begin
    @testset "`push_atom!` for `AtomicVelocity`" begin
        v = [AtomicVelocity("H", [0.140374E-04, -0.333683E-04, 0.231834E-04])]
        push_atom!(v, "S", "N")
        @test [x.atom for x in v] == ["H", "S", "N"]
    end # testset
    @testset "" begin
        a = ["H", "O"]
        v = [
            [0.140374E-04, -0.333683E-04, 0.231834E-04],
            [-0.111125E-05, -0.466724E-04, 0.972745E-05],
        ]
        card = AtomicVelocitiesCard(StructArray{AtomicVelocity}((a, v)))
        push_atom!(card, "S", "N")
        @test [x.atom for x in card.data] == ["H", "O", "S", "N"]
    end # testset
end # testset

@testset "Test `append_atom!`" begin
    @testset "`append_atom!` to `AtomicVelocity`" begin
        v = [AtomicVelocity("H", [0.140374E-04, -0.333683E-04, 0.231834E-04])]
        append_atom!(v, ["S", "N"])
        @test [x.atom for x in v] == ["H", "S", "N"]
    end # testset
    @testset "`append_atom!` to `AtomicVelocitiesCard`" begin
        a = ["H", "O"]
        v = [
            [0.140374E-04, -0.333683E-04, 0.231834E-04],
            [-0.111125E-05, -0.466724E-04, 0.972745E-05],
        ]
        card = AtomicVelocitiesCard(StructArray{AtomicVelocity}((a, v)))
        append_atom!(card, ["S", "N"])
        @test [x.atom for x in card.data] == ["H", "O", "S", "N"]
    end # testset
end # testset

@testset "Constructing `RefCellParametersCard`" begin
    option = "bohr"
    data = [
        12 0 0
        0 5 0
        0 0 5
    ]
    init = RefCellParametersCard(option, data)
    @test_throws AssertionError @set init.option = "alat"
    @test_throws AssertionError @set init.option = "crystal" # Allowed options are angstrom, bohr
    @test_throws AssertionError @set init.data = [ # Matrix size should be (3, 3)
        1 2
        3 4
    ]
    @test_throws AssertionError @set init.data = [
        1 2 3 4
        5 6 7 8
        4 3 2 1
        8 7 6 5
    ]
    @test RefCellParametersCard(data).option == "bohr"
end # testset

end # module CP
