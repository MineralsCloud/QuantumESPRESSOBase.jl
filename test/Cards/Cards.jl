using Test

using Setfield
using StructArrays: StructArray

using QuantumESPRESSOBase
using QuantumESPRESSOBase.Cards.PWscf

@testset "Constructing `AtomicSpecies`" begin
    # Data from https://github.com/QEF/q-e/blob/7be27df/PW/examples/gatefield/run_example#L128.
    x = AtomicSpecies("S", 32.066, "S.pz-n-rrkjus_psl.0.1.UPF")
    @test_throws AssertionError @set x.atom = "sulfur"
    @test_throws InexactError @set x.mass = 1im
    @test x == AtomicSpecies('S', 32.066, "S.pz-n-rrkjus_psl.0.1.UPF")
    y = AtomicSpecies('S')  # Incomplete initialization
    @test_throws UndefRefError y == AtomicSpecies("S")
    @test_throws UndefRefError y.pseudopot
    @test_throws AssertionError y.atom = "sulfur"
    y.mass, y.pseudopot = 32.066, "S.pz-n-rrkjus_psl.0.1.UPF"
    @test x == y  # Constructing `AtomicSpecies` in 3 steps is equivalent to a one-time construction
    @test AtomicSpecies(AtomicPosition(
        'S',
        [0.500000000, 0.288675130, 1.974192764],
    )).atom == "S"
    @test AtomicSpecies(
        AtomicPosition('S', [0.500000000, 0.288675130, 1.974192764]),
        32.066,
        "S.pz-n-rrkjus_psl.0.1.UPF",
    ) == x
end # testset

@testset "Test constructing `AtomicSpeciesCard` from `StructArray`s" begin
    atoms = ["Al", "As"]
    masses = [24590.7655930491, 68285.4024548272]
    pseudopotentials = ["Al.pbe-n-kjpaw_psl.1.0.0.UPF", "As.pbe-n-kjpaw_psl.1.0.0.UPF"]
    init = AtomicSpeciesCard(StructArray{AtomicSpecies}((atoms, masses, pseudopotentials)))
    @test init.data == [
        AtomicSpecies("Al", 24590.7655930491, "Al.pbe-n-kjpaw_psl.1.0.0.UPF"),
        AtomicSpecies("As", 68285.4024548272, "As.pbe-n-kjpaw_psl.1.0.0.UPF"),
    ]
    push!(init.data, AtomicSpecies("Si", 25591.1924913552, "Si.pbe-n-kjpaw_psl.1.0.0.UPF"))
    @test init.data == [
        AtomicSpecies("Al", 24590.7655930491, "Al.pbe-n-kjpaw_psl.1.0.0.UPF"),
        AtomicSpecies("As", 68285.4024548272, "As.pbe-n-kjpaw_psl.1.0.0.UPF"),
        AtomicSpecies("Si", 25591.1924913552, "Si.pbe-n-kjpaw_psl.1.0.0.UPF"),
    ]
    @testset "Mutual construction" begin
        @test map(x -> x.atom, AtomicPositionsCard("alat", init).data) == ["Al", "As", "Si"]
    end # testset
    @testset "Test `pseudopot_format`" begin
        @test unique(pseudopot_format.(init.data)) == [UnifiedPseudopotentialFormat()]
    end # testset
    @testset "Test `to_qe`" begin
        @test to_qe(init) ==
              "ATOMIC_SPECIES\n     Al     24590.7655930491 Al.pbe-n-kjpaw_psl.1.0.0.UPF\n     As     68285.4024548272 As.pbe-n-kjpaw_psl.1.0.0.UPF\n     Si     25591.1924913552 Si.pbe-n-kjpaw_psl.1.0.0.UPF"
        @test to_qe(init; delim = "", indent = "") ==
              "ATOMIC_SPECIES\n Al    24590.7655930491Al.pbe-n-kjpaw_psl.1.0.0.UPF\n As    68285.4024548272As.pbe-n-kjpaw_psl.1.0.0.UPF\n Si    25591.1924913552Si.pbe-n-kjpaw_psl.1.0.0.UPF"
    end # testset
end # testset

@testset "Constructing `AtomicPosition`" begin
    # Data from https://github.com/QEF/q-e/blob/7be27df/PW/examples/gatefield/run_example#L129-L132.
    x = AtomicPosition("S", [0.500000000, 0.288675130, 1.974192764])
    @test_throws AssertionError @set x.atom = "sulfur"
    @test_throws InexactError @set x.pos = [1im, 2im, 3im]
    @test_throws AssertionError x.pos = [1, 2, 3, 4]
    @test_throws AssertionError x.if_pos = [1, 2, 3]
    @test_throws AssertionError x.if_pos = [1, 0, 1, 1]
    @test x.if_pos == [1, 1, 1]
    @test x == AtomicPosition('S', [0.500000000, 0.288675130, 1.974192764])
    y = AtomicPosition('S')  # Incomplete initialization
    @test_throws UndefRefError y == AtomicPosition("S")
    @test_throws UndefRefError y.pos
    @test_throws UndefRefError y.if_pos
    @test_throws AssertionError y.atom = "sulfur"
    y.pos, y.if_pos = [0.500000000, 0.288675130, 1.974192764], [1, 1, 1]
    @test x == y  # Constructing `AtomicSpecies` in 3 steps is equivalent to a one-time construction
    @test AtomicPosition(AtomicSpecies('S', 32.066, "S.pz-n-rrkjus_psl.0.1.UPF")).atom ==
          "S"
    @test AtomicPosition(
        AtomicSpecies('S', 32.066, "S.pz-n-rrkjus_psl.0.1.UPF"),
        [0.500000000, 0.288675130, 1.974192764],
    ) == x
end # testset

@testset "Test constructing `AtomicPositionsCard` from `StructArray`s" begin
    # Data from https://github.com/QEF/q-e/blob/7be27df/PW/examples/gatefield/run_example#L129-L132.
    atoms = ["S", "Mo", "S"]
    positions = [
        [0.500000000, 0.288675130, 1.974192764],
        [0.000000000, 0.577350270, 2.462038339],
        [0.000000000, -0.577350270, 2.950837559],
    ]
    init = AtomicPositionsCard(
        "alat",
        StructArray{AtomicPosition}((atoms, positions, [[1, 1, 1], [1, 1, 1], [1, 1, 1]])),
    )
    @test init.data == [
        AtomicPosition("S", [0.500000000, 0.288675130, 1.974192764]),
        AtomicPosition("Mo", [0.000000000, 0.577350270, 2.462038339]),
        AtomicPosition("S", [0.000000000, -0.577350270, 2.950837559]),
    ]
    @testset "Mutual construction" begin
        @test map(x -> x.atom, AtomicSpeciesCard(init).data) == ["S", "Mo", "S"]
    end # testset
    @testset "Test `to_qe`" begin
        @test to_qe(init) ==
              "ATOMIC_POSITIONS { alat }\n      S    0.500000000    0.288675130    1.974192764\n     Mo    0.000000000    0.577350270    2.462038339\n      S    0.000000000   -0.577350270    2.950837559"
        @test to_qe(init; delim = "", indent = "") ==
              "ATOMIC_POSITIONS { alat }\n  S   0.500000000   0.288675130   1.974192764\n Mo   0.000000000   0.577350270   2.462038339\n  S   0.000000000  -0.577350270   2.950837559"
    end # testset
end # testset

@testset "Construct `MonkhorstPackGrid` incorrectly" begin
    @test_throws AssertionError MonkhorstPackGrid([4, 4, 4, 4], [1, 1, 1])
    @test_throws AssertionError MonkhorstPackGrid([4, 4, 0], [1, 1, 1])
    @test_throws AssertionError MonkhorstPackGrid([4, 4, 4], [1, 1, 1, 1])
    @test_throws AssertionError MonkhorstPackGrid([4, 4, 0], [1, 1, 2])
end # testset

@testset "Construct `SpecialKPoint`" begin
    @test_throws AssertionError SpecialKPoint([1 / 4, 1 / 4], 1 / 2)
    @test SpecialKPoint([3 / 4, 1 / 4, 1 / 4], 1 / 2) ==
          SpecialKPoint(3 / 4, 1 / 4, 1 / 4, 1 / 2)
end # testset

@testset "Construct `KPointsCard` incorrectly" begin
    @test_throws AssertionError KPointsCard("automatic", GammaPoint())
    @test_throws AssertionError KPointsCard(
        "automatic",
        [
            SpecialKPoint([3 / 4, 1 / 4, 1 / 4], 1 / 2),
            SpecialKPoint([1 / 4, 1 / 4, 1 / 4], 1 / 2),
        ],
    )
    @test_throws AssertionError KPointsCard(
        "gamma",
        MonkhorstPackGrid([4, 4, 4], [1, 1, 1]),
    )
    @test_throws AssertionError KPointsCard(
        "gamma",
        [
            SpecialKPoint([3 / 4, 1 / 4, 1 / 4], 1 / 2),
            SpecialKPoint([1 / 4, 1 / 4, 1 / 4], 1 / 2),
        ],
    )
    @test_throws AssertionError KPointsCard(
        "tpiba",
        MonkhorstPackGrid([4, 4, 4], [1, 1, 1]),
    )
    @test_throws AssertionError KPointsCard("tpiba", GammaPoint())
end # testset
