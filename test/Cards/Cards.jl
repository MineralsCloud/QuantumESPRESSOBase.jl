using Test

using Setfield
using StructArrays: StructArray

using QuantumESPRESSOBase
using QuantumESPRESSOBase.Cards.CP
using QuantumESPRESSOBase.Cards.PWscf

@testset "Constructing `AtomicSpecies`" begin
    # Data from https://github.com/QEF/q-e/blob/7be27df/PW/examples/gatefield/run_example#L128.
    x = PWscf.AtomicSpecies("S", 32.066, "S.pz-n-rrkjus_psl.0.1.UPF")
    @test_throws AssertionError @set x.atom = "sulfur"
    @test_throws InexactError @set x.mass = 1im
    @test x == PWscf.AtomicSpecies('S', 32.066, "S.pz-n-rrkjus_psl.0.1.UPF")
    y = PWscf.AtomicSpecies('S')  # Incomplete initialization
    @test_throws UndefRefError y == PWscf.AtomicSpecies("S")
    @test_throws UndefRefError y.pseudopot
    @test_throws AssertionError y.atom = "sulfur"
    x.atom = 'S'  # Setting `atom` with a `Char` still works
    @test x.atom == "S"
    @test_throws ErrorException y.mss = 12.0  # Given a wrong field name
    y.mass, y.pseudopot = 32.066, "S.pz-n-rrkjus_psl.0.1.UPF"
    @test x == y  # Constructing `AtomicSpecies` in 3 steps is equivalent to a one-time construction
    @test PWscf.AtomicSpecies(PWscf.AtomicPosition(
        'S',
        [0.500000000, 0.288675130, 1.974192764],
    )).atom == "S"
    @test PWscf.AtomicSpecies(
        PWscf.AtomicPosition('S', [0.500000000, 0.288675130, 1.974192764]),
        32.066,
        "S.pz-n-rrkjus_psl.0.1.UPF",
    ) == x
end # testset

@testset "Test constructing `AtomicSpeciesCard` from `StructArray`s" begin
    atoms = ["Al", "As"]
    masses = [24590.7655930491, 68285.4024548272]
    pseudopotentials = ["Al.pbe-n-kjpaw_psl.1.0.0.UPF", "As.pbe-n-kjpaw_psl.1.0.0.UPF"]
    init = PWscf.AtomicSpeciesCard(StructArray{PWscf.AtomicSpecies}((atoms, masses, pseudopotentials)))
    @test init.data == [
        PWscf.AtomicSpecies("Al", 24590.7655930491, "Al.pbe-n-kjpaw_psl.1.0.0.UPF"),
        PWscf.AtomicSpecies("As", 68285.4024548272, "As.pbe-n-kjpaw_psl.1.0.0.UPF"),
    ]
    push!(init.data, PWscf.AtomicSpecies("Si", 25591.1924913552, "Si.pbe-n-kjpaw_psl.1.0.0.UPF"))
    @test init.data == [
        PWscf.AtomicSpecies("Al", 24590.7655930491, "Al.pbe-n-kjpaw_psl.1.0.0.UPF"),
        PWscf.AtomicSpecies("As", 68285.4024548272, "As.pbe-n-kjpaw_psl.1.0.0.UPF"),
        PWscf.AtomicSpecies("Si", 25591.1924913552, "Si.pbe-n-kjpaw_psl.1.0.0.UPF"),
    ]
    @testset "Mutual construction" begin
        @test map(x -> x.atom, PWscf.AtomicPositionsCard("alat", init).data) == ["Al", "As", "Si"]
    end # testset
    @testset "Test `pseudopot_format`" begin
        @test unique(PWscf.pseudopot_format.(init.data)) == [PWscf.UnifiedPseudopotentialFormat()]
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
    x = PWscf.AtomicPosition("S", [0.500000000, 0.288675130, 1.974192764])
    @test_throws AssertionError @set x.atom = "sulfur"
    @test_throws InexactError @set x.pos = [1im, 2im, 3im]
    @test_throws AssertionError x.pos = [1, 2, 3, 4]
    @test_throws AssertionError x.if_pos = [1, 2, 3]
    @test_throws AssertionError x.if_pos = [1, 0, 1, 1]
    x.atom = 'S'  # Setting `atom` with a `Char` still works
    @test x.atom == "S"
    @test_throws ErrorException x.posi = [0.1, 0.2, 0.3]  # Given a wrong field name
    @test x.if_pos == [1, 1, 1]
    @test x == PWscf.AtomicPosition('S', [0.500000000, 0.288675130, 1.974192764])
    y = PWscf.AtomicPosition('S')  # Incomplete initialization
    @test_throws UndefRefError y == PWscf.AtomicPosition("S")
    @test_throws UndefRefError y.pos
    @test_throws UndefRefError y.if_pos
    @test_throws AssertionError y.atom = "sulfur"
    y.pos, y.if_pos = [0.500000000, 0.288675130, 1.974192764], [1, 1, 1]
    @test x == y  # Constructing `AtomicSpecies` in 3 steps is equivalent to a one-time construction
    @test PWscf.AtomicPosition(PWscf.AtomicSpecies('S', 32.066, "S.pz-n-rrkjus_psl.0.1.UPF")).atom ==
          "S"
    @test PWscf.AtomicPosition(
        PWscf.AtomicSpecies('S', 32.066, "S.pz-n-rrkjus_psl.0.1.UPF"),
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
    init = PWscf.AtomicPositionsCard(
        "alat",
        StructArray{PWscf.AtomicPosition}((atoms, positions, [[1, 1, 1], [1, 1, 1], [1, 1, 1]])),
    )
    @test init.data == [
        PWscf.AtomicPosition("S", [0.500000000, 0.288675130, 1.974192764]),
        PWscf.AtomicPosition("Mo", [0.000000000, 0.577350270, 2.462038339]),
        PWscf.AtomicPosition("S", [0.000000000, -0.577350270, 2.950837559]),
    ]
    @testset "Mutual construction" begin
        @test map(x -> x.atom, PWscf.AtomicSpeciesCard(init).data) == ["S", "Mo", "S"]
    end # testset
    @testset "Test `to_qe`" begin
        @test to_qe(init) ==
              "ATOMIC_POSITIONS { alat }\n      S    0.500000000    0.288675130    1.974192764\n     Mo    0.000000000    0.577350270    2.462038339\n      S    0.000000000   -0.577350270    2.950837559"
        @test to_qe(init; delim = "", indent = "") ==
              "ATOMIC_POSITIONS { alat }\n  S   0.500000000   0.288675130   1.974192764\n Mo   0.000000000   0.577350270   2.462038339\n  S   0.000000000  -0.577350270   2.950837559"
    end # testset
end # testset

@testset "Test `push_atom!`" begin
    x = PWscf.AtomicSpecies("S", 32.066, "S.pz-n-rrkjus_psl.0.1.UPF")
    v = [x]
    @test PWscf.push_atom!(v, "H", "O")[1].atom == "S"
    @test PWscf.push_atom!(v, "H", "O")[2].atom == "H"
    @test PWscf.push_atom!(v, "H", "O")[3].atom == "O"
    atoms = ["Al", "As"]
    masses = [24590.7655930491, 68285.4024548272]
    pseudopotentials = ["Al.pbe-n-kjpaw_psl.1.0.0.UPF", "As.pbe-n-kjpaw_psl.1.0.0.UPF"]
    init = PWscf.AtomicSpeciesCard(StructArray{PWscf.AtomicSpecies}((atoms, masses, pseudopotentials)))
    @test PWscf.push_atom!(init, "S", "N").data[1].atom == "Al"
    @test PWscf.push_atom!(init, "S", "N").data[2].atom == "As"
    @test PWscf.push_atom!(init, "S", "N").data[3].atom == "S"
    @test PWscf.push_atom!(init, "S", "N").data[4].atom == "N"
    x = PWscf.AtomicPosition("S", [0.500000000, 0.288675130, 1.974192764])
    v = [x]
    @test PWscf.push_atom!(v, "H", "O")[1].atom == "S"
    @test PWscf.push_atom!(v, "H", "O")[2].atom == "H"
    @test PWscf.push_atom!(v, "H", "O")[3].atom == "O"
    atoms = ["S", "Mo", "S"]
    positions = [
        [0.500000000, 0.288675130, 1.974192764],
        [0.000000000, 0.577350270, 2.462038339],
        [0.000000000, -0.577350270, 2.950837559],
    ]
    init = PWscf.AtomicPositionsCard(
        "alat",
        StructArray{PWscf.AtomicPosition}((atoms, positions, [[1, 1, 1], [1, 1, 1], [1, 1, 1]])),
    )
    @test PWscf.push_atom!(init, "H", "O").data[1].atom == "S"
    @test PWscf.push_atom!(init, "H", "O").data[2].atom == "Mo"
    @test PWscf.push_atom!(init, "H", "O").data[3].atom == "S"
    @test PWscf.push_atom!(init, "H", "O").data[4].atom == "H"
    @test PWscf.push_atom!(init, "H", "O").data[5].atom == "O"
end # testset

@testset "Test `append_atom!`" begin
    x = PWscf.AtomicSpecies("S", 32.066, "S.pz-n-rrkjus_psl.0.1.UPF")
    v = [x]
    @test PWscf.append_atom!(v, ["H", "O"])[1].atom == "S"
    @test PWscf.append_atom!(v, ["H", "O"])[2].atom == "H"
    @test PWscf.append_atom!(v, ["H", "O"])[3].atom == "O"
    atoms = ["Al", "As"]
    masses = [24590.7655930491, 68285.4024548272]
    pseudopotentials = ["Al.pbe-n-kjpaw_psl.1.0.0.UPF", "As.pbe-n-kjpaw_psl.1.0.0.UPF"]
    init = PWscf.AtomicSpeciesCard(StructArray{PWscf.AtomicSpecies}((atoms, masses, pseudopotentials)))
    @test PWscf.append_atom!(init, ["S", "N"]).data[1].atom == "Al"
    @test PWscf.append_atom!(init, ["S", "N"]).data[2].atom == "As"
    @test PWscf.append_atom!(init, ["S", "N"]).data[3].atom == "S"
    @test PWscf.append_atom!(init, ["S", "N"]).data[4].atom == "N"
    x = PWscf.AtomicPosition("S", [0.500000000, 0.288675130, 1.974192764])
    v = [x]
    @test PWscf.append_atom!(v, ["H", "O"])[1].atom == "S"
    @test PWscf.append_atom!(v, ["H", "O"])[2].atom == "H"
    @test PWscf.append_atom!(v, ["H", "O"])[3].atom == "O"
    atoms = ["S", "Mo", "S"]
    positions = [
        [0.500000000, 0.288675130, 1.974192764],
        [0.000000000, 0.577350270, 2.462038339],
        [0.000000000, -0.577350270, 2.950837559],
    ]
    init = PWscf.AtomicPositionsCard(
        "alat",
        StructArray{PWscf.AtomicPosition}((atoms, positions, [[1, 1, 1], [1, 1, 1], [1, 1, 1]])),
    )
    @test PWscf.append_atom!(init, ["H", "O"]).data[1].atom == "S"
    @test PWscf.append_atom!(init, ["H", "O"]).data[2].atom == "Mo"
    @test PWscf.append_atom!(init, ["H", "O"]).data[3].atom == "S"
    @test PWscf.append_atom!(init, ["H", "O"]).data[4].atom == "H"
    @test PWscf.append_atom!(init, ["H", "O"]).data[5].atom == "O"
end # testset

@testset "Constructing `CellParametersCard`" begin
    #Data from https://gitlab.com/QEF/q-e/blob/master/NEB/examples/neb1.in
    option = "bohr"
    data = [
        12 0  0
        0  5  0
        0  0  5
    ]
    init = PWscf.CellParametersCard(option, data)
    @test_throws AssertionError @set init.option = "ala"
    @test_throws AssertionError @set init.option = "crystal" # Allowed options are alat, angstrom, bohr
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
    @test PWscf.CellParametersCard(data).option == "alat" # default option is alat
end

@testset "Constructing `AtomicForce`" begin
    x = PWscf.AtomicForce("H", [1, 2, 3])
    @test_throws AssertionError @set x.force = [1, 2]
    @test_throws AssertionError @set x.force = [1, 2, 3, 4]
end

@testset "Test constructing `AtomicForce` from `StructArray`" begin
    atoms = ["H", "O", "H"]
    forces = [
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9],
    ]
    init = PWscf.AtomicForcesCard(StructArray{PWscf.AtomicForce}((atoms, forces)))
    @test init.data == [
        PWscf.AtomicForce("H", [1, 2, 3]),
        PWscf.AtomicForce("O", [4, 5, 6]),
        PWscf.AtomicForce("H", [7, 8, 9]),
    ]
end

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

@testset "Constructing `AtomicVelocity`" begin
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
    atomicspecies = AtomicSpecies("H", 1.00794000, "H.pbe-rrkjus_psl.1.0.0.UPF")
    velocity = [
         0.140374E-04, 
         -0.333683E-04, 
         0.231834E-04
        ]
    @test AtomicVelocity(atomicspecies).atom == "H"
    @test AtomicVelocity(atomicspecies, velocity) ==  AtomicVelocity("H", velocity)
    atomicpositions = AtomicPosition('H',[0.500000000, 0.288675130, 1.974192764])
    @test AtomicVelocity(atomicpositions).atom == "H"
    @test AtomicVelocity(atomicpositions, velocity) == AtomicVelocity("H", velocity)
end # testset

@testset "Test constructing `AtomicVelocitiesCard` from `StructArray`" begin
    atoms = ["H", "O"]
    velocities = [
        [0.140374E-04, -0.333683E-04, 0.231834E-04],
        [-0.111125E-05, -0.466724E-04,  0.972745E-05],
    ]
    init = AtomicVelocitiesCard(
        StructArray{AtomicVelocity}((atoms, velocities))
    )
    @test init.data == [
        AtomicVelocity("H", [0.140374E-04, -0.333683E-04, 0.231834E-04]),
        AtomicVelocity("O", [-0.111125E-05, -0.466724E-04, 0.972745E-05]),
    ]
end # testset

@testset "Test `push_atom!`" begin
    x = AtomicVelocity("H", [0.140374E-04, -0.333683E-04, 0.231834E-04])
    v = [x]
    @test push_atom!(v, "S", "N")[1].atom == "H"
    @test push_atom!(v, "S", "N")[2].atom == "S"
    @test push_atom!(v, "S", "N")[3].atom == "N"
    atoms = ["H", "O"]
    velocities = [
        [0.140374E-04, -0.333683E-04, 0.231834E-04],
        [-0.111125E-05, -0.466724E-04, 0.972745E-05],
    ]
    init = AtomicVelocitiesCard(
        StructArray{AtomicVelocity}((atoms, velocities))
    )
    @test push_atom!(init, "S", "N").data[1].atom == "H"
    @test push_atom!(init, "S", "N").data[2].atom == "O"
    @test push_atom!(init, "S", "N").data[3].atom == "S"
    @test push_atom!(init, "S", "N").data[4].atom == "N"
end # testset

@testset "Test `append_atom!`" begin
    x = AtomicVelocity("H", [0.140374E-04, -0.333683E-04, 0.231834E-04])
    v = [x]
    @test append_atom!(v, ["S", "N"])[1].atom == "H"
    @test append_atom!(v, ["S", "N"])[2].atom == "S"
    @test append_atom!(v, ["S", "N"])[3].atom == "N"
    atoms = ["H", "O"]
    velocities = [
        [0.140374E-04, -0.333683E-04, 0.231834E-04],
        [-0.111125E-05, -0.466724E-04, 0.972745E-05],
    ]
    init = AtomicVelocitiesCard(
        StructArray{AtomicVelocity}((atoms, velocities))
    )
    @test append_atom!(init, ["S", "N"]).data[1].atom == "H"
    @test append_atom!(init, ["S", "N"]).data[2].atom == "O"
    @test append_atom!(init, ["S", "N"]).data[3].atom == "S"
    @test append_atom!(init, ["S", "N"]).data[4].atom == "N"
end # testset

@testset "Constructing `RefCellParametersCard`"  begin
    option = "bohr"
    data = [
        12 0  0
        0  5  0
        0  0  5
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
