module PWscf

using Test

using Pseudopotentials: UnifiedPseudopotentialFormat, pseudopot_format
using Setfield
using StructArrays: StructArray

using QuantumESPRESSOBase
using QuantumESPRESSOBase.Inputs.PWscf

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
    x.atom = 'S'  # Setting `atom` with a `Char` still works
    @test x.atom == "S"
    @test_throws ErrorException y.mss = 12.0  # Given a wrong field name
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
    a = ["Al", "As"]
    m = [24590.7655930491, 68285.4024548272]
    pp = ["Al.pbe-n-kjpaw_psl.1.0.0.UPF", "As.pbe-n-kjpaw_psl.1.0.0.UPF"]
    card = AtomicSpeciesCard(StructArray{AtomicSpecies}((a, m, pp)))
    @test card.data == [
        AtomicSpecies("Al", 24590.7655930491, "Al.pbe-n-kjpaw_psl.1.0.0.UPF"),
        AtomicSpecies("As", 68285.4024548272, "As.pbe-n-kjpaw_psl.1.0.0.UPF"),
    ]
    push!(card.data, AtomicSpecies("Si", 25591.1924913552, "Si.pbe-n-kjpaw_psl.1.0.0.UPF"))
    @test card.data == [
        AtomicSpecies("Al", 24590.7655930491, "Al.pbe-n-kjpaw_psl.1.0.0.UPF"),
        AtomicSpecies("As", 68285.4024548272, "As.pbe-n-kjpaw_psl.1.0.0.UPF"),
        AtomicSpecies("Si", 25591.1924913552, "Si.pbe-n-kjpaw_psl.1.0.0.UPF"),
    ]
    @testset "Mutual construction" begin
        @test map(x -> x.atom, AtomicPositionsCard("alat", card).data) == ["Al", "As", "Si"]
    end # testset
    @testset "Test `pseudopot_format`" begin
        @test unique(pseudopot_format.(card.data)) == [UnifiedPseudopotentialFormat()]
    end # testset
    @testset "Test `qestring`" begin
        @test qestring(card) ==
              "ATOMIC_SPECIES\n     Al     24590.7655930491 Al.pbe-n-kjpaw_psl.1.0.0.UPF\n     As     68285.4024548272 As.pbe-n-kjpaw_psl.1.0.0.UPF\n     Si     25591.1924913552 Si.pbe-n-kjpaw_psl.1.0.0.UPF"
        @test qestring(card; delim = "", indent = "") ==
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
    x.atom = 'S'  # Setting `atom` with a `Char` still works
    @test x.atom == "S"
    @test_throws ErrorException x.posi = [0.1, 0.2, 0.3]  # Given a wrong field name
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
    a = ["S", "Mo", "S"]
    pos = [
        [0.500000000, 0.288675130, 1.974192764],
        [0.000000000, 0.577350270, 2.462038339],
        [0.000000000, -0.577350270, 2.950837559],
    ]
    card = AtomicPositionsCard(
        "alat",
        StructArray{AtomicPosition}((a, pos, [[1, 1, 1], [1, 1, 1], [1, 1, 1]])),
    )
    @test card.data == [
        AtomicPosition("S", [0.500000000, 0.288675130, 1.974192764]),
        AtomicPosition("Mo", [0.000000000, 0.577350270, 2.462038339]),
        AtomicPosition("S", [0.000000000, -0.577350270, 2.950837559]),
    ]
    @testset "Mutual construction" begin
        @test map(x -> x.atom, AtomicSpeciesCard(card).data) == ["S", "Mo", "S"]
    end # testset
    @testset "Test `qestring`" begin
        @test qestring(card) ==
              "ATOMIC_POSITIONS { alat }\n      S    0.500000000    0.288675130    1.974192764\n     Mo    0.000000000    0.577350270    2.462038339\n      S    0.000000000   -0.577350270    2.950837559"
        @test qestring(card; delim = "", indent = "") ==
              "ATOMIC_POSITIONS { alat }\n  S   0.500000000   0.288675130   1.974192764\n Mo   0.000000000   0.577350270   2.462038339\n  S   0.000000000  -0.577350270   2.950837559"
    end # testset
end # testset

@testset "Test `push_atom!`" begin
    @testset "`push_atom!` to `AtomicSpecies`" begin
        v = [AtomicSpecies("S", 32.066, "S.pz-n-rrkjus_psl.0.1.UPF")]
        push_atom!(v, "H", "O")
        @test [x.atom for x in v] == ["S", "H", "O"]
    end # testset
    @testset "`push_atom!` to `AtomicSpeciesCard`" begin
        a = ["Al", "As"]
        m = [24590.7655930491, 68285.4024548272]
        pp = ["Al.pbe-n-kjpaw_psl.1.0.0.UPF", "As.pbe-n-kjpaw_psl.1.0.0.UPF"]
        card = AtomicSpeciesCard(StructArray{AtomicSpecies}((a, m, pp)))
        push_atom!(card, "S", "N")
        @test [x.atom for x in card.data] == ["Al", "As", "S", "N"]
    end # testset
    @testset "`push_atom!` to `AtomicPosition`" begin
        v = [AtomicPosition("S", [0.500000000, 0.288675130, 1.974192764])]
        push_atom!(v, "H", "O")
        @test [x.atom for x in v] == ["S", "H", "O"]
    end # testset
    @testset "`push_atom!` to `AtomicPositionsCard`" begin
        a = ["S", "Mo", "S"]
        pos = [
            [0.500000000, 0.288675130, 1.974192764],
            [0.000000000, 0.577350270, 2.462038339],
            [0.000000000, -0.577350270, 2.950837559],
        ]
        card = AtomicPositionsCard(
            "alat",
            StructArray{AtomicPosition}((a, pos, [[1, 1, 1], [1, 1, 1], [1, 1, 1]])),
        )
        push_atom!(card, "H", "O")
        @test [x.atom for x in card.data] == ["S", "Mo", "S", "H", "O"]
    end # testset
end # testset

@testset "Test `append_atom!`" begin
    @testset "`append_atom!` to `AtomicSpecies`" begin
        v = [AtomicSpecies("S", 32.066, "S.pz-n-rrkjus_psl.0.1.UPF")]
        append_atom!(v, ["H", "O"])
        @test [x.atom for x in v] == ["S", "H", "O"]
    end # testset
    @testset "`append_atom!` to `AtomicSpeciesCard`" begin
        a = ["Al", "As"]
        m = [24590.7655930491, 68285.4024548272]
        pp = ["Al.pbe-n-kjpaw_psl.1.0.0.UPF", "As.pbe-n-kjpaw_psl.1.0.0.UPF"]
        card = AtomicSpeciesCard(StructArray{AtomicSpecies}((a, m, pp)))
        append_atom!(card, ["S", "N"])
        @test [x.atom for x in card.data] == ["Al", "As", "S", "N"]
    end # testset
    @testset "`append_atom!` to `AtomicPosition`" begin
        v = [AtomicPosition("S", [0.500000000, 0.288675130, 1.974192764])]
        append_atom!(v, ["H", "O"])
        @test [x.atom for x in v] == ["S", "H", "O"]
    end # testset
    @testset "`append_atom!` to `AtomicPositionsCard`" begin
        a = ["S", "Mo", "S"]
        pos = [
            [0.500000000, 0.288675130, 1.974192764],
            [0.000000000, 0.577350270, 2.462038339],
            [0.000000000, -0.577350270, 2.950837559],
        ]
        card = AtomicPositionsCard(
            "alat",
            StructArray{AtomicPosition}((a, pos, [[1, 1, 1], [1, 1, 1], [1, 1, 1]])),
        )
        append_atom!(card, ["H", "O"])
        @test [x.atom for x in card.data] == ["S", "Mo", "S", "H", "O"]
    end # testset
end # testset

@testset "Constructing `CellParametersCard`" begin
    #Data from https://gitlab.com/QEF/q-e/blob/master/NEB/examples/neb1.in
    option = "bohr"
    data = [
        12 0 0
        0 5 0
        0 0 5
    ]
    card = CellParametersCard(option, data)
    @test_throws AssertionError @set card.option = "ala"
    @test_throws AssertionError @set card.option = "crystal" # Allowed options are alat, angstrom, bohr
    @test_throws AssertionError @set card.data = [ # Matrix size should be (3, 3)
        1 2
        3 4
    ]
    @test_throws AssertionError @set card.data = [
        1 2 3 4
        5 6 7 8
        4 3 2 1
        8 7 6 5
    ]
    @test CellParametersCard(data).option == "alat" # default option is alat
end

@testset "Constructing `AtomicForce`" begin
    x = AtomicForce("H", [1, 2, 3])
    @test_throws AssertionError @set x.force = [1, 2]
    @test_throws AssertionError @set x.force = [1, 2, 3, 4]
end

@testset "Test constructing `AtomicForce` from `StructArray`" begin
    a = ["H", "O", "H"]
    f = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    card = AtomicForcesCard(StructArray{AtomicForce}((a, f)))
    @test card.data == [
        AtomicForce("H", [1, 2, 3]),
        AtomicForce("O", [4, 5, 6]),
        AtomicForce("H", [7, 8, 9]),
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

end # module PWscf

@testset "Test constructing a `PWInput`: silicon" begin
    # This example is from https://github.com/QEF/q-e/blob/master/PW/examples/example01/run_example.
    for diago in ("david", "cg", "ppcg")
        control = ControlNamelist(
            tstress = true,
            tprnfor = true,
            outdir = raw"$TMP_DIR/",
            prefix = "silicon",
            pseudo_dir = raw"$PSEUDO_DIR/",
        )
        system =
            SystemNamelist(ibrav = 2, celldm = [10.2], nat = 2, ntyp = 1, ecutwfc = 18.0)
        electrons = ElectronsNamelist(conv_thr = 1.0e-8, diagonalization = "$diago")
        atomic_species = AtomicSpeciesCard([AtomicSpecies("Si", 28.086, "Si.pz-vbc.UPF")])
        atomic_positions = AtomicPositionsCard([
            AtomicPosition("Si", [0.0, 0.0, 0.0]),
            AtomicPosition("Si", [0.25, 0.25, 0.25]),
        ])
        k_points = KPointsCard(
            "tpiba",
            [
                0.125 0.125 0.125 1.0
                0.125 0.125 0.375 3.0
                0.125 0.125 0.625 3.0
                0.125 0.125 0.875 3.0
                0.125 0.375 0.375 3.0
                0.125 0.375 0.625 6.0
                0.125 0.375 0.875 6.0
                0.125 0.625 0.625 3.0
                0.375 0.375 0.375 1.0
                0.375 0.375 0.625 3.0
            ],
        )
        object = PWInput(
            control = control,
            system = system,
            electrons = electrons,
            atomic_species = atomic_species,
            atomic_positions = atomic_positions,
            k_points = k_points,
            cell_parameters = nothing,
        )
    end
end # testset

@testset "Test constructing a `PWInput`: silicon bands" begin
    # This example is from https://github.com/QEF/q-e/blob/master/PW/examples/example01/run_example.
    for diago in ("david", "cg", "ppcg")
        control = ControlNamelist(
            calculation = "bands",
            pseudo_dir = raw"$PSEUDO_DIR/",
            outdir = raw"$TMP_DIR/",
            prefix = "silicon",
        )
        system = SystemNamelist(
            ibrav = 2,
            celldm = [10.2],
            nat = 2,
            ntyp = 1,
            ecutwfc = 18.0,
            nbnd = 8,
        )
        electrons = ElectronsNamelist(diagonalization = "$diago")
        atomic_species = AtomicSpeciesCard([AtomicSpecies("Si", 28.086, "Si.pz-vbc.UPF")])
        atomic_positions = AtomicPositionsCard([
            AtomicPosition("Si", [0.0, 0.0, 0.0]),
            AtomicPosition("Si", [0.25, 0.25, 0.25]),
        ])
        k_points = KPointsCard([
            0.0 0.0 0.0 1.0
            0.0 0.0 0.1 1.0
            0.0 0.0 0.2 1.0
            0.0 0.0 0.3 1.0
            0.0 0.0 0.4 1.0
            0.0 0.0 0.5 1.0
            0.0 0.0 0.6 1.0
            0.0 0.0 0.7 1.0
            0.0 0.0 0.8 1.0
            0.0 0.0 0.9 1.0
            0.0 0.0 1.0 1.0
            0.0 0.0 0.0 1.0
            0.0 0.1 0.1 1.0
            0.0 0.2 0.2 1.0
            0.0 0.3 0.3 1.0
            0.0 0.4 0.4 1.0
            0.0 0.5 0.5 1.0
            0.0 0.6 0.6 1.0
            0.0 0.7 0.7 1.0
            0.0 0.8 0.8 1.0
            0.0 0.9 0.9 1.0
            0.0 1.0 1.0 1.0
            0.0 0.0 0.0 1.0
            0.1 0.1 0.1 1.0
            0.2 0.2 0.2 1.0
            0.3 0.3 0.3 1.0
            0.4 0.4 0.4 1.0
            0.5 0.5 0.5 1.0
        ])
        object = PWInput(
            control = control,
            system = system,
            electrons = electrons,
            atomic_species = atomic_species,
            atomic_positions = atomic_positions,
            k_points = k_points,
            cell_parameters = nothing,
        )
    end
end # testset

@testset "Test constructing a `PWInput`: aluminium" begin
    # This example is from https://github.com/QEF/q-e/blob/master/PW/examples/example01/run_example.
    for diago in ("david", "cg", "ppcg")
        control = ControlNamelist(
            calculation = "scf",
            restart_mode = "from_scratch",
            pseudo_dir = raw"$PSEUDO_DIR/",
            outdir = raw"$TMP_DIR/",
            prefix = "al",
            tprnfor = true,
            tstress = true,
        )
        system = SystemNamelist(
            ibrav = 2,
            celldm = [7.50],
            nat = 1,
            ntyp = 1,
            ecutwfc = 15.0,
            occupations = "smearing",
            smearing = "marzari-vanderbilt",
            degauss = 0.05,
        )
        electrons = ElectronsNamelist(diagonalization = "$diago", mixing_beta = 0.7)
        atomic_species = AtomicSpeciesCard([AtomicSpecies("Al", 26.98, "Al.pz-vbc.UPF")])
        atomic_positions = AtomicPositionsCard([AtomicPosition("Al", [0.0, 0.0, 0.0])])
        k_points = KPointsCard([
            0.0625000 0.0625000 0.0625000 1.00
            0.0625000 0.0625000 0.1875000 3.00
            0.0625000 0.0625000 0.3125000 3.00
            0.0625000 0.0625000 0.4375000 3.00
            0.0625000 0.0625000 0.5625000 3.00
            0.0625000 0.0625000 0.6875000 3.00
            0.0625000 0.0625000 0.8125000 3.00
            0.0625000 0.0625000 0.9375000 3.00
            0.0625000 0.1875000 0.1875000 3.00
            0.0625000 0.1875000 0.3125000 6.00
            0.0625000 0.1875000 0.4375000 6.00
            0.0625000 0.1875000 0.5625000 6.00
            0.0625000 0.1875000 0.6875000 6.00
            0.0625000 0.1875000 0.8125000 6.00
            0.0625000 0.1875000 0.9375000 6.00
            0.0625000 0.3125000 0.3125000 3.00
            0.0625000 0.3125000 0.4375000 6.00
            0.0625000 0.3125000 0.5625000 6.00
            0.0625000 0.3125000 0.6875000 6.00
            0.0625000 0.3125000 0.8125000 6.00
            0.0625000 0.3125000 0.9375000 6.00
            0.0625000 0.4375000 0.4375000 3.00
            0.0625000 0.4375000 0.5625000 6.00
            0.0625000 0.4375000 0.6875000 6.00
            0.0625000 0.4375000 0.8125000 6.00
            0.0625000 0.4375000 0.9375000 6.00
            0.0625000 0.5625000 0.5625000 3.00
            0.0625000 0.5625000 0.6875000 6.00
            0.0625000 0.5625000 0.8125000 6.00
            0.0625000 0.6875000 0.6875000 3.00
            0.0625000 0.6875000 0.8125000 6.00
            0.0625000 0.8125000 0.8125000 3.00
            0.1875000 0.1875000 0.1875000 1.00
            0.1875000 0.1875000 0.3125000 3.00
            0.1875000 0.1875000 0.4375000 3.00
            0.1875000 0.1875000 0.5625000 3.00
            0.1875000 0.1875000 0.6875000 3.00
            0.1875000 0.1875000 0.8125000 3.00
            0.1875000 0.3125000 0.3125000 3.00
            0.1875000 0.3125000 0.4375000 6.00
            0.1875000 0.3125000 0.5625000 6.00
            0.1875000 0.3125000 0.6875000 6.00
            0.1875000 0.3125000 0.8125000 6.00
            0.1875000 0.4375000 0.4375000 3.00
            0.1875000 0.4375000 0.5625000 6.00
            0.1875000 0.4375000 0.6875000 6.00
            0.1875000 0.4375000 0.8125000 6.00
            0.1875000 0.5625000 0.5625000 3.00
            0.1875000 0.5625000 0.6875000 6.00
            0.1875000 0.6875000 0.6875000 3.00
            0.3125000 0.3125000 0.3125000 1.00
            0.3125000 0.3125000 0.4375000 3.00
            0.3125000 0.3125000 0.5625000 3.00
            0.3125000 0.3125000 0.6875000 3.00
            0.3125000 0.4375000 0.4375000 3.00
            0.3125000 0.4375000 0.5625000 6.00
            0.3125000 0.4375000 0.6875000 6.00
            0.3125000 0.5625000 0.5625000 3.00
            0.4375000 0.4375000 0.4375000 1.00
            0.4375000 0.4375000 0.5625000 3.00
        ])
        object = PWInput(
            control = control,
            system = system,
            electrons = electrons,
            atomic_species = atomic_species,
            atomic_positions = atomic_positions,
            k_points = k_points,
            cell_parameters = nothing,
        )
    end
end # testset