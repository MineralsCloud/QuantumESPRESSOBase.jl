module PWscf

using Test: @testset, @test, @test_throws
using CrystallographyBase: ReciprocalPoint
using Setfield: @set
using StructArrays: StructArray

using QuantumESPRESSOBase
using QuantumESPRESSOBase.PWscf

@testset "Test if `==` is working" begin
    @test ControlNamelist() == ControlNamelist()
    @test SystemNamelist() == SystemNamelist()
    @test ElectronsNamelist() == ElectronsNamelist()
    @test IonsNamelist() == IonsNamelist()
    @test CellNamelist() == CellNamelist()
    @test DosNamelist() == DosNamelist()
    @test BandsNamelist() == BandsNamelist()
end

@testset "Construct `AtomicSpecies`" begin
    # Data from https://github.com/QEF/q-e/blob/7be27df/PW/examples/gatefield/run_example#L128.
    x = AtomicSpecies("S", 32.066, "S.pz-n-rrkjus_psl.0.1.UPF")
    @test_throws AssertionError @set x.atom = "sulfur"
    @test_throws InexactError @set x.mass = 1im
    @test x == AtomicSpecies('S', 32.066, "S.pz-n-rrkjus_psl.0.1.UPF")
    @test x == AtomicSpecies(
        AtomicPosition('S', [0.500000000, 0.288675130, 1.974192764]),
        32.066,
        "S.pz-n-rrkjus_psl.0.1.UPF",
    )
end

@testset "Construct `AtomicSpeciesCard` from a `StructArray`" begin
    a = ["Al", "As"]
    m = [24590.7655930491, 68285.4024548272]
    pp = ["Al.pbe-n-kjpaw_psl.1.0.0.UPF", "As.pbe-n-kjpaw_psl.1.0.0.UPF"]
    card = AtomicSpeciesCard(StructArray{AtomicSpecies}((a, m, pp)))
    @test card == AtomicSpeciesCard([
        AtomicSpecies("Al", 24590.7655930491, "Al.pbe-n-kjpaw_psl.1.0.0.UPF"),
        AtomicSpecies("As", 68285.4024548272, "As.pbe-n-kjpaw_psl.1.0.0.UPF"),
    ])
    push!(card.data, AtomicSpecies("Si", 25591.1924913552, "Si.pbe-n-kjpaw_psl.1.0.0.UPF"))
    @test card == AtomicSpeciesCard([
        AtomicSpecies("Al", 24590.7655930491, "Al.pbe-n-kjpaw_psl.1.0.0.UPF"),
        AtomicSpecies("As", 68285.4024548272, "As.pbe-n-kjpaw_psl.1.0.0.UPF"),
        AtomicSpecies("Si", 25591.1924913552, "Si.pbe-n-kjpaw_psl.1.0.0.UPF"),
    ])
end

@testset "Construct `AtomicPosition`" begin
    # Data from https://github.com/QEF/q-e/blob/7be27df/PW/examples/gatefield/run_example#L129-L132.
    x = AtomicPosition("S", [0.500000000, 0.288675130, 1.974192764])
    @test_throws AssertionError @set x.atom = "sulfur"
    @test_throws InexactError @set x.pos = [1im, 2im, 3im]
    @test_throws ErrorException x.posi = [0.1, 0.2, 0.3]  # Given a wrong field name
    @test x.if_pos == [1, 1, 1]
    @test x == AtomicPosition('S', [0.500000000, 0.288675130, 1.974192764])
    @test x == AtomicPosition(
        AtomicSpecies('S', 32.066, "S.pz-n-rrkjus_psl.0.1.UPF"),
        [0.500000000, 0.288675130, 1.974192764],
    )
end

@testset "Construct `AtomicPositionsCard` from a `StructArray`" begin
    # Data from https://github.com/QEF/q-e/blob/7be27df/PW/examples/gatefield/run_example#L129-L132.
    a = ["S", "Mo", "S"]
    pos = [
        [0.500000000, 0.288675130, 1.974192764],
        [0.000000000, 0.577350270, 2.462038339],
        [0.000000000, -0.577350270, 2.950837559],
    ]
    card = AtomicPositionsCard(
        StructArray{AtomicPosition}((a, pos, [[1, 1, 1], [1, 1, 1], [1, 1, 1]])), "alat"
    )
    @test card == AtomicPositionsCard(
        [
            AtomicPosition("S", [0.500000000, 0.288675130, 1.974192764]),
            AtomicPosition("Mo", [0.000000000, 0.577350270, 2.462038339]),
            AtomicPosition("S", [0.000000000, -0.577350270, 2.950837559]),
        ],
        "alat",
    )
end

@testset "Construct `CellParametersCard`" begin
    #Data from https://gitlab.com/QEF/q-e/blob/master/NEB/examples/neb1.in
    option = "bohr"
    data = [
        12 0 0
        0 5 0
        0 0 5
    ]
    card = CellParametersCard(data, option)
    @test_throws AssertionError @set card.option = "ala"
    @test_throws AssertionError @set card.option = "crystal" # Allowed options are alat, angstrom, bohr
    @test_throws DimensionMismatch @set card.data = [ # Matrix size should be (3, 3)
        1 2
        3 4
    ]
    @test_throws DimensionMismatch @set card.data = [
        1 2 3 4
        5 6 7 8
        4 3 2 1
        8 7 6 5
    ]
    @test CellParametersCard(data).option == "alat" # default option is alat
end

@testset "Construct `AtomicForce`" begin
    x = AtomicForce("H", [1, 2, 3])
    @test_throws DimensionMismatch @set x.force = [1, 2]
    @test_throws DimensionMismatch @set x.force = [1, 2, 3, 4]
end

@testset "Construct `AtomicForce` from a `StructArray`" begin
    a = ["H", "O", "H"]
    f = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    card = AtomicForcesCard(StructArray{AtomicForce}((a, f)))
    @test card == AtomicForcesCard([
        AtomicForce("H", [1, 2, 3]),
        AtomicForce("O", [4, 5, 6]),
        AtomicForce("H", [7, 8, 9]),
    ])
end

@testset "Construct a `PWInput`: silicon" begin
    # This example is from https://github.com/QEF/q-e/blob/master/PW/examples/example01/run_example.
    for diago in ("david", "cg", "ppcg")
        control = ControlNamelist(;
            tstress=true, tprnfor=true, outdir="./", prefix="silicon", pseudo_dir="pseudo/"
        )
        system = SystemNamelist(; ibrav=2, celldm=[10.2], nat=2, ntyp=1, ecutwfc=18.0)
        electrons = ElectronsNamelist(; conv_thr=1.0e-8, diagonalization=diago)
        atomic_species = AtomicSpeciesCard([AtomicSpecies("Si", 28.086, "Si.pz-vbc.UPF")])
        atomic_positions = AtomicPositionsCard([
            AtomicPosition("Si", [0.0, 0.0, 0.0]), AtomicPosition("Si", [0.25, 0.25, 0.25])
        ])
        k_points = SpecialPointsCard([
            ReciprocalPoint(0.125, 0.125, 0.125, 1),
            ReciprocalPoint(0.125, 0.125, 0.375, 3),
            ReciprocalPoint(0.125, 0.125, 0.625, 3),
            ReciprocalPoint(0.125, 0.125, 0.875, 3),
            ReciprocalPoint(0.125, 0.375, 0.375, 3),
            ReciprocalPoint(0.125, 0.375, 0.625, 6),
            ReciprocalPoint(0.125, 0.375, 0.875, 6),
            ReciprocalPoint(0.125, 0.625, 0.625, 3),
            ReciprocalPoint(0.375, 0.375, 0.375, 1),
            ReciprocalPoint(0.375, 0.375, 0.625, 3),
        ])
        input = PWInput(;
            control=control,
            system=system,
            electrons=electrons,
            atomic_species=atomic_species,
            atomic_positions=atomic_positions,
            k_points=k_points,
        )
        @test input.electrons.diagonalization == diago
        @test input == PWInput(;
            control=deepcopy(control),
            system=deepcopy(system),
            electrons=deepcopy(electrons),
            atomic_species=deepcopy(atomic_species),
            atomic_positions=deepcopy(atomic_positions),
            k_points=deepcopy(k_points),
        )
        @test listpotentials(input) == ["Si.pz-vbc.UPF"]
        if !Sys.iswindows()
            @test getpseudodir(input) == joinpath(@__DIR__, "pseudo/")
            @test getxmldir(input) == joinpath(@__DIR__, "silicon.save")
        end
    end
end

@testset "Construct a `PWInput`: silicon bands" begin
    # This example is from https://github.com/QEF/q-e/blob/master/PW/examples/example01/run_example.
    for diago in ("david", "cg", "ppcg")
        control = ControlNamelist(;
            calculation="bands", pseudo_dir="pseudo/", outdir="./", prefix="silicon"
        )
        system = SystemNamelist(;
            ibrav=2, celldm=[10.2], nat=2, ntyp=1, ecutwfc=18.0, nbnd=8
        )
        electrons = ElectronsNamelist(; diagonalization=diago)
        atomic_species = AtomicSpeciesCard([AtomicSpecies("Si", 28.086, "Si.pz-vbc.UPF")])
        atomic_positions = AtomicPositionsCard([
            AtomicPosition("Si", [0.0, 0.0, 0.0]), AtomicPosition("Si", [0.25, 0.25, 0.25])
        ])
        k_points = SpecialPointsCard([
            ReciprocalPoint(0.0, 0.0, 0.0, 1),
            ReciprocalPoint(0.0, 0.0, 0.1, 1),
            ReciprocalPoint(0.0, 0.0, 0.2, 1),
            ReciprocalPoint(0.0, 0.0, 0.3, 1),
            ReciprocalPoint(0.0, 0.0, 0.4, 1),
            ReciprocalPoint(0.0, 0.0, 0.5, 1),
            ReciprocalPoint(0.0, 0.0, 0.6, 1),
            ReciprocalPoint(0.0, 0.0, 0.7, 1),
            ReciprocalPoint(0.0, 0.0, 0.8, 1),
            ReciprocalPoint(0.0, 0.0, 0.9, 1),
            ReciprocalPoint(0.0, 0.0, 1.0, 1),
            ReciprocalPoint(0.0, 0.0, 0.0, 1),
            ReciprocalPoint(0.0, 0.1, 0.1, 1),
            ReciprocalPoint(0.0, 0.2, 0.2, 1),
            ReciprocalPoint(0.0, 0.3, 0.3, 1),
            ReciprocalPoint(0.0, 0.4, 0.4, 1),
            ReciprocalPoint(0.0, 0.5, 0.5, 1),
            ReciprocalPoint(0.0, 0.6, 0.6, 1),
            ReciprocalPoint(0.0, 0.7, 0.7, 1),
            ReciprocalPoint(0.0, 0.8, 0.8, 1),
            ReciprocalPoint(0.0, 0.9, 0.9, 1),
            ReciprocalPoint(0.0, 1.0, 1.0, 1),
            ReciprocalPoint(0.0, 0.0, 0.0, 1),
            ReciprocalPoint(0.1, 0.1, 0.1, 1),
            ReciprocalPoint(0.2, 0.2, 0.2, 1),
            ReciprocalPoint(0.3, 0.3, 0.3, 1),
            ReciprocalPoint(0.4, 0.4, 0.4, 1),
            ReciprocalPoint(0.5, 0.5, 0.5, 1),
        ])
        input = PWInput(;
            control=control,
            system=system,
            electrons=electrons,
            atomic_species=atomic_species,
            atomic_positions=atomic_positions,
            k_points=k_points,
        )
        @test input.electrons.diagonalization == diago
        # Test whether equality holds for different constructions of `PWInput`
        @test input == PWInput(;
            control=deepcopy(control),
            system=deepcopy(system),
            electrons=deepcopy(electrons),
            atomic_species=deepcopy(atomic_species),
            atomic_positions=deepcopy(atomic_positions),
            k_points=deepcopy(k_points),
        )
        @test listpotentials(input) == ["Si.pz-vbc.UPF"]
        if !Sys.iswindows()
            @test getpseudodir(input) == joinpath(@__DIR__, "pseudo/")
            @test getxmldir(input) == joinpath(@__DIR__, "silicon.save")
        end
    end
end

@testset "Construct a `PWInput`: aluminium" begin
    # This example is from https://github.com/QEF/q-e/blob/master/PW/examples/example01/run_example.
    for diago in ("david", "cg", "ppcg")
        control = ControlNamelist(;
            calculation="scf",
            restart_mode="from_scratch",
            pseudo_dir="pseudo/",
            outdir="./",
            prefix="al",
            tprnfor=true,
            tstress=true,
        )
        system = SystemNamelist(;
            ibrav=2,
            celldm=[7.50],
            nat=1,
            ntyp=1,
            ecutwfc=15.0,
            occupations="smearing",
            smearing="marzari-vanderbilt",
            degauss=0.05,
        )
        electrons = ElectronsNamelist(; diagonalization=diago, mixing_beta=0.7)
        atomic_species = AtomicSpeciesCard([AtomicSpecies("Al", 26.98, "Al.pz-vbc.UPF")])
        atomic_positions = AtomicPositionsCard([AtomicPosition("Al", [0.0, 0.0, 0.0])])
        k_points = SpecialPointsCard([
            ReciprocalPoint(0.0625, 0.0625, 0.0625, 1),
            ReciprocalPoint(0.0625, 0.0625, 0.1875, 3),
            ReciprocalPoint(0.0625, 0.0625, 0.3125, 3),
            ReciprocalPoint(0.0625, 0.0625, 0.4375, 3),
            ReciprocalPoint(0.0625, 0.0625, 0.5625, 3),
            ReciprocalPoint(0.0625, 0.0625, 0.6875, 3),
            ReciprocalPoint(0.0625, 0.0625, 0.8125, 3),
            ReciprocalPoint(0.0625, 0.0625, 0.9375, 3),
            ReciprocalPoint(0.0625, 0.1875, 0.1875, 3),
            ReciprocalPoint(0.0625, 0.1875, 0.3125, 6),
            ReciprocalPoint(0.0625, 0.1875, 0.4375, 6),
            ReciprocalPoint(0.0625, 0.1875, 0.5625, 6),
            ReciprocalPoint(0.0625, 0.1875, 0.6875, 6),
            ReciprocalPoint(0.0625, 0.1875, 0.8125, 6),
            ReciprocalPoint(0.0625, 0.1875, 0.9375, 6),
            ReciprocalPoint(0.0625, 0.3125, 0.3125, 3),
            ReciprocalPoint(0.0625, 0.3125, 0.4375, 6),
            ReciprocalPoint(0.0625, 0.3125, 0.5625, 6),
            ReciprocalPoint(0.0625, 0.3125, 0.6875, 6),
            ReciprocalPoint(0.0625, 0.3125, 0.8125, 6),
            ReciprocalPoint(0.0625, 0.3125, 0.9375, 6),
            ReciprocalPoint(0.0625, 0.4375, 0.4375, 3),
            ReciprocalPoint(0.0625, 0.4375, 0.5625, 6),
            ReciprocalPoint(0.0625, 0.4375, 0.6875, 6),
            ReciprocalPoint(0.0625, 0.4375, 0.8125, 6),
            ReciprocalPoint(0.0625, 0.4375, 0.9375, 6),
            ReciprocalPoint(0.0625, 0.5625, 0.5625, 3),
            ReciprocalPoint(0.0625, 0.5625, 0.6875, 6),
            ReciprocalPoint(0.0625, 0.5625, 0.8125, 6),
            ReciprocalPoint(0.0625, 0.6875, 0.6875, 3),
            ReciprocalPoint(0.0625, 0.6875, 0.8125, 6),
            ReciprocalPoint(0.0625, 0.8125, 0.8125, 3),
            ReciprocalPoint(0.1875, 0.1875, 0.1875, 1),
            ReciprocalPoint(0.1875, 0.1875, 0.3125, 3),
            ReciprocalPoint(0.1875, 0.1875, 0.4375, 3),
            ReciprocalPoint(0.1875, 0.1875, 0.5625, 3),
            ReciprocalPoint(0.1875, 0.1875, 0.6875, 3),
            ReciprocalPoint(0.1875, 0.1875, 0.8125, 3),
            ReciprocalPoint(0.1875, 0.3125, 0.3125, 3),
            ReciprocalPoint(0.1875, 0.3125, 0.4375, 6),
            ReciprocalPoint(0.1875, 0.3125, 0.5625, 6),
            ReciprocalPoint(0.1875, 0.3125, 0.6875, 6),
            ReciprocalPoint(0.1875, 0.3125, 0.8125, 6),
            ReciprocalPoint(0.1875, 0.4375, 0.4375, 3),
            ReciprocalPoint(0.1875, 0.4375, 0.5625, 6),
            ReciprocalPoint(0.1875, 0.4375, 0.6875, 6),
            ReciprocalPoint(0.1875, 0.4375, 0.8125, 6),
            ReciprocalPoint(0.1875, 0.5625, 0.5625, 3),
            ReciprocalPoint(0.1875, 0.5625, 0.6875, 6),
            ReciprocalPoint(0.1875, 0.6875, 0.6875, 3),
            ReciprocalPoint(0.3125, 0.3125, 0.3125, 1),
            ReciprocalPoint(0.3125, 0.3125, 0.4375, 3),
            ReciprocalPoint(0.3125, 0.3125, 0.5625, 3),
            ReciprocalPoint(0.3125, 0.3125, 0.6875, 3),
            ReciprocalPoint(0.3125, 0.4375, 0.4375, 3),
            ReciprocalPoint(0.3125, 0.4375, 0.5625, 6),
            ReciprocalPoint(0.3125, 0.4375, 0.6875, 6),
            ReciprocalPoint(0.3125, 0.5625, 0.5625, 3),
            ReciprocalPoint(0.4375, 0.4375, 0.4375, 1),
            ReciprocalPoint(0.4375, 0.4375, 0.5625, 3),
        ])
        input = PWInput(;
            control=control,
            system=system,
            electrons=electrons,
            atomic_species=atomic_species,
            atomic_positions=atomic_positions,
            k_points=k_points,
        )
        @test input.electrons.diagonalization == diago
        @test input == PWInput(;
            control=deepcopy(control),
            system=deepcopy(system),
            electrons=deepcopy(electrons),
            atomic_species=deepcopy(atomic_species),
            atomic_positions=deepcopy(atomic_positions),
            k_points=deepcopy(k_points),
        )
        @test input.k_points == SpecialPointsCard([
            ReciprocalPoint([0.0625, 0.0625, 0.0625], 1.0),
            ReciprocalPoint([0.0625, 0.0625, 0.1875], 3.0),
            ReciprocalPoint([0.0625, 0.0625, 0.3125], 3.0),
            ReciprocalPoint([0.0625, 0.0625, 0.4375], 3.0),
            ReciprocalPoint([0.0625, 0.0625, 0.5625], 3.0),
            ReciprocalPoint([0.0625, 0.0625, 0.6875], 3.0),
            ReciprocalPoint([0.0625, 0.0625, 0.8125], 3.0),
            ReciprocalPoint([0.0625, 0.0625, 0.9375], 3.0),
            ReciprocalPoint([0.0625, 0.1875, 0.1875], 3.0),
            ReciprocalPoint([0.0625, 0.1875, 0.3125], 6.0),
            ReciprocalPoint([0.0625, 0.1875, 0.4375], 6.0),
            ReciprocalPoint([0.0625, 0.1875, 0.5625], 6.0),
            ReciprocalPoint([0.0625, 0.1875, 0.6875], 6.0),
            ReciprocalPoint([0.0625, 0.1875, 0.8125], 6.0),
            ReciprocalPoint([0.0625, 0.1875, 0.9375], 6.0),
            ReciprocalPoint([0.0625, 0.3125, 0.3125], 3.0),
            ReciprocalPoint([0.0625, 0.3125, 0.4375], 6.0),
            ReciprocalPoint([0.0625, 0.3125, 0.5625], 6.0),
            ReciprocalPoint([0.0625, 0.3125, 0.6875], 6.0),
            ReciprocalPoint([0.0625, 0.3125, 0.8125], 6.0),
            ReciprocalPoint([0.0625, 0.3125, 0.9375], 6.0),
            ReciprocalPoint([0.0625, 0.4375, 0.4375], 3.0),
            ReciprocalPoint([0.0625, 0.4375, 0.5625], 6.0),
            ReciprocalPoint([0.0625, 0.4375, 0.6875], 6.0),
            ReciprocalPoint([0.0625, 0.4375, 0.8125], 6.0),
            ReciprocalPoint([0.0625, 0.4375, 0.9375], 6.0),
            ReciprocalPoint([0.0625, 0.5625, 0.5625], 3.0),
            ReciprocalPoint([0.0625, 0.5625, 0.6875], 6.0),
            ReciprocalPoint([0.0625, 0.5625, 0.8125], 6.0),
            ReciprocalPoint([0.0625, 0.6875, 0.6875], 3.0),
            ReciprocalPoint([0.0625, 0.6875, 0.8125], 6.0),
            ReciprocalPoint([0.0625, 0.8125, 0.8125], 3.0),
            ReciprocalPoint([0.1875, 0.1875, 0.1875], 1.0),
            ReciprocalPoint([0.1875, 0.1875, 0.3125], 3.0),
            ReciprocalPoint([0.1875, 0.1875, 0.4375], 3.0),
            ReciprocalPoint([0.1875, 0.1875, 0.5625], 3.0),
            ReciprocalPoint([0.1875, 0.1875, 0.6875], 3.0),
            ReciprocalPoint([0.1875, 0.1875, 0.8125], 3.0),
            ReciprocalPoint([0.1875, 0.3125, 0.3125], 3.0),
            ReciprocalPoint([0.1875, 0.3125, 0.4375], 6.0),
            ReciprocalPoint([0.1875, 0.3125, 0.5625], 6.0),
            ReciprocalPoint([0.1875, 0.3125, 0.6875], 6.0),
            ReciprocalPoint([0.1875, 0.3125, 0.8125], 6.0),
            ReciprocalPoint([0.1875, 0.4375, 0.4375], 3.0),
            ReciprocalPoint([0.1875, 0.4375, 0.5625], 6.0),
            ReciprocalPoint([0.1875, 0.4375, 0.6875], 6.0),
            ReciprocalPoint([0.1875, 0.4375, 0.8125], 6.0),
            ReciprocalPoint([0.1875, 0.5625, 0.5625], 3.0),
            ReciprocalPoint([0.1875, 0.5625, 0.6875], 6.0),
            ReciprocalPoint([0.1875, 0.6875, 0.6875], 3.0),
            ReciprocalPoint([0.3125, 0.3125, 0.3125], 1.0),
            ReciprocalPoint([0.3125, 0.3125, 0.4375], 3.0),
            ReciprocalPoint([0.3125, 0.3125, 0.5625], 3.0),
            ReciprocalPoint([0.3125, 0.3125, 0.6875], 3.0),
            ReciprocalPoint([0.3125, 0.4375, 0.4375], 3.0),
            ReciprocalPoint([0.3125, 0.4375, 0.5625], 6.0),
            ReciprocalPoint([0.3125, 0.4375, 0.6875], 6.0),
            ReciprocalPoint([0.3125, 0.5625, 0.5625], 3.0),
            ReciprocalPoint([0.4375, 0.4375, 0.4375], 1.0),
            ReciprocalPoint([0.4375, 0.4375, 0.5625], 3.0),
        ])
        @test listpotentials(input) == ["Al.pz-vbc.UPF"]
        if !Sys.iswindows()
            @test getpseudodir(input) == joinpath(@__DIR__, "pseudo/")
            @test getxmldir(input) == joinpath(@__DIR__, "al.save")
        end
    end
end

end
