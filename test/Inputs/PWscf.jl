using Test

using QuantumESPRESSOBase.Namelists.PWscf
using QuantumESPRESSOBase.Cards.PWscf
using QuantumESPRESSOBase.Inputs.PWscf

@testset "Test constructing a `PWscfInput`: silicon" begin
    # This example is from https://github.com/QEF/q-e/blob/master/PW/examples/example01/run_example.
    for diago in ("david", "cg", "ppcg")
        control = ControlNamelist(
            tstress = true,
            tprnfor = true,
            outdir = raw"$TMP_DIR/",
            prefix = "silicon",
            pseudo_dir = raw"$PSEUDO_DIR/"
        )
        system = SystemNamelist(ibrav = 2, celldm = [10.2], nat = 2, ntyp = 1, ecutwfc = 18.0)
        electrons = ElectronsNamelist(conv_thr = 1.0e-8, diagonalization = "$diago")
        atomic_species = AtomicSpeciesCard([AtomicSpecies("Si", 28.086, "Si.pz-vbc.UPF")])
        atomic_positions = AtomicPositionsCard(data = [
            AtomicPosition("Si", [0.0, 0.0, 0.0]),
            AtomicPosition("Si", [0.25, 0.25, 0.25])
        ])
        k_points = KPointsCard(
            "tpiba", [
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
            ]
        )
        object = PWscfInput(
            control = control,
            system = system,
            electrons = electrons,
            atomic_species = atomic_species,
            atomic_positions = atomic_positions,
            k_points = k_points,
            cell_parameters = nothing
        )
    end
end # testset

@testset "Test constructing a `PWscfInput`: silicon bands" begin
    # This example is from https://github.com/QEF/q-e/blob/master/PW/examples/example01/run_example.
    for diago in ("david", "cg", "ppcg")
        control = ControlNamelist(
            calculation = "bands",
            pseudo_dir = raw"$PSEUDO_DIR/",
            outdir = raw"$TMP_DIR/",
            prefix = "silicon"
        )
        system = SystemNamelist(
            ibrav = 2,
            celldm = [10.2],
            nat = 2,
            ntyp = 1,
            ecutwfc = 18.0,
            nbnd = 8
        )
        electrons = ElectronsNamelist(diagonalization = "$diago")
        atomic_species = AtomicSpeciesCard([AtomicSpecies("Si", 28.086, "Si.pz-vbc.UPF")])
        atomic_positions = AtomicPositionsCard(data = [
            AtomicPosition("Si", [0.0, 0.0, 0.0]),
            AtomicPosition("Si", [0.25, 0.25, 0.25]),
        ])
        k_points = KPointsCard(data = [
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
        object = PWscfInput(
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

@testset "Test constructing a `PWscfInput`: aluminium" begin
    for diago in ("david", "cg", "ppcg")
        control = ControlNamelist(
            calculation = "scf",
            restart_mode = "from_scratch",
            pseudo_dir = raw"$PSEUDO_DIR/",
            outdir = raw"$TMP_DIR/",
            prefix = "al",
            tprnfor = true,
            tstress = true
        )
        system = SystemNamelist(
            ibrav = 2,
            celldm = [7.50],
            nat = 1,
            ntyp = 1,
            ecutwfc = 15.0,
            occupations = "smearing",
            smearing = "marzari-vanderbilt",
            degauss = 0.05
        )
        electrons = ElectronsNamelist(diagonalization = "$diago", mixing_beta = 0.7)
        atomic_species = AtomicSpeciesCard([AtomicSpecies("Al", 26.98, "Al.pz-vbc.UPF")])
        atomic_positions = AtomicPositionsCard(data = [AtomicPosition(
            "Al",
            [0.0, 0.0, 0.0],
        )])
        k_points = KPointsCard(data = [
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
        object = PWscfInput(
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
