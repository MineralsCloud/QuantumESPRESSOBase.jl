using Test

using QuantumESPRESSOBase.Namelists.PWscf
using QuantumESPRESSOBase.Cards.PWscf
using QuantumESPRESSOBase.Inputs.PWscf

@testset "Test constructing a `PWscfInput`: silicon" begin
# This example is from https://github.com/QEF/q-e/blob/master/PW/examples/example01/run_example.
    control = ControlNamelist(
        tstress = true,
        tprnfor = true,
        outdir = raw"$TMP_DIR/",
        prefix = "silicon",
        pseudo_dir = raw"$PSEUDO_DIR/"
    )
    system = SystemNamelist(
        ibrav = 2,
        celldm = [10.2],
        nat = 2,
        ntyp = 1,
        ecutwfc = 18.0
    )
    electrons = ElectronsNamelist(
        conv_thr = 1.0e-8,
        diagonalization = raw"$diago"
    )
    atomic_species = AtomicSpeciesCard([AtomicSpecies("Si", 28.086, "Si.pz-vbc.UPF")])
    atomic_positions = AtomicPositionsCard(
        data = [
            AtomicPosition("Si", [0.0, 0.0, 0.0]),
            AtomicPosition("Si", [0.25, 0.25, 0.25])
        ]
    )
    k_points = KPointsCard(
        "tpiba",
        [
            SpecialKPoint(0.125, 0.125, 0.125, 1.0),
            SpecialKPoint(0.125, 0.125, 0.375, 3.0),
            SpecialKPoint(0.125, 0.125, 0.625, 3.0),
            SpecialKPoint(0.125, 0.125, 0.875, 3.0),
            SpecialKPoint(0.125, 0.375, 0.375, 3.0),
            SpecialKPoint(0.125, 0.375, 0.625, 6.0),
            SpecialKPoint(0.125, 0.375, 0.875, 6.0),
            SpecialKPoint(0.125, 0.625, 0.625, 3.0),
            SpecialKPoint(0.375, 0.375, 0.375, 1.0),
            SpecialKPoint(0.375, 0.375, 0.625, 3.0)
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
end # testset
