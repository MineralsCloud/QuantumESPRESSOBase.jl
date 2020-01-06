using Test

using StructArrays: StructArray

using QuantumESPRESSOBase.Cards.PWscf

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
    @testset "Test `pseudopot_format`" begin
        @test unique(pseudopot_format.(init.data)) == [UnifiedPseudopotentialFormat()]
    end # testset
end # testset

@testset "Test constructing `AtomicSpeciesCard` from `StructArray`s" begin
    # Data from https://github.com/QEF/q-e/blob/7be27df/PW/examples/gatefield/run_example#L129-L132.
    atoms = ["S", "Mo", "S"]
    positions = [
        [0.500000000, 0.288675130, 1.974192764],
        [0.000000000, 0.577350270, 2.462038339],
        [0.000000000, -0.577350270, 2.950837559],
    ]
    @test AtomicPositionsCard("alat", StructArray{AtomicPosition}((
        atoms,
        positions,
        [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
    ))).data == [
        AtomicPosition("S", [0.500000000, 0.288675130, 1.974192764]),
        AtomicPosition("Mo", [0.000000000, 0.577350270, 2.462038339]),
        AtomicPosition("S", [0.000000000, -0.577350270, 2.950837559]),
    ]
end # testset
