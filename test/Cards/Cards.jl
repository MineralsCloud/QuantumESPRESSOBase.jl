using Test

using StructArrays

using QuantumESPRESSOBase.Cards.PWscf

@testset "Test constructing `AtomicSpeciesCard` from `StructArray`s" begin
    atoms = ["Al", "As"]
    masses = [24590.7655930491, 68285.4024548272]
    pseudopotentials = ["Al.pbe-n-kjpaw_psl.1.0.0.UPF", "As.pbe-n-kjpaw_psl.1.0.0.UPF"]
    @teest AtomicSpeciesCard(StructArray{AtomicSpecies}((
        atoms,
        masses,
        pseudopotentials,
    ))).data == AtomicSpecies[
        AtomicSpecies("Al", 24590.7655930491, "Al.pbe-n-kjpaw_psl.1.0.0.UPF"),
        AtomicSpecies("As", 68285.4024548272, "As.pbe-n-kjpaw_psl.1.0.0.UPF"),
    ]
end # testset
