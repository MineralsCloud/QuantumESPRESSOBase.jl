module PHonon

using Test

using QuantumESPRESSOBase.Inputs.PHonon:
    PhNamelist, Q2rNamelist, MatdynNamelist, QPointsCard, PhInput, Q2rInput, MatdynInput

@testset "Construct a `PhInput`: silicon" begin
    ph = PhNamelist(;
        verbosity = "high",
        fildyn = "dyn",
        outdir = "./tmp",
        prefix = "silicon",
        ldisp = true,
        tr2_ph = 1e-14,
        nq1 = 2,
        nq2 = 2,
        nq3 = 2,
        amass = [28.086],
    )
    q_points = QPointsCard(
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
    input = PhInput(ph, q_points)
    # Test struct equality
    @test input == PhInput(deepcopy(ph), deepcopy(q_points))
end

@testset "Construct a `Q2rInput`" begin
    q2r = Q2rNamelist(; fildyn = "dyn", zasr = "crystal", flfrc = "fc.out")
    input = Q2rInput(q2r)
    # Test struct equality
    @test input == Q2rInput(deepcopy(q2r))
end

@testset "Construct a `MatdynInput`: silicon" begin
    matdyn = MatdynNamelist(;
        asr = "crystal",
        amass = [28.086],
        flfrc = "fc.out",
        flfrq = "freq.out",
        flvec = "modes.out",
        dos = true,
        q_in_band_form = false,
        nk1 = 8,
        nk2 = 8,
        nk3 = 8,
    )
    q_points = QPointsCard(
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
    input = MatdynInput(matdyn, q_points)
    # Test struct equality
    @test input == MatdynInput(deepcopy(matdyn), deepcopy(q_points))
end

end
