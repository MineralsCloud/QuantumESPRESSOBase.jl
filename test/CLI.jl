module CLI

using Test

using QuantumESPRESSOBase.CLI: PWExec, PhExec, Q2rExec, MatdynExec, scriptify, binpath

@testset "Create `QuantumESPRESSOExec`s" begin
    @test binpath(PWExec) == "pw.x"
    @test binpath(PhExec) == "ph.x"
    @test binpath(Q2rExec) == "q2r.x"
    @test binpath(MatdynExec) == "matdyn.x"
end

@testset "Scriptify" begin
    @test_throws UndefKeywordError scriptify(PWExec(nimage = 2, npool = 3, ntg = 4))
end

end
