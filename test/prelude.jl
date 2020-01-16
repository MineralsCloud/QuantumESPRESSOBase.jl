using Test

using QuantumESPRESSOBase

# Data are calcuted from http://calistry.org/viz/direct-and-reciprocal-lattice-visualizer
@testset "Test `reciprocal lattice` " begin
    @testset "`ibrav = 1` simple cubic" begin
        ibrav = 1
        celldm = [1.0, nothing]
        @test round.(reciprocal_lattice(ibrav, celldm), digits = 5) == [
            6.28319 0.0 0.0
            0.0 6.28319 0.0
            0.0 0.0 6.28319
        ]
    end
    @testset "`ibrav = 2` face centered cubic" begin
        ibrav = 2
        celldm = [1.0, nothing]
        @test round.(reciprocal_lattice(ibrav, celldm), digits = 5) == [
            -6.28319 6.28319 -6.28319
            -6.28319 6.28319 6.28319
            6.28319 6.28319 -6.28319
        ]
    end
    @testset "`ibrav = 3` body centered cubic" begin
        ibrav = 3
        celldm = [1.0, nothing]
        @test round.(reciprocal_lattice(ibrav, celldm), digits = 5) == [
            6.28319 -6.28319 0.0
            0.0 6.28319 -6.28319
            6.28319 0.0 6.28319
        ]
    end
    @testset "`ibrav = -3` body centered cubic, more symmetry axis" begin
        ibrav = -3
        celldm = [1.0, nothing]
        @test round.(reciprocal_lattice(ibrav, celldm), digits = 5) == [
            0.0 6.28319 6.28319
            6.28319 0.0 6.28319
            6.28319 6.28319 0.0
        ]
    end
    @testset "`ibrav = 4` Hexagonal and trigonal P" begin
        ibrav = 4
        celldm = [1.0, nothing, 2.0]
        @test round.(reciprocal_lattice(ibrav, celldm), digits = 5) == [
            6.28319 0.0 0.0
            3.6276 7.2552 -0.0
            -0.0 0.0 3.14159
        ]
    end
    @testset "`ibrav = 5 Trigonal R`, 3fold axis c" begin
        ibrav = 5
        celldm = [1.0, nothing, nothing, 0.5]
        @test round.(reciprocal_lattice(ibrav, celldm), digits = 5) == [
            6.28319 0.0 -6.28319
            -3.6276 7.2552 -3.6276
            2.5651 2.5651 2.5651
        ]
    end
    @testset "`ibrav = -5`, 3fold axis <111>" begin
        ibrav = -5
        celldm = [1.0, nothing, nothing, 0.5]
        @test round.(reciprocal_lattice(ibrav, celldm), digits = 5) == [
            -4.44288 4.44288 4.44288
            4.44288 -4.44288 4.44288
            4.44288 4.44288 -4.44288
        ]
    end
    @testset "`ibrav = 6` Tetragonal P" begin
        ibrav = 6
        celldm = [1.0, nothing, 1.0]
        @test round.(reciprocal_lattice(ibrav, celldm), digits = 5) == [
            6.28319 0.0 0.0
            0.0 6.28319 0.0
            0.0 0.0 6.28319
        ]
    end
    @testset "`ibrav = 7` Tetragonal I" begin
        ibrav = 7
        celldm = [1.0, nothing, 1.0]
        @test round.(reciprocal_lattice(ibrav, celldm), digits = 5) == [
            6.28319 0.0 -6.28319
            -6.28319 6.28319 0.0
            0.0 6.28319 6.28319
        ]
    end
    @testset "`ibrav = 8` Orthorhombic P" begin
        ibrav = 8
        celldm = [1.0, 2.0, 3.0, nothing]
        @test round.(reciprocal_lattice(ibrav, celldm), digits = 5) == [
            6.28319 0.0 0.0
            0.0 3.14159 0.0
            0.0 0.0 2.0944
        ]
    end
    @testset "`ibrav = 9` Orthorhombic base-centered" begin
        ibrav = 9
        celldm = [1.0, 2.0, 3.0, nothing]
        @test round.(reciprocal_lattice(ibrav, celldm), digits = 5) == [
            6.28319 -6.28319 0.0
            3.14159 3.14159 -0.0
            -0.0 0.0 2.0944
        ]
    end
    @testset "`ibrav = -9`" begin
        ibrav = -9
        celldm = [1.0, 2.0, 3.0, nothing]
        @test round.(reciprocal_lattice(ibrav, celldm), digits = 5) == [
            6.28319 6.28319 -0.0
            -3.14159 3.14159 0.0
            0.0 -0.0 2.0944
        ]
    end
    @testset "`ibrav = 91` Orthorhombic one-face base-centered A-type" begin
        ibrav = 91
        celldm = [1.0, 2.0, 3.0, nothing]
        @test round.(reciprocal_lattice(ibrav, celldm), digits = 5) == [
            6.28319 0.0 -0.0
            -0.0 3.14159 3.14159
            0.0 -2.0944 2.0944
        ]
    end
    @testset "`ibrav = 10` Orthorhombic face-centered" begin
        ibrav = 10
        celldm = [1.0, 2.0, 3.0, nothing]
        @test round.(reciprocal_lattice(ibrav, celldm), digits = 5) == [
            6.28319 6.28319 -6.28319
            -3.14159 3.14159 3.14159
            2.0944 -2.0944 2.0944
        ]
    end
    @testset "`ibrav = 11` Orthorhombic body-centered" begin
        ibrav = 11
        celldm = [1.0, 2.0, 3.0, nothing]
        @test round.(reciprocal_lattice(ibrav, celldm), digits = 5) == [
            6.28319 6.28319 -6.28319
            -3.14159 3.14159 3.14159
            2.0944 -2.0944 2.0944
        ]
    end
    @testset "`ibrav = 12` Monoclinic P, unique axis c" begin
        ibrav = 12
        celldm = [1.0, 2.0, 3.0, 0.5, nothing]
        @test round.(reciprocal_lattice(ibrav, celldm), digits = 5) == [
            6.28319 0.0 0.0
            -3.6276 3.6276 0.0
            0.0 0.0 2.0944
        ]
    end
    @testset "`ibrav = -12` Monoclinic P, unique axis b" begin
        ibrav = -12
        celldm = [1.0, 2.0, 3.0, nothing, 0.5]
        @test round.(reciprocal_lattice(ibrav, celldm), digits = 5) == [
            6.28319 0.0 0.0
            0.0 3.14159 0.0
            -3.6276 0.0 2.4184
        ]
    end
    @testset "`ibrav = 13` Monoclinic base-centered(unique axis b) " begin
        ibrav = 13
        celldm = [1.0, 2.0, 3.0, 0.5, nothing]
        @test round.(reciprocal_lattice(ibrav, celldm), digits = 5) == [
            6.28319 -0.0 6.28319
            -3.6276 3.6276 -3.6276
            -2.0944 0.0 2.0944
        ]
    end
    @testset "`ibrav = -13` Monoclinic base-centered(unique axis c) " begin
        ibrav = -13
        celldm = [1.0, 2.0, 3.0, nothing, 0.5]
        @test round.(reciprocal_lattice(ibrav, celldm), digits = 5) == [
            6.28319 6.28319 -0.0
            -3.14159 3.14159 0.0
            -3.6276 -3.6276 2.4184
        ]
    end
    @testset "`ibrav = 14` Triclinic" begin
        ibrav = 14
        celldm = [1.0, 2.0, 3.0, 0.2, 0.3, 0.4, nothing]
        @test round.(reciprocal_lattice(ibrav, celldm), digits = 5) == [
            6.28319 0.0 0.0
            -2.74221 3.42776 0.0
            -1.73232 -0.31497 2.20477
        ]
    end
end
