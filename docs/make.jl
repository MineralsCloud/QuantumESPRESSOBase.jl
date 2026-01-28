using QuantumESPRESSOBase
using Documenter

DocMeta.setdocmeta!(
    QuantumESPRESSOBase,
    :DocTestSetup,
    :(using QuantumESPRESSOBase, QuantumESPRESSOBase.PWscf, QuantumESPRESSOBase.PHonon);
    recursive=true,
)

makedocs(;
    modules=[QuantumESPRESSOBase],
    authors="singularitti <singularitti@outlook.com> and contributors",
    sitename="QuantumESPRESSOBase.jl",
    format=Documenter.HTML(;
        canonical="https://MineralsCloud.github.io/QuantumESPRESSOBase.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Manual" => [
            "Installation Guide" => "man/installation.md",
            "Troubleshooting" => "man/troubleshooting.md",
        ],
        "Reference" => Any[
            "Public API" => "lib/public.md",
            "Internals" => map(
                s -> "lib/internals/$(s)",
                sort(readdir(joinpath(@__DIR__, "src/lib/internals"))),
            ),
        ],
        "Developer Docs" => [
            "Contributing" => "developers/contributing.md",
            "Style Guide" => "developers/style-guide.md",
            "Design Principles" => "developers/design-principles.md",
        ],
    ],
)

deploydocs(; repo="github.com/MineralsCloud/QuantumESPRESSOBase.jl", devbranch="main")
