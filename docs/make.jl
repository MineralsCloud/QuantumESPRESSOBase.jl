using QuantumESPRESSOBase
using Documenter

DocMeta.setdocmeta!(QuantumESPRESSOBase, :DocTestSetup, :(using QuantumESPRESSOBase); recursive=true)

makedocs(;
    modules=[QuantumESPRESSOBase],
    authors="singularitti <singularitti@outlook.com> and contributors",
    repo="https://github.com/MineralsCloud/QuantumESPRESSOBase.jl/blob/{commit}{path}#{line}",
    sitename="QuantumESPRESSOBase.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://MineralsCloud.github.io/QuantumESPRESSOBase.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Manual" => [
            "Installation Guide" => "installation.md",
        ],
        "Public API" => [
            "`QuantumESPRESSOBase` module" => "api/api.md",
            "`Inputs` module" => "api/Inputs/Inputs.md",
            "`Inputs.PWscf` module" => "api/Inputs/PWscf.md",
            "`Inputs.PHonon` module" => "api/Inputs/PHonon.md",
        ],
        "Developer Docs" => [
            "Contributing" => "developers/contributing.md",
            "Style Guide" => "developers/style-guide.md",
            "Design Principles" => "developers/design-principles.md",
        ],
        "Troubleshooting" => "troubleshooting.md",
    ],
)

deploydocs(;
    repo="github.com/MineralsCloud/QuantumESPRESSOBase.jl",
    devbranch="main",
)
