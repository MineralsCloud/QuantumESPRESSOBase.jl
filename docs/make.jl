using QuantumESPRESSOBase
using Documenter

DocMeta.setdocmeta!(QuantumESPRESSOBase, :DocTestSetup, :(using QuantumESPRESSOBase); recursive=true)

makedocs(;
    modules=[QuantumESPRESSOBase],
    authors="Qi Zhang <singularitti@outlook.com>",
    repo="https://github.com/MineralsCloud/QuantumESPRESSOBase.jl/blob/{commit}{path}#{line}",
    sitename="QuantumESPRESSOBase.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://MineralsCloud.github.io/QuantumESPRESSOBase.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Manual" => ["Installation" => "install.md", "Development" => "develop.md"],
        "API by module" => [
            "`QuantumESPRESSOBase` module" => "api/api.md",
            "`Inputs` module" => "api/Inputs/Inputs.md",
            "`Inputs.PWscf` module" => "api/Inputs/PWscf.md",
            "`Inputs.PHonon` module" => "api/Inputs/PHonon.md",
        ],
    ],
)

deploydocs(;
    repo="github.com/MineralsCloud/QuantumESPRESSOBase.jl",
)
