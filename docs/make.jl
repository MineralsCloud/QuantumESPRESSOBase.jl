using QuantumESPRESSOBase
using Documenter

makedocs(;
    modules=[QuantumESPRESSOBase],
    authors="Qi Zhang <singularitti@outlook.com>",
    repo="https://github.com/MineralsCloud/QuantumESPRESSOBase.jl/blob/{commit}{path}#L{line}",
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
            "`Inputs` module" => "api/Inputs.md",
        ],
    ],
)

deploydocs(;
    repo="github.com/MineralsCloud/QuantumESPRESSOBase.jl",
)
