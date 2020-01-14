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
        "Manual" => Any[
            "Installation" => "install.md",
        ],
        "Modules" => Any[
            "`QuantumESPRESSOBase` module" => "api/api.md",
            "`Namelists` module" => "api/Namelists.md",
            "`Cards` module" => "api/Cards.md",
            "`Inputs` module" => "api/Inputs.md",
        ],
    ],
)

deploydocs(;
    repo="github.com/MineralsCloud/QuantumESPRESSOBase.jl",
)
