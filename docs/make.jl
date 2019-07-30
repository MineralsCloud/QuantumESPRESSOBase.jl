using Documenter, QuantumESPRESSOBase

makedocs(;
    modules=[QuantumESPRESSOBase],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/MineralsCloud/QuantumESPRESSOBase.jl/blob/{commit}{path}#L{line}",
    sitename="QuantumESPRESSOBase.jl",
    authors="Qi Zhang <singularitti@outlook.com>",
    assets=String[],
)

deploydocs(;
    repo="github.com/MineralsCloud/QuantumESPRESSOBase.jl",
)
