<div align="center">
  <img src="https://raw.githubusercontent.com/MineralsCloud/QuantumESPRESSOBase.jl/master/docs/src/assets/logo.png" height="200"><br>
</div>

# QuantumESPRESSOBase

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://MineralsCloud.github.io/QuantumESPRESSOBase.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://MineralsCloud.github.io/QuantumESPRESSOBase.jl/dev)
[![Build Status](https://github.com/MineralsCloud/QuantumESPRESSOBase.jl/workflows/CI/badge.svg)](https://github.com/MineralsCloud/QuantumESPRESSOBase.jl/actions)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/MineralsCloud/QuantumESPRESSOBase.jl?svg=true)](https://ci.appveyor.com/project/singularitti/QuantumESPRESSOBase-jl)
[![Build Status](https://api.cirrus-ci.com/github/MineralsCloud/QuantumESPRESSOBase.jl.svg)](https://cirrus-ci.com/github/MineralsCloud/QuantumESPRESSOBase.jl)
[![pipeline status](https://gitlab.com/singularitti/QuantumESPRESSOBase.jl/badges/master/pipeline.svg)](https://gitlab.com/singularitti/QuantumESPRESSOBase.jl/-/pipelines)
[![Coverage](https://codecov.io/gh/MineralsCloud/QuantumESPRESSOBase.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/MineralsCloud/QuantumESPRESSOBase.jl)
[![PkgEval](https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/Q/QuantumESPRESSOBase.svg)](https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/report.html)
[![GitHub license](https://img.shields.io/github/license/MineralsCloud/QuantumESPRESSOBase.jl)](https://github.com/MineralsCloud/QuantumESPRESSOBase.jl/blob/master/LICENSE)

[`QuantumESPRESSOBase.jl`](https://github.com/MineralsCloud/QuantumESPRESSOBase.jl) declares
basic data types and methods for manipulating crystal structures, generating input files for
[Quantum ESPRESSO](https://www.quantum-espresso.org/), error checking before running, etc.
It is written purely in language [`Julia`](https://julialang.org/).

Please [cite this package](https://doi.org/10.1016/j.cpc.2022.108515) as:

Q. Zhang, C. Gu, J. Zhuang et al., `express`: extensible, high-level workflows for swifter *ab initio* materials modeling, *Computer Physics Communications*, 108515, doi: https://doi.org/10.1016/j.cpc.2022.108515.

The BibTeX format is: 

```bibtex
@article{ZHANG2022108515,
  title    = {express: extensible, high-level workflows for swifter ab initio materials modeling},
  journal  = {Computer Physics Communications},
  pages    = {108515},
  year     = {2022},
  issn     = {0010-4655},
  doi      = {https://doi.org/10.1016/j.cpc.2022.108515},
  url      = {https://www.sciencedirect.com/science/article/pii/S001046552200234X},
  author   = {Qi Zhang and Chaoxuan Gu and Jingyi Zhuang and Renata M. Wentzcovitch},
  keywords = {automation, workflow, high-level, high-throughput, data lineage}
}
```

We also have an [arXiv prepint](https://arxiv.org/abs/2109.11724).
