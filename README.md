<div align="center">
  <img src="https://raw.githubusercontent.com/MineralsCloud/QuantumESPRESSOBase.jl/master/docs/src/assets/logo.png" height="200"><br>
</div>

# QuantumESPRESSOBase

| **Documentation** | [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://MineralsCloud.github.io/QuantumESPRESSOBase.jl/stable/) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://MineralsCloud.github.io/QuantumESPRESSOBase.jl/dev/)                                                                                                                                                                                                                                                                                                                                     |
| :---------------: | :--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Build Status**  | [![Build Status](https://github.com/MineralsCloud/QuantumESPRESSOBase.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/MineralsCloud/QuantumESPRESSOBase.jl/actions/workflows/CI.yml?query=branch%3Amain) [![Build Status](https://ci.appveyor.com/api/projects/status/github/MineralsCloud/QuantumESPRESSOBase.jl?svg=true)](https://ci.appveyor.com/project/MineralsCloud/QuantumESPRESSOBase-jl)[![Build Status](https://api.cirrus-ci.com/github/MineralsCloud/QuantumESPRESSOBase.jl.svg)](https://cirrus-ci.com/github/MineralsCloud/QuantumESPRESSOBase.jl) |
|   **Coverage**    | [![Coverage](https://github.com/MineralsCloud/QuantumESPRESSOBase.jl/badges/main/coverage.svg)](https://github.com/MineralsCloud/QuantumESPRESSOBase.jl/commits/main) [![Coverage](https://codecov.io/gh/MineralsCloud/QuantumESPRESSOBase.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/MineralsCloud/QuantumESPRESSOBase.jl)                                                                                                                                                                                                                                                  |
|    **Others**     | [![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle) [![License](https://img.shields.io/github/license/MineralsCloud/QuantumESPRESSOBase.jl)](https://github.com/MineralsCloud/QuantumESPRESSOBase.jl/blob/main/LICENSE) [![QuantumESPRESSOBase Downloads][https://shields.io/endpoint?url=https://pkgs.genieframework.com/api/v1/badge/QuantumESPRESSOBase]][https://pkgs.genieframework.com?packages=QuantumESPRESSOBase]                                                                                            |

The code, which is [hosted on GitHub](https://github.com/MineralsCloud/QuantumESPRESSOBase.jl), is tested
using various continuous integration services for its validity.

This repository is created and maintained by
[@singularitti](https://github.com/singularitti), and contributions are highly welcome.

## Package features

[`QuantumESPRESSOBase.jl`](https://github.com/MineralsCloud/QuantumESPRESSOBase.jl) declares
basic data types and methods for manipulating crystal structures, generating input files for
[Quantum ESPRESSO](https://www.quantum-espresso.org/), error checking before running, etc.
It is written purely in the language [Julia](https://julialang.org/).

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

The code is [hosted on GitHub](https://github.com/MineralsCloud/QuantumESPRESSOBase.jl),
with some continuous integration services to test its validity.

This repository is created and maintained by [@singularitti](https://github.com/singularitti).
You are very welcome to contribute.

## Installation

The package can be installed with the Julia package manager.
From [the Julia REPL](https://docs.julialang.org/en/v1/stdlib/REPL/), type `]` to enter
the [Pkg mode](https://docs.julialang.org/en/v1/stdlib/REPL/#Pkg-mode) and run:

```julia-repl
pkg> add QuantumESPRESSOBase
```

Or, equivalently, via [`Pkg.jl`](https://pkgdocs.julialang.org/v1/):

```julia
julia> import Pkg; Pkg.add("QuantumESPRESSOBase")
```

## Documentation

- [**STABLE**][docs-stable-url] — **documentation of the most recently tagged version.**
- [**DEV**][docs-dev-url] — _documentation of the in-development version._

## Project status

The package is developed for and tested against Julia `v1.6` and above on Linux, macOS, and
Windows.

## Questions and contributions

You can post usage questions on
[our discussion page](https://github.com/MineralsCloud/QuantumESPRESSOBase.jl/discussions).

We welcome contributions, feature requests, and suggestions. If you encounter any problems,
please open an [issue](https://github.com/MineralsCloud/QuantumESPRESSOBase.jl/issues).
The [Contributing](@ref) page has
a few guidelines that should be followed when opening pull requests and contributing code.

## Star History

[![Star History Chart](https://api.star-history.com/svg?repos=MineralsCloud/QuantumESPRESSOBase.jl&type=Date)](https://star-history.com/#MineralsCloud/QuantumESPRESSOBase.jl&Date)
