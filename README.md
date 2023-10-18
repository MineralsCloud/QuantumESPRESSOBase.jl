<div align="center">
  <img src="https://raw.githubusercontent.com/MineralsCloud/QuantumESPRESSOBase.jl/master/docs/src/assets/logo.png" height="200"><br>
</div>

# QuantumESPRESSOBase

|                                 **Documentation**                                  |                                                                                                 **Build Status**                                                                                                 |                                                                  **Others**                                                                   |
| :--------------------------------------------------------------------------------: | :--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------: | :-------------------------------------------------------------------------------------------------------------------------------------------: |
| [![Stable][docs-stable-img]][docs-stable-url] [![Dev][docs-dev-img]][docs-dev-url] | [![Build Status][gha-img]][gha-url] [![Build Status][appveyor-img]][appveyor-url] [![Build Status][cirrus-img]][cirrus-url] [![pipeline status][gitlab-img]][gitlab-url] [![Coverage][codecov-img]][codecov-url] | [![GitHub license][license-img]][license-url] [![Code Style: Blue][style-img]][style-url] [![QuantumESPRESSOBase Downloads][downloads-img]][downloads-url] |

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://singularitti.github.io/Spglib.jl/stable
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://singularitti.github.io/Spglib.jl/dev
[gha-img]: https://github.com/singularitti/Spglib.jl/workflows/CI/badge.svg
[gha-url]: https://github.com/singularitti/Spglib.jl/actions
[appveyor-img]: https://ci.appveyor.com/api/projects/status/github/singularitti/Spglib.jl?svg=true
[appveyor-url]: https://ci.appveyor.com/project/singularitti/Spglib-jl
[cirrus-img]: https://api.cirrus-ci.com/github/singularitti/Spglib.jl.svg
[cirrus-url]: https://cirrus-ci.com/github/singularitti/Spglib.jl
[gitlab-img]: https://gitlab.com/singularitti/Spglib.jl/badges/main/pipeline.svg
[gitlab-url]: https://gitlab.com/singularitti/Spglib.jl/-/pipelines
[codecov-img]: https://codecov.io/gh/singularitti/Spglib.jl/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/singularitti/Spglib.jl
[license-img]: https://img.shields.io/github/license/singularitti/Spglib.jl
[license-url]: https://github.com/singularitti/Spglib.jl/blob/main/LICENSE
[style-img]: https://img.shields.io/badge/code%20style-blue-4495d1.svg
[style-url]: https://github.com/invenia/BlueStyle
[downloads-img]: https://shields.io/endpoint?url=https://pkgs.genieframework.com/api/v1/badge/QuantumESPRESSOBase
[downloads-url]: https://pkgs.genieframework.com?packages=QuantumESPRESSOBase

[QuantumESPRESSOBase.jl](https://github.com/MineralsCloud/QuantumESPRESSOBase.jl) declares
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
From the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```
pkg> add QuantumESPRESSOBase
```

Or, equivalently, via the [`Pkg` API](https://pkgdocs.julialang.org/v1/getting-started/):

```julia
julia> import Pkg; Pkg.add("QuantumESPRESSOBase")
```

## Documentation

- [**STABLE**][docs-stable-url] — **documentation of the most recently tagged version.**
- [**DEV**][docs-dev-url] — _documentation of the in-development version._

## Project status

The package is tested against, and being developed for, Julia `1.6` and above on Linux,
macOS, and Windows.

## Questions and contributions

You are welcome to post usage questions on [our discussion page][discussions-url].

Contributions are very welcome, as are feature requests and suggestions. Please open an
[issue][issues-url] if you encounter any problems. The [Contributing](@ref) page has
guidelines that should be followed when opening pull requests and contributing code.

[discussions-url]: https://github.com/MineralsCloud/QuantumESPRESSOBase.jl/discussions
[issues-url]: https://github.com/MineralsCloud/QuantumESPRESSOBase.jl/issues

## Star History

[![Star History Chart](https://api.star-history.com/svg?repos=MineralsCloud/QuantumESPRESSOBase.jl&type=Date)](https://star-history.com/#MineralsCloud/QuantumESPRESSOBase.jl&Date)
