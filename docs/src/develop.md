# How to develop this package by yourself

## Download the project

Similar to what we have mentioned in section "[Installation](@ref)", instead of running

```julia
julia> Pkg.add(PackageSpec(url="https://github.com/MineralsCloud/QuantumESPRESSOBase.jl.git"))
```

run

```julia
julia> Pkg.dev(PackageSpec(url="https://github.com/MineralsCloud/QuantumESPRESSOBase.jl.git"))
```

Then the package will be cloned to your local machine at a path. On macOS, by default is
located at `~/.julia/dev/QuantumESPRESSOBase` unless you modify the `JULIA_DEPOT_PATH`
environment variable. (See [Julia's official
documentation](http://docs.julialang.org/en/v1/manual/environment-variables/#JULIA_DEPOT_PATH-1)
on how to do this.) In the following text, we will call it `PKGROOT`.

## [Instantiate the project](@id instantiating)

Go to `PKGROOT`, start a new Julia session and run

```julia
julia> using Pkg; Pkg.instantiate()
```

## How to build docs

After [instantiating](@ref) the project, go to `PKGROOT`, run (without the `$` prompt)

```bash
$ julia --color=yes docs/make.jl
```

in your terminal. In a while a folder `PKGROOT/docs/build` will appear. Open
`PKGROOT/docs/build/index.html` with your favorite browser and have fun!

## How to run tests

After [instantiating](@ref) the project, go to `PKGROOT`, run (without the `$` prompt)

```bash
$ julia --color=yes test/runtests.jl
```

in your terminal.
