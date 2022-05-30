# Installation

## Install Julia

First, you should install Julia. We recommend downloading it from
[its official website](https://julialang.org/downloads/). Versions higher than `v1.3`,
especially `v1.6`, are strongly recommended. This package may not work on `v0.7` and below.
Please follow the detailed instructions on its website if you have to
[build Julia from source](https://github.com/JuliaLang/julia/blob/master/doc/build/build.md).
Some computing centers provide preinstalled Julia. Please contact your administrator for
more information in that case.

If you have [Homebrew](https://brew.sh) installed, open
`Terminal.app` and type

```shell
$ brew install --cask julia  # Install a binary, works only on macOS
```

or

```shell
$ brew install julia  # Build from source, works on all operating systems
```

On Linux, the best way to install Julia is to use the *Generic Linux Binaries*. The
[JILL](https://github.com/abelsiqueira/jill) script does this for you. Run

```shell
$ bash -ci "$(curl -fsSL https://raw.githubusercontent.com/abelsiqueira/jill/master/jill.sh)"
```

to install Julia in `$HOME/.local/bin`. This script also has a Python version,
[JILL.py](https://github.com/johnnychen94/jill.py). It can also be used on macOS.

If you want to install multiple Julia versions in the same operating system,
a suggested way is to use a version manager such as
[`asdf`](https://asdf-vm.com/guide/introduction.html).
First, [install `asdf`](https://asdf-vm.com/guide/getting-started.html#_3-install-asdf).
Then, run

```shell
$ asdf install julia 1.6.2   # or other versions of Julia
$ asdf global julia 1.6.2
```

to install Julia and
[set `v1.6.2` as a global version](https://asdf-vm.com/guide/getting-started.html#_6-set-a-version).

## Install `QuantumESPRESSOBase`

Now I am using [macOS](https://en.wikipedia.org/wiki/MacOS) as a standard
platform to explain the following steps:

1. Open `Terminal.app`, and type `julia` to start an interactive session (known as
   [REPL](https://docs.julialang.org/en/v1/stdlib/REPL/)).

2. Run the following commands and wait for them to finish:

   ```julia
   julia> using Pkg

   julia> Pkg.add("QuantumESPRESSOBase")
   ```

3. Run

   ```julia
   julia> using QuantumESPRESSOBase
   ```

   and have fun!

4. While using, please keep this Julia session alive. Restarting might recompile
   the package and cost some time.

If you want to install the latest in development (maybe buggy) version of `QuantumESPRESSOBase`, type

```julia
julia> using Pkg

julia> pkg"add QuantumESPRESSOBase#master"
```

in the second step instead.

## Uninstall and reinstall `QuantumESPRESSOBase`

1. To uninstall, in a Julia session, run

   ```julia
   julia> Pkg.rm("QuantumESPRESSOBase"); Pkg.gc()
   ```

2. Press `ctrl+d` to quit the current session. Start a new Julia session and
   reinstall `QuantumESPRESSOBase`.
