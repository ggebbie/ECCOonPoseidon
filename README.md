# ECCOonPoseidon

Analysis of ECCO version 4 release 4 using output from Poseidon @ WHOI. A scientific project that relies heavily on the [ECCOtour](https://github.com/ggebbie/ECCOtour.jl), MITgcmTools, and MeshArrays Julia packages.

# Getting started

Clone with `git clone https://github.com/ggebbie/ECCOonPoseidon` \
Or `git clone git@github.com:ggebbie/ECCOonPoseidon.git` with SSH keys.

# Requirements
- Julia 1.6
- matplotlib

Install a Julia-specific version of Python with Matplotlib through the Julia REPL:\
`ENV["PYTHON"] = ""`\
`# make sure python command is available, not just python 3`\
`using Pkg`\
`pkg"add PyCall"`\
`Pkg.build("PyCall")`\
`using Conda`\
`Conda.add("matplotlib")`

- Add DrWatson to your default enviroment.

One method from the Julia REPL:\
`]` \
`add DrWatson`\
Backspace to leave package manager\
`;` enter shell mode \
`cd BASEDIR/ECCOonPoseidon` where BASEDIR is the location of the project
Backspace to leave shell \
`]` \
`activate .` \
`instantiate`

# Directory structure
- `scripts`: production-ready scripts
- `src`: definitions of poseidon-specific variables and functions

# Running scripts

Scripts are available for preprocessing, postprocessing, and scientific analysis. From a shell, use the following commands to run a script:\
`cd BASEDIR/ECCOonPoseidon/scripts`, where `BASEDIR` is where you cloned the repository.

Preprocess the forcing fields with:\
`julia filter_interannual.jl southpac`, one example where the script takes one argument (surface region).

Postprocess the experimental output with:\
`julia postprocess.jl nointerannual`, one example where the script takes one argument (experiment name).

- Preprocessing \
filter_interannual.jl 

- Postprocessing\
map2regularpolesDepth.jl \
mdsio2regularpoles.jl \
netcdf2regularpoles.jl \
regularpoles2sigma1.jl \
state2sigma1.jl 

- Scientific Analysis \
argotrends.jl \
experiment_divergence.jl \
plot_divergence.jl \
trends.jl \
trends_plot.jl

# Reproducibility

This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> ECCOonPoseidon

It is authored by G Jake Gebbie.

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia REPL and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.
