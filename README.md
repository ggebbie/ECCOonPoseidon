# ECCOonPoseidon

Analysis of ECCO version 4 release 4 using output from Poseidon @ WHOI. A scientific project that relies heavily on the [ECCOtour](https://github.com/ggebbie/ECCOtour.jl), MITgcmTools, and MeshArrays Julia packages.

# Getting started

Clone with `git clone https://github.com/ggebbie/ECCOonPoseidon` \
Or `git clone git@github.com:ggebbie/ECCOonPoseidon.git` with SSH keys.

# Recommendations

- Julia 1.10

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

- `GeoPythonPlot` for plotting

`GeoPythonPlot.jl` uses matplotlib, cmocean, and cartopy. It should automatically build a Conda/Python environment within Julia for you. If you find this wasteful or duplicative, direct the python environment to an existing system-wide version of python with these already installed:
`ENV["PYTHON"]="python/directory/on/your/machine"` and rebuild GeoPythonPlot. 


# Directory structure
- `scripts`: production-ready scripts (may require updating)
- `src`: definitions of *poseidon-specific* variables and functions (other variables and functions belong in `ECCOtour.jl`)

# Running scripts

Scripts are available for preprocessing, postprocessing, and scientific analysis. From a shell, use the following commands to run a script:\
`cd BASEDIR/ECCOonPoseidon/scripts`, where `BASEDIR` is where you cloned the repository.

Preprocess the forcing fields with:\
`julia --project=@. filter_interannual.jl southpac`, one example where the script takes one argument (surface region).

Postprocess the experimental output with:\
`julia --project=@. postprocess.jl nointerannual`, one example where the script takes one argument (experiment name).

For long jobs, it is worth using nohup following:
`nohup julia --project=@. scripts/postprocess.jl nointerannual > postprocess_nointerannual_23nov2021.out &`

- Preprocessing \
filter_interannual.jl \
filter_interannual_basinmask.jl 

- Postprocessing\
state2sigmaregularpoles.jl \
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

# Workflow for updating code

1. check out main branch on ECCOonPoseidon.
2. activate and instantiate in Julia REPL
3. get julia to latest version using `juliaup update` (should be julia 1.10.2 at time of writing)
4. batch run with e.g., 36 threads: `~/.juliaup/bin/julia -t 36 --project=@. state2sigmaregularpoles.jl `

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
