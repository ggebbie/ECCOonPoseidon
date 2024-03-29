#  First steps: 1. go into Reemergence project directory. 2. go into julia REPL package mode with `]`. 3. `activate .` 4. Backspace to return to command mode in REPL.
# Modeled after "filter_interannual.jl"
# Solve for linear trends at all points.

<<<<<<< HEAD
include("../../src/intro.jl")

using Revise
using ECCOtour
using ECCOonPoseidon
using MeshArrays
using MITgcmTools
using Statistics
using Distributions
using LinearAlgebra
using StatsBase
=======
include("../src/intro.jl")

using Revise
using ECCOtour, ECCOonPoseidon,  MeshArrays, MITgcmTools
using Statistics, PyPlot, Distributions, LinearAlgebra, StatsBase
>>>>>>> origin/main

include(srcdir("config_exp.jl"))

# list of experiments on poseidon
runpath,diagpath = listexperiments(exprootdir());

# abbreviations for each experiment for labels, etc.
shortnames = expnames()

## SELECT EXPERIMENTS TO ANALYZE #################################
# manually choose from available experiments listed above.
#exps = ("iter129_bulkformula","nointerannual")

# to do all experiments:
exps = keys(shortnames)
#################################################################

# assume monthly averages, ECCOv4r4
tstart = 1992 + 1/24
tend = 2018
tecco = range(tstart,step=1/12,stop=2018)
nt = length(tecco)

<<<<<<< HEAD
# eliminate spinup years
tstart = 2004 + 1/24
tend = 2018
teccocrop = range(tstart,step=1/12,stop=2018)
ntcrop = length(teccocrop)

# get weight matrix for finding trends
E,F = trend_matrices(teccocrop)
Ecrop = vcat(zeros(nt - ntcrop,size(E,2)),E)
Fcrop = (Ecrop'*Ecrop)\Ecrop' # least squares estimator
=======
# get weight matrix for finding trends
E,F = trend_matrices(tecco)
>>>>>>> origin/main

# pre-allocate β, linear trends
nz = length(z)

# cycle through all chosen experiments
for expt in exps

<<<<<<< HEAD
    local β = trend_theta(diagpath[expt],tecco,γ,Fcrop)
=======
    β = trend_theta(diagpath[expt],tecco,γ,F)
>>>>>>> origin/main
#        trend_theta!(β,diagpath[expt],tecco,γ,F)
    println(maximum(β))
    # save β for each experiment
    # make an output directory for each experiment
    outdir = datadir("trends",expt)
    !isdir(outdir) ? mkpath(outdir) : nothing;
    Toutname = datadir(outdir,"DthetaDt.data")
    # save to file before overwritten next time step.
    γ.write(Toutname,β)
end
