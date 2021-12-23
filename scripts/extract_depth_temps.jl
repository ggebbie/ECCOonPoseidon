using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools,
    PyPlot, JLD2, DrWatson

include(srcdir("config_exp.jl"))


# ϕ = lat, λ = lon
(ϕ,λ) = latlonC(γ)
area = readarea(γ)

runpath,diagpath = listexperiments(exprootdir());

# abbreviations for each experiment for labels, etc.
shortnames = expnames()
marks = expsymbols()
nexp = length(shortnames) # number of experiments
 
fileroot = "state_3d_set1"
dryval = 0.0
iswet(x) = x != dryval

