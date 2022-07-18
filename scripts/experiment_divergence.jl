# interactive analysis of divergence of runs, originally 10-Feb-2021
# recreate MATLAB difference plots
# use files directly on poseidon through cifs mount

include("../src/intro.jl")

using Revise 
using MITgcmTools, MeshArrays, ECCOonPoseidon, ECCOtour
using Statistics, JLD2

# save output to file?
output2file = true

include(srcdir("config_exp.jl"))

# list of experiments on poseidon
runpath,diagpath = listexperiments(exprootdir())

# abbreviations for each experiment for labels, etc.
shortnames = expnames()

## SELECT EXPERIMENTS TO COMPARE #################################
# baseline experiment
expbase = "iter129_bulkformula"
# comparison experiment(s)
#expcompare = "iter129_fluxforced"
#expcompare = "iter0_bulkformula"
#expcompare = "nointerannual"
#expcompare = "nosfcadjust"
expcompare = "noinitadjust"
#################################################################

# Improve code here to compare to multiple experiments.
#nexp = length(expcompare) # number of experiments

# name of file inside diagspath
fileroot = "state_3d_set1"
# first filter for files corresponding to state_3d_set1
filelist = searchdir(diagpath[expbase],fileroot) 
 # second filter for files corresponding to the "data" of state_3d_set1
datafilelist  = filter(x -> occursin("data",x),filelist)
nt = length(datafilelist)

# Improve code here to read meta file, make variable selection transparent.
nc = 2 # do first two of three properties in state_3d_set1

# pre-allocate
nz = length(z)
σx = Array{Float32, 2}(undef, nt, nz*nc)
# xmed = similar(σ); median was slow so I dropped it.
xbar = similar(σx); xmax = similar(σx); xmin = similar(σx); absxbar = similar(xbar);

tt = 0
for fname in datafilelist
    global tt += 1

    #print timestamp
    year,month = timestamp_monthly_v4r4(tt)
    # need to get all three tracers read.
    @time xbase = γ.read(diagpath[expbase]*fname,MeshArray(γ,Float32,nz*nc))
    @time x = γ.read(diagpath[expcompare]*fname,MeshArray(γ,Float32,nz*nc))
    print(size(xbase))
    print(size(xbase[1]))
    x = x - xbase # turn into a difference 
    for zc = 1:nz*nc
        xbar[tt,zc],xmax[tt,zc],xmin[tt,zc],σx[tt,zc],absxbar[tt,zc] = faststats(x[:,zc])
    end
end

# save output as JLD2
if output2file
    outfile = datadir("faststats_"*shortnames[expbase]*"_vs_"*shortnames[expcompare]*".jld2")
    @save outfile absxbar xbar σx xmax xmin z 
end

# Call this to plot the results
include("plot_divergence.jl")
