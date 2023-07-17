# map θ, S, p to sigma 1 surfaces.
# This is a wrapper routine to read files on poseidon.
# ggebbie, 1-Apr-2021
#julia --project=@. scripts/state2sigma1.jl 
include("../src/intro.jl")

using Revise # for interactive use
using MITgcmTools, MeshArrays, Statistics
using ECCOtour, ECCOonPoseidon
sig2dir(expt::String) = rundir(expt)*"sigma2/"

# using JLD2, Dierckx, Interpolations
include(srcdir("config_exp.jl"))
runpath,diagpath = listexperiments(exprootdir())

# DEFINE THE LIST OF SIGMA1 VALUES.
sig2grid = sigma2grid()

## specific for state
# the state_3d monthly-average diagnostic output
TSroot = "state_3d_set1" # 1: θ, 2: S
RProot = ("trsp_3d_set2","trsp_3d_set3") # uvel, vvel, gm psix, gm psiy, rhoanoma, phihyd

splorder = 3 # spline order

expt = "iter129_bulkformula"
# first filter for state_3d_set1
filelist = searchdir(diagpath[expt],TSroot)

# second filter for "data"
datafilelist  = filter(x -> occursin("data",x),filelist)

# make an output directory for each expteriment
pathout = sig2dir(expt)
!isdir(pathout) && mkdir(pathout)
nt = length(datafilelist)
siggrid = sigma2grid()
for (tt, datafile) in enumerate(datafilelist)
    
    #print timestamp
    year,month = timestamp_monthly_v4r4(tt)

    # eliminate suffix
    fileroots = Vector{String}()
    fileroot = rstrip(datafile,['.','d','a','t','a'])
    push!(fileroots,fileroot)
    for root in RProot
        push!(fileroots,root*fileroot[14:end]) # a better way than "14"?
    end
    # Read from filelist, map to sigma-1, write to file

    mdsio2sigma(diagpath[expt],pathout,fileroots,γ, pstdz, siggrid, 2000.0, "sigma2";
    splorder=splorder) 

end
