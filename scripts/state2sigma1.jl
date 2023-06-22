# map θ, S, p to sigma 1 surfaces.
# This is a wrapper routine to read files on poseidon.
# ggebbie, 1-Apr-2021
#julia --project=@. scripts/state2sigma1.jl 
include("../src/intro.jl")

using Revise # for interactive use
using MITgcmTools, MeshArrays, Statistics
using ECCOtour, ECCOonPoseidon
# using JLD2, Dierckx, Interpolations

include(srcdir("config_exp.jl"))


runpath,diagpath = listexperiments(exprootdir())

# DEFINE THE LIST OF SIGMA1 VALUES.
sig1grid = sigma1grid()

## specific for state
# the state_3d monthly-average diagnostic output
TSroot = "state_3d_set1" # 1: θ, 2: S
RProot = ("trsp_3d_set1","state_3d_set2") # uvel, vvel, gm psix, gm psiy, rhoanoma, phihyd

splorder = 3 # spline order

expt = "noinitadjust"
# first filter for state_3d_set1
filelist = searchdir(diagpath[expt],TSroot)

# second filter for "data"
datafilelist  = filter(x -> occursin("data",x),filelist)

# make an output directory for each expteriment
pathout = sig1dir(expt)
!isdir(pathout) && mkdir(pathout)
nt = length(datafilelist)

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
    mdsio2sigma1(diagpath[expt],pathout,fileroots,γ,pstdz,sig1grid;splorder=splorder)

    fileprefix1 = pathout
    filesuffix1 = "_on_sigma1"*fileroots[1][14:end]*".data"
    test = joinpath(fileprefix1,"THETA"*filesuffix1)
    
end
