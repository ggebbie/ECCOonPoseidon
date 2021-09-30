# map θ, S, p to sigma 1 surfaces.
# This is a wrapper routine to read files on poseidon.
# ggebbie, 1-Apr-2021

include("../src/intro.jl")

using Revise # for interactive use
using MITgcmTools, MeshArrays, Statistics
using ECCOtour, ECCOonPoseidon
# using JLD2, Dierckx, Interpolations

include(srcdir("config_exp.jl"))

# DEFINE THE LIST OF SIGMA1 VALUES.
sig1grid = sigma1grid()

## specific for state
# the state_3d monthly-average diagnostic output
TSroot = "state_3d_set1" # 1: θ, 2: S
RProot = "state_3d_set2" # 1:rhoanoma, 2 phihyd

splorder = 100 # spline order

# first filter for state_3d_set1
filelist = searchdir(diagpath,TSroot)

# second filter for "data"
datafilelist  = filter(x -> occursin("data",x),filelist)

filelist2 = searchdir(diagpath,RProot) 
datafilelist2  = filter(x -> occursin("data",x),filelist)

# make an output directory for each expteriment
!isdir(path_out) ? mkdir(path_out) : nothing;
nt = length(datafilelist)
    
global tt = 0
for datafile in datafilelist
    tt += 1

    #print timestamp
    year,month = timestamp_monthly_v4r4(tt)

    # eliminate suffix
    fileroot = rstrip(datafile,['.','d','a','t','a'])
    fileroot2 = RProot*fileroot[14:end] # a better way?
    fileroots = (fileroot,fileroot2)
    
    # Read from filelist, map to sigma-1, write to file
    mdsio2sigma1(diagpath,path_out,fileroots,γ,pstdz,sig1grid,splorder)
end
