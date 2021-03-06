# map regularpoles output to sigma 1 surfaces.
# ggebbie, 25-May-2021

include("../src/intro.jl")

using Revise 
using MITgcmTools, MeshArrays, Statistics
using ECCOtour, ECCOonPoseidon #,JLD2 , Dierckx, Interpolations

include(srcdir("config_exp.jl"))

# the state_3d monthly-average diagnostic output on regularpoles grid
# some are not available
#varroot = ("EVEL","EVELMASS","EVELSTAR","NVEL","NVELMASS","NVELSTAR","PHIHYD","RHOAnoma","SALT","THETA")
varroot = ("PHIHYD","RHOAnoma","SALT","THETA")

sig1grid = sigma1grid()

# if splorder is larger than the number of points in a profile,
# then it will default to linear interpolation
#splorder = 3 # spline order
splorder = 100 # spline order

# first filter for state_3d_set1
regpolesroot = regpolesdir(expt)
!isdir(regpolesroot) ?  mkpath(regpolesroot) : nothing

# get specific file names, one for each variable
fileroots = Dict{String,String}()

for tt = 1:312 # must be a better way
    for fldname in varroot
        push!(fileroots,fldname => regpolesroot*fldname*"/"*fldname)
    end

    #print timestamp
    year,month = timestamp_monthly_v4r4(tt)
    if month < 10
        tstamp = "_"*string(year)*"_0"*string(month)*".nc"
    else
        tstamp = "_"*string(year)*"_"*string(month)*".nc"
    end

    # add the tstamp
    ncfilenames = fileroots
    for (kk,vv) in ncfilenames
        ncfilenames[kk] = vv.*tstamp
    end

    # Read from filelist, map to sigma-1, write to file
    netcdf2sigma1(regpolesroot,regpolesroot,ncfilenames,γ,pstdz,sig1grid,splorder)
end
