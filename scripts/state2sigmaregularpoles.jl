# map θ, S, p to sigma 1 surfaces.
# This is a wrapper routine to read files on poseidon.
# ggebbie, 1-Apr-2021

include("../src/intro.jl")

using Revise # for interactive use
using MITgcmTools
using MeshArrays
using Statistics
using ECCOtour
using ECCOonPoseidon
# using JLD2, Dierckx, Interpolations

include(srcdir("config_exp.jl"))
include(srcdir("config_regularpoles.jl"))

# DEFINE THE LIST OF SIGMA1 VALUES.
# 120 sigma values
sig1grid = sigma1grid("mixed layer")

## specific for state
# the state_3d monthly-average diagnostic output
TSroot = "state_3d_set1" # 1: θ, 2: S
RProot = ("trsp_3d_set1","state_3d_set2") # uvel, vvel, gm psix, gm psiy, rhoanoma, phihyd

splorder = 1 # spline order
eos_mitgcm = "JMD95"

# first filter for state_3d_set1
filelist = searchdir(diagpath,TSroot)

# second filter for "data"
datafilelist  = filter(x -> occursin("data",x),filelist)

# make an output directory for each expteriment
!isdir(path_out) && mkdir(path_out)

# attributes of sigma-1 vertical coordinate
# for writing to NetCDF
sigmaatts = Dict("longname" => "Sigma-1", "units" => "kg/m^3 - 1000")

tt = 0
Threads.@threads for datafile in datafilelist
    global tt += 1

    #print timestamp
    year,month = timestamp_monthly_v4r4(tt)
    if month < 10
        filesuffix = "_on_sigma1_"*string(year)*"_0"*string(month)*".nc"
    else 
        filesuffix = "_on_sigma1_"*string(year)*"_"*string(month)*".nc"
    end

    # eliminate suffix
    fileroots = Vector{String}()
    fileroot = rstrip(datafile,['.','d','a','t','a'])
    push!(fileroots,fileroot)

    # put all file roots into a tuple
    for root in RProot
        push!(fileroots,root*fileroot[14:end]) # a better way than "14"?
    end
    
    # Read from filelist, map to sigma-1, write to file
    @time varsσ = mdsio2sigma1(diagpath,
        path_out,
        fileroots,
        γ,
        pstdz,
        sig1grid,
        splorder=splorder,
        linearinterp=true,
        eos=eos_mitgcm,
        writefiles = false)

    # transfer to regularpoles grid
    @time varsσregpoles = regularpoles(varsσ,γ,rp_params)

    # write to NetCDF
    @time writeregularpoles(varsσregpoles,γ,
        pathout, filesuffix, filelog,
        rp_params,
        gridatts)
        # filelog, rp_params.λC,
        # gridatts.lon, rp_params.ϕC,
        # gridatts.lat,
        # sig1grid, sigmaatts)

end
