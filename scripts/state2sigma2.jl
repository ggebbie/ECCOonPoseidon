# map θ, S, p to sigma 1 surfaces.
# This is a wrapper routine to read files on poseidon.
# ggebbie, 1-Apr-2021
#julia --project=@. scripts/state2sigma1.jl 
include("../src/intro.jl")

using Revise # for interactive use
using MITgcmTools, MeshArrays, Statistics
using ECCOtour, ECCOonPoseidon, PyCall

sig2dirs = Dict()
sig2dirs["iter0_bulkformula"] = vastrundir("iter0_bulkformula")*"sigma2/"
sig2dirs["iter129_bulkformula"] = vastrundir("iter129_bulkformula")*"sigma2/"

sig2dirs["only_init"] = vastrundir("nosfcadjust_exps", "run_adjust_init")*"sigma2/"
sig2dirs["only_kappa"] = vastrundir("nosfcadjust_exps", "run_adjust_kappa")*"sigma2/"
sig2dirs["only_sfc"] = vastrundir("nooceanadjust", "run_noadjusts")*"sigma2/"
sig2dirs["only_buoyancy"] = vastrundir("nooceanadjust", "run_noadjusts_nowind")*"sigma2/"
sig2dirs["only_wind"] = vastrundir("nooceanadjust", "run_noadjusts_nobuoyancy")*"sigma2/"

sig2dir(expt::String) = sig2dirs[expt]
# using JLD2, Dierckx, Interpolations
include(srcdir("config_exp.jl"))
runpath,diagpath = listexperiments(exprootdir())

# define the sigma grid you wish to interpolate onto
sig2grid = ECCOtour.sigma2grid()

## specific for state
# the state_3d monthly-average diagnostic output
TSroot = "state_3d_set1" # 1: θ, 2: S
# RProot = ("trsp_3d_set2","trsp_3d_set3") # uvel, vvel, gm psix, gm psiy, rhoanoma, phihyd

splorder = 3 # spline order

include(srcdir("plot_and_dir_config.jl"))
searchdir(diagpath["iter129_bulkformula"],TSroot)
for expt in keys(sig2dirs)
    # first filter for state_3d_set1
    filelist = searchdir(diagpath[expt],TSroot)

    # second filter for "data"
    datafilelist  = filter(x -> occursin("data",x),filelist)

    # make an output directory for each expteriment
    pathout = sig2dir(expt)
    !isdir(pathout) && mkdir(pathout)
    nt = length(datafilelist)
    siggrid = sig2grid

    for (tt, datafile) in enumerate(datafilelist)
        
        #print timestamp
        year,month = timestamp_monthly_v4r4(tt)

        # eliminate suffix
        fileroots = Vector{String}()
        fileroot = rstrip(datafile,['.','d','a','t','a'])
        push!(fileroots,fileroot)
        
        # for root in RProot
        #     push!(fileroots,root*fileroot[14:end]) # a better way than "14"?
        # end
        # Read from filelist, map to sigma-1, write to file

        mdsio2sigma(diagpath[expt],pathout,fileroots,γ, pstdz, siggrid, 2000.0, "sigma2";
        splorder=splorder) 

    end
end