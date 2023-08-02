# map θ, S, p to sigma 2 surfaces.
# ggebbie, 1-Apr-2021
#julia --project=@. scripts/state2sigma1.jl 
include("../src/intro.jl")
include("../src/OHC_helper.jl")
using Revise # for interactive use
using MITgcmTools, MeshArrays, Statistics
using ECCOtour, ECCOonPoseidon
using .OHC_helper
include(srcdir("config_exp.jl"))

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

area = readarea(γ)
ocean_mask = OHC_helper.wet_pts(Γ)
cell_depths = OHC_helper.get_cell_depths(ocean_mask, ΔzF, Γ.hFacC); 
cell_volumes = Float32.(OHC_helper.get_cell_volumes(area, cell_depths));

inv_RAU = cell_depths .* Γ.DYG #area of "u-face"
inv_RAU = inv_ma(inv_RAU .* Γ.hFacW)

inv_RAV = cell_depths .* Γ.DXG #area of "v-face"
inv_RAV = inv_ma(inv_RAV .* Γ.hFacS)

inv_RAC = inv_ma(deepcopy(Γ.RAC)) #area of "v-face"

expt = "iter129_bulkformula"
# first filter for state_3d_set1
filelist = searchdir(diagpath[expt],TSroot)

# second filter for "data"
datafilelist  = filter(x -> occursin("data",x),filelist)

# make an output directory for each expteriment
pathout = sig2dir(expt)
!isdir(pathout) && mkdir(pathout)
nt = length(datafilelist)

for (tt, datafile) in enumerate(datafilelist[1:30])
    
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
    @time θbudg2sigma2(inv_RAC, inv_RAU, inv_RAV, diagpath[expt],pathout,fileroots,γ,pstdz,sig2grid;splorder=splorder)
    #outputs ADV_TH DFE_TH DFr_TH ADVr_TH
end
