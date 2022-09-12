#analysis should complete within 12-13 minutes 
#using 12 threads 
# julia --threads=6 --project=@. ./extract_heat_budget_native.jl

include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, BenchmarkTools
using .OHC_helper
using Plots
using ColorSchemes

include(srcdir("config_exp.jl"))


#using . is faster 
(ϕ,λ) = latlonC(γ)
area = readarea(γ)

runpath,diagpath = listexperiments(exprootdir());

# abbreviations for each experiment for labels, etc.
ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)
marks = expsymbols()
nexp = length(shortnames) # number of experiments
 
fileroot = "trsp_3d_set1"
dryval = 0.0
iswet(x) = x != dryval
ocean_mask = OHC_helper.wet_pts(Γ)
region = "PAC"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; region)

tecco = 1992+1/24:1/12:2018
#pre-allocate
nt = length(tecco)
function filter_heat_budget(diagpath::Dict{String, String}, expname::String, γ::gcmgrid)
    filelist = searchdir(diagpath[expname],"trsp_3d_set2") # first filter for state_3d_set1
    datafilelist_H  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    filelist = searchdir(diagpath[expname],"trsp_3d_set3") # first filter for state_3d_set1
    datafilelist_R  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    nt = length(datafilelist_R)
    Threads.@threads for tt = 1:nt
        # d = Dict{String, MeshArray}()
        fnameH = datafilelist_H[tt]
        fnameR = datafilelist_R[tt]
        HBUDGH = extract_θHbudget(expname,diagpath, γ, fnameH)
        HBUDGR = extract_θRbudget(expname,diagpath, γ, fnameR)
    end
    GC.safepoint()
end
@time for expname in keys(shortnames)
    println(expname)
    filter_heat_budget(diagpath, expname, γ)
end
