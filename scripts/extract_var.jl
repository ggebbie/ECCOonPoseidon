# this script will compile and save a basinwide average of OHC 
# for a particular basin (basin_name) 

include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise 
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools,
    PyPlot, JLD2, DrWatson, Statistics, JLD2
using .OHC_helper
include(srcdir("config_exp.jl"))

do_constant_density=true 
do_total=true 
output2file = true

(ϕ,λ) = latlonC(γ) #each tile has dims=(λ, ϕ)
area = readarea(γ)

runpath,diagpath = listexperiments(exprootdir())
shortnames = expnames()
marks = expsymbols()
nexp = length(shortnames) # number of experiments
fileroot = "trsp_3d_set2" #DFxE_TH (1), DFyE_TH(2), ADVx_TH(3), ADVy_TH(4)
filename = "ADVy_TH"
outdir = "ECCO_vars/"
# nc = 2 #THETA vert flux is in here explicit 
nc = 4 #THETA vert. flux id implicit
for (keys,values) in shortnames
    expname = keys
    # if value ∉ ["129ff", "noIA"]
    if values ∉ ["129ff", "noIA"]
        println(expname)
        nz = 50
        filelist = searchdir(diagpath[expname],fileroot) # 1st filter for state_3d_set1
        datafilelist  = filter(x -> occursin("data",x),filelist) # 2nd filter for "data"
        @time var_exp = extract_3d_var(γ, nc, diagpath, expname, datafilelist, nz)  
        @save datadir(outdir * filename*"_"*expname*".jld2") var_exp
    end
end

