include("../src/intro.jl")
include("../src/OHC_helper.jl")
using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, JLD2,
NCDatasets, NetCDF, Printf
using .OHC_helper
using PyCall
@pyimport seaborn as sns
sns.set()
pygui(false)
cm = pyimport("cmocean.cm")
colorway = cm.balance;

include(srcdir("config_exp.jl"))
(ϕ,λ) = latlonC(γ)
area = readarea(γ)
tecco = collect(Float64, 1992+1/24:1/12:2018)

runpath,diagpath = listexperiments(exprootdir())
shortnames = expnames()
ignore_list = ["noIA", "129ff"]
shortnames = reduce_dict(shortnames, ignore_list)

marks = expsymbols()
nexp = length(shortnames) # number of experiments

ocean_mask = wet_pts(Γ)
msk = Float32.(ocean_mask)
cell_depths = get_cell_depths(msk, ΔzF, Γ.hFacC)
cell_depths_32 = Float32.(cell_depths)
cell_volumes = get_cell_volumes(area, cell_depths);
cell_volumes_Inf = Float32.(cell_depths)
cell_volumes_Inf[findall(cell_volumes .== 0)] = Inf32


ocean_mask = wet_pts(Γ)
msk = ocean_mask
smush_depths =  Float32.(smush(cell_depths))
smush_depths[findall(smush_depths .== 0)] = Inf32
inv_depths = 1 ./ smush_depths


uplvl = -1e3; botlvl = -4e3
lvls = findall( botlvl .<= z[:].<= uplvl)
suffix = "1to4"

""" compute volumes average heat budget"""
filedir = "ECCO_vars/"
filenameA = "THETA_AFLX_";
filenameE = "THETA_DFLXExplicit_";
filenameI = "THETA_DFLXImplicit_";
filenamesDHs = ["DFxE_TH_", "DFyE_TH_"]
filenamesAHs = ["ADVx_TH_", "ADVy_TH_"]


function compute_budget(expname::String, lvls::Vector{Int64}, γ::gcmgrid, tecco::Vector{Float64}, 
                        filedir::String, filenamesAHs::Vector{String}, filenamesDHs::Vector{String}, 
                        filenameA::String, filenameE::String, filenameI::String) 
    outdir = "ECCO_vars/"
    savename = "THETA_BUDG"
    nt = length(tecco)
    θ_budget = Dict()

    Diffxpath = datadir(filedir*filenamesDHs[1]*expname*".jld2")
    Diffypath = datadir(filedir*filenamesDHs[2]*expname*".jld2")
    @time θ_budget["DiffH"] = load_and_calc_UV_conv(Diffxpath, Diffypath);

    Advxpath = datadir(filedir*filenamesAHs[1]*expname*".jld2")
    Advypath = datadir(filedir*filenamesAHs[2]*expname*".jld2")
    @time θ_budget["AdvH"] = load_and_calc_UV_conv(Advxpath, Advypath)

    @time var_exp = load_object_compress(datadir(filedir*filenameA*expname*".jld2")); 
    @time θ_budget["AdvR"] = diff_ma_vec(var_exp, γ);
    @time var_exp .= load_object_compress(datadir(filedir*filenameE*expname*".jld2"));
    @time θ_budget["DiffZ"] = diff_ma_vec(var_exp, γ);
    @time var_exp .= load_object_compress(datadir(filedir*filenameI*expname*".jld2"))
    @time θ_budget["DiffZ"] .+= diff_ma_vec(var_exp, γ);

    var_exp = θ_budget
    println("Saving budget")
    @time save(datadir(outdir * savename *"_"*expname*"_full.jld2"), 
    Dict("var_exp" => var_exp))
    GC.gc(true)
end    

for expname in keys(shortnames)
    @time compute_budget(expname, lvls, γ, tecco, 
    filedir, filenamesAHs, filenamesDHs, filenameA, 
    filenameE, filenameI);
end
""" vertically integrate theta (should use eTan """


outdir = "ECCO_vars/"

filedir = "ECCO_vars/"
filename = "THETAs"
filename2 = "ETAN"

savename = "STHETA"
for (key,values) in shortnames
    expname = key; println(key)
    @time θ = load_object(datadir(filedir*filename*"_"*expname*".jld2")); 
    @time ETAN = load_object(datadir(filedir*filename2*"_"*expname*".jld2"));  
    # var_exp  = Vector{MeshArrays.gcmarray{Float32, 1, Matrix{Float32}}}(undef, 0)
    θ = slice_mavec(θ, lvls, γ)
    for tt in 1:length(θ)
        s1 = @views ((ETAN[tt] .* inv_depths) )
        s1 .+= 1
        s1 .*= msk 
        sθ = @views θ[tt] .* s1
        jldsave(datadir(filedir * expname * "/" * 
                      savename*"_" * suffix * "_"*
                      string(tt)*"_.jld2"),
                      sθ = sθ)
    end
    
    @time GC.gc(true)
end

