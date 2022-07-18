
# depthlbl = string(abs(round(z[lvl_idx], digits = 2)))


# this script will compile and save a basinwide average of OHC 
# for a particular basin (basin_name) 

include("../src/intro.jl")
include("../src/OHC_helper.jl")
using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools,
PyPlot, JLD2, DrWatson, Statistics, JLD2,
GoogleDrive,NCDatasets, NetCDF, Printf
using .OHC_helper
using PyPlot   # important!
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
uplvl = -2000; botlvl = -3000
tt3km = findall( botlvl .<= z[:].<= uplvl)
lvls = tt3km
function slice_mavec(mavec, lvls)
    temp = Vector{MeshArrays.gcmarray{Float32, 2, Matrix{Float32}}}(undef, 0)
    for tt in 1:length(mavec)
        push!(temp, mavec[tt][:, lvls])
    end
    return temp
end

""" vertically integrate theta (should use eTan """
outdir = "ECCO_vars/"
savename = "Deep_THETA"

filedir = "ECCO_vars/"
filename = "THETAs"
ocean_mask = wet_pts(Γ)
msk = ocean_mask
for (keys,values) in shortnames
    expname = keys; println(keys)

    @load datadir(filedir*filename*"_"*expname*".jld2") var_exp; θ = var_exp
    temp = slice_mavec(var_exp, lvls)
    var_exp = temp 
    @save datadir(outdir * savename*"_"*expname*".jld2") var_exp
    GC.gc()
end

""" compute volumes average heat budget"""
outdir = "ECCO_vars/"
filedir = "ECCO_vars/"
filenameA = "THETA_AFLX_";
filenameE = "THETA_DFLXExplicit_";
filenameI = "THETA_DFLXImplicit_";
filenamesDHs = ["DFxE_TH_", "DFyE_TH_"]
filenamesAHs = ["ADVx_TH_", "ADVy_TH_"]
# expname = "noinitadjust"
# expname = "nosfcadjust"

savename = "THETA_BUDG"
for (key,values) in shortnames
    expname = key; println(expname)
    θ_budget = Dict()
    θ_budget["AdvH"]  = Vector{MeshArrays.gcmarray{Float32, 2, Matrix{Float32}}}(undef, 0)
    θ_budget["AdvR"]  = Vector{MeshArrays.gcmarray{Float32, 2, Matrix{Float32}}}(undef, 0)
    θ_budget["DiffH"] = Vector{MeshArrays.gcmarray{Float32, 2, Matrix{Float32}}}(undef, 0)
    θ_budget["DiffZ"] = Vector{MeshArrays.gcmarray{Float32, 2, Matrix{Float32}}}(undef, 0)
    @time @load datadir(filedir*filenamesDHs[1]*expname*".jld2") var_exp; 
    Difx = slice_mavec(var_exp, lvls); var_exp = nothing 
    @time @load datadir(filedir*filenamesDHs[2]*expname*".jld2") var_exp; 
    Dify = slice_mavec(var_exp, lvls); var_exp = nothing 
    println("Computing diffusion horizontal convergences...")
    @time DiffH = calc_UV_conv(Difx,  Dify); Difx = nothing; Dify = nothing;
    GC.gc()

    @time @load datadir(filedir*filenamesAHs[1]*expname*".jld2") var_exp; 
    Advx = slice_mavec(var_exp, lvls); var_exp = nothing 
    @time @load datadir(filedir*filenamesAHs[2]*expname*".jld2") var_exp; 
    Advy = slice_mavec(var_exp, lvls); var_exp = nothing 
    println("Computing advection horizontal convergences...")
    @time AdvH = calc_UV_conv(Advx,  Advy); Advx = nothing; Advy = nothing; 
    GC.gc()

    @time @load datadir(filedir*filenameA*expname*".jld2") var_exp; 
    AdvR = slice_mavec(var_exp, lvls); AdvRp1 = slice_mavec(var_exp, lvls.+1); var_exp = nothing
    @time @load datadir(filedir*filenameE*expname*".jld2") var_exp;
    DifER = slice_mavec(var_exp, lvls); DifERp1 = slice_mavec(var_exp, lvls.+1); var_exp = nothing
    @time @load datadir(filedir*filenameI*expname*".jld2") var_exp;
    DifIR = slice_mavec(var_exp, lvls); DifIRp1 = slice_mavec(var_exp, lvls.+1); var_exp = nothing
    GC.gc()

    @time for tt in 1:length(tecco)
        temp1 = @views (DifER[tt] .- DifERp1[tt]) 
        temp2 = @views (DifIR[tt] .- DifIRp1[tt]) 
        temp3 = @views temp1 .+ temp2
        temp4 = @views (AdvR[tt]  .- AdvRp1[tt]) 
        push!(θ_budget["AdvH"], AdvH[tt] )
        push!(θ_budget["DiffH"], DiffH[tt] )  
        push!(θ_budget["DiffZ"], temp3)
        push!(θ_budget["AdvR"], temp4)
    end
    # z = θ_budget["AdvH"][1] .+ θ_budget["DiffH"][1] .- θ_budget["DiffZ"][1] .- θ_budget["AdvR"][1]
    var_exp = θ_budget
    println("Saving budget")
    @time @save datadir(outdir * savename *"_"*expname*".jld2") var_exp
end


#check budget 
# H = -DiffR .- AdvR .+ DifsH .+ AdvH
#piecuch uses -H
println(H[1][6, 200] / cell_volumes[1, 30][6, 200]) #budget looks closed when normalized by volume! 
