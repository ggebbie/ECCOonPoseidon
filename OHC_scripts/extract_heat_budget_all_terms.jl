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
uplvl = Inf; botlvl = -Inf
tt3km = findall( botlvl .<= z[:].<= uplvl)
lvls = tt3km
lvls = [10, 11, 12]
outdir = "ECCO_vars/"
savename = "Surface_STHETA"

filedir = "ECCO_vars/"
filename = "THETAs"
filename2 = "ETAN"

ocean_mask = wet_pts(Γ)
msk = ocean_mask
smush_depths =  Float32.(smush(cell_depths))
smush_depths[findall(smush_depths .== 0)] = Inf32
inv_depths = 1 ./ smush_depths

""" vertically integrate theta (should use eTan """
for (key,values) in shortnames
    expname = key; println(key)
    @time @load datadir(filedir*filename*"_"*expname*".jld2") var_exp; θ = var_exp
    @time @load datadir(filedir*filename2*"_"*expname*".jld2") var_exp; ETAN = var_exp
    var_exp  = Vector{MeshArrays.gcmarray{Float32, 1, Matrix{Float32}}}(undef, 0)
    temp = slice_mavec(θ, lvls, γ)
    @inbounds for tt in 1:length(temp)
        s1 = @views ((ETAN[tt] .* inv_depths) )
        s2 = (s1 .+ 1)
        s3 = @views s2 .* msk 
        temp2 = temp[tt] .* s3
        push!(var_exp, temp2)
    end
    @save datadir(outdir * savename*"_"*expname*".jld2") var_exp
    @time GC.gc(true)
end

""" compute volumes average heat budget"""
filedir = "ECCO_vars/"
filenameA = "THETA_AFLX_";
filenameE = "THETA_DFLXExplicit_";
filenameI = "THETA_DFLXImplicit_";
filenamesDHs = ["DFxE_TH_", "DFyE_TH_"]
filenamesAHs = ["ADVx_TH_", "ADVy_TH_"]

# expnamez = ["iter129_bulkformula", "iter0_bulkformula", "nosfcadjust", "noinitadjust"]
# expname = expnamez[1] #loops cause major memory leaks! 
# println(expname)

function diff_ma_vec(var_ts::Vector{MeshArrays.gcmarray{Float32, 2, Matrix{Float32}}}, 
    lvls::Vector{Int64}, γ::gcmgrid; tecco = tecco)
    temp_vec = Vector{MeshArrays.gcmarray{Float32, 2, Matrix{Float32}}}(undef, length(tecco))
    ma_vec = slice_mavec(var_ts, lvls, γ); ma_vecp1 = slice_mavec(var_ts, lvls.+1, γ);
    @inbounds for tt in 1:length(tecco)
        temp_vec[tt] = @views ma_vec[tt] .- ma_vecp1[tt]
    end
    return temp_vec
end

# function compute_budget(expname::String, lvls::Vector{Int64}, γ::gcmgrid, tecco::Vector{Float64}, 
#                         filedir::String, filenamesAHs::Vector{String}, filenamesDHs::Vector{String}, 
#                         filenameA::String, filenameE::String, filenameI::String) 
    #MeshArray(γ,Float32,50 * 312)
    expname = "iter129_bulkformula"
    outdir = "ECCO_vars/"
    savename = "THETA_BUDG"

    θ_budget = Dict()
    θ_budget["AdvH"]  = Vector{MeshArrays.gcmarray{Float32, 2, Matrix{Float32}}}(undef, length(tecco))
    θ_budget["AdvR"]  = Vector{MeshArrays.gcmarray{Float32, 2, Matrix{Float32}}}(undef, length(tecco))
    θ_budget["DiffH"] = Vector{MeshArrays.gcmarray{Float32, 2, Matrix{Float32}}}(undef, length(tecco))
    θ_budget["DiffZ"] = Vector{MeshArrays.gcmarray{Float32, 2, Matrix{Float32}}}(undef, length(tecco))

    Diffxpath = datadir(filedir*filenamesDHs[1]*expname*".jld2")
    Diffypath = datadir(filedir*filenamesDHs[2]*expname*".jld2")
    @time θ_budget["DiffH"] = OHC_helper.load_and_calc_UV_conv(Diffxpath, Diffypath, lvls, γ);

    Advxpath = datadir(filedir*filenamesAHs[1]*expname*".jld2")
    Advypath = datadir(filedir*filenamesAHs[2]*expname*".jld2")
    @time θ_budget["AdvH"] = OHC_helper.load_and_calc_UV_conv(Advxpath, Advypath, lvls, γ)

    @time var_exp = load_object_compress(datadir(filedir*filenameA*expname*".jld2")); 
    @time θ_budget["AdvR"] = diff_ma_vec(var_exp, lvls, γ)
    # AdvR = slice_mavec(var_exp, lvls, γ); AdvRp1 = slice_mavec(var_exp, lvls.+1, γ); var_exp = nothing
    @time var_exp .= load_object_compress(datadir(filedir*filenameE*expname*".jld2"));
    @time θ_budget["DiffZ"] = diff_ma_vec(var_exp, lvls, γ)
    # DifER = slice_mavec(var_exp, lvls, γ); DifERp1 = slice_mavec(var_exp, lvls.+1, γ); var_exp = nothing
    @time var_exp .= load_object_compress(datadir(filedir*filenameI*expname*".jld2"))
    @time θ_budget["DiffZ"] .+= diff_ma_vec(var_exp, lvls, γ)
    # DifIR = slice_mavec(var_exp, lvls, γ); DifIRp1 = slice_mavec(var_exp, lvls.+1, γ); var_exp = nothing
    @time GC.gc(true)

    @time @inbounds for tt in 1:length(tecco)
        θ_budget["AdvH"][tt] =  AdvH[tt] 
        θ_budget["DiffH"][tt] =  DiffH[tt]
        θ_budget["DiffZ"][tt] =  @views divDifER[tt] .+ divDifIR[tt]
        θ_budget["AdvR"][tt] =  divAdvR[tt]
    end
    var_exp = θ_budget
    println("Saving budget")
    @time save(datadir(outdir * savename *"_"*expname*"_test.jld2"), 
    Dict("var_exp" => var_exp))
    GC.gc(true)
end    

for (key,values) in shortnames
    expname = key
    @time compute_budget(expname, lvls, γ, tecco, 
    filedir, filenamesAHs, filenamesDHs, filenameA, 
    filenameE, filenameI);
end


#check budget 
# H = -DiffR .- AdvR .+ DifsH .+ AdvH
#piecuch uses -H
# println(H[1][6, 200] / cell_volumes[1, 30][6, 200]) #budget looks closed when normalized by volume! 

