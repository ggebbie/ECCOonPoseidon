include("../../src/intro.jl")

using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, 
Printf, LinearAlgebra, LaTeXStrings, PyCall, DataFrames, CSV, DSP
import PyPlot as plt
import NaNMath as nm

using CSV

PDO = CSV.read("/home/ameza/ECCOonPoseidon/Native_scripts/winds/annpdo.csv",  DataFrame,
header=false)

@pyimport matplotlib.animation as anim

cm = pyimport("cmocean.cm");
@pyimport seaborn as sns;

sns.set_theme(context = "notebook", style = "ticks",
              palette = sns.color_palette("colorblind"));
              
include(srcdir("config_exp.jl"))
runpath,diagpath = listexperiments(exprootdir());

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

region = "NPAC"; 
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not")
reg_mask = LLCcropC(PAC_msk,γ)
function extract_tau(expname, γ)
    filelist = searchdir(diagpath[expname],"state_2d_set1") # first filter for state_3d_set1
    τdatafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    nt = length(τdatafilelist)
    tauX = zeros(size(reg_mask)..., nt)
    tauY = zeros(size(reg_mask)..., nt)
    curlτ = zeros(size(reg_mask)..., nt)

    for tt = 1:nt
        println("year ",Int(floor(tecco[tt]))," month ",((tt-1)%12)+1)
        Tname = τdatafilelist[tt]
        
        τx, τy = extract_ocnTAU(diagpath, expname , Tname, γ)
        τcurl = MeshArrays.curl(τx, τy, Γ)
        curlτ[:, :, tt].= LLCcropC(τcurl, γ)
        τxC, τyC = velocity2center(τx, τy, Γ)
        τE, τN = rotate_uv(τxC, τyC, Γ)
        tauX[:, :, tt] .= LLCcropC(τE, γ)
        tauY[:, :, tt] .= LLCcropC(τN, γ)
    end
    return tauX, tauY, curlτ
end

τexps = Dict()
expname = "iter129_bulkformula"
τx_129, τy_129, curlτ129 = extract_tau(expname, γ);
τexps[expname] = (τx = τx_129, τy = τy_129, curlτ = curlτ129 .* 1)

expname = "iter0_bulkformula"
τx_0, τy_0, curlτ0 = extract_tau(expname, γ);
τexps[expname] = (τx = τx_0, τy = τy_0, curlτ = curlτ0 .* 1)

jldsave(datadir("native/oceτ_timeseries_regular_grid.jld2"), τ = τexps)
