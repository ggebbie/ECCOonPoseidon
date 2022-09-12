include("../src/intro.jl")
include("../src/OHC_helper.jl")
using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, JLD2,
NCDatasets, Printf, 
DataFrames, LaTeXStrings,
Plots
gr()
import NaNMath as nm
using .OHC_helper
import PyPlot as plt
using PyCall
@pyimport seaborn as sns;
@pyimport pandas as pd;
colors =  sns.color_palette("deep")[1:4]
sns.set_theme(context = "poster", font_scale = 1.5,
              palette = sns.color_palette("deep"));#sns.set_context("talk")
pygui(false)
cm = pyimport("cmocean.cm");colorway = cm.balance;
include(srcdir("config_exp.jl"))
(ϕ,λ) = latlonC(γ)
area = readarea(γ)
const tecco = collect(Float64, 1992+1/24:1/12:2018)
runpath,diagpath = listexperiments(exprootdir())
ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)
labels_L = [L"\theta^{129}", L"\theta^{\Delta F}",L"\theta^{\Delta T}", L"\theta^{0}"]


suffix = "2to3"
uplvl = -2e3; botlvl = -3e3
lvls = findall( botlvl .<= z[:].<= uplvl)
region = "NPAC"; 
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; region)
ocean_mask = wet_pts(Γ)
msk = PAC_msk;
cell_depths = get_cell_depths(msk, ΔzF, Γ.hFacC)
cell_volumes = get_cell_volumes(area, cell_depths);

""" get the bottom heating """
GTF1 = get_geothermalheating(γ)
GTF = 1.82e-12 #derived by guess and check 
"""Compare vertical convergence of pot. temperature"""
filedir = "ECCO_vars/"
# fname1 = "Deep_THETA"

fname1 = "STHETA"
fname2 = "THETA_BUDG"

θC, θCAR, θCAH, θCDH, θCDR, θz = Dict(), Dict(), Dict(), Dict(), Dict(), Dict()
θ_depths, θ_zonal,θ_budget_deltas = Dict(), Dict(), Dict()
θz_std = Dict()

budg_names = ["AdvR", "AdvH", "DiffH", "DiffR"]
crop_vols = cell_volumes[:, lvls]
tot_vol = Float32(sum(cell_volumes[:, lvls]))
vol_weight(x) = Float32(sum(x .* msk) / tot_vol)
#x =  h5open(datadir(filedir*fname2*"_"*expname *".jld2"), "r")

function filter_heat_budget_terms!(expname::String, lvls::Vector{Int64},
                            crop_vols::MeshArrays.gcmarray{T, 2, Matrix{T}},
                            nt :: Int64, GTF :: Float64,
                            θz::Dict, θC::Dict, θCAH::Dict, θCAR::Dict, 
                            θCDH::Dict, θCDR::Dict, θ_depths::Dict, θ_zonal::Dict) where T<:Real
    nlev = length(lvls)
    for var in [θz, θz_std,θC, θCAH, θCAR, θCDH, θCDR]
        var[expname] = Float32[]
    end
    θ_depths[expname] = zeros(Float32, nlev, nt)
    θ_zonal[expname] = zeros(nlev, 270, nt)
    θ_contribution_TS = Dict(var => zeros(nlev, nt) 
                            for var in budg_names)

    println("constructing time series")
    @time for tt in 1:nt
        fdir1 = filedir * expname * "/" * fname1 * "_" * suffix * "_"*"$tt" *"_.jld2"
        fdir2 = filedir * expname * "/" * fname2 * "_" * suffix * "_"*"$tt" *".jld2"
        jldopen(datadir(fdir2), "r") do HBUDG
            AdvH  = vol_weight(HBUDG["AdvH"]);push!(θCAH[expname],AdvH)
            AdvR  = vol_weight(-HBUDG["AdvR"]);push!(θCAR[expname],AdvR)
            DiffH = vol_weight(HBUDG["DiffH"]);push!(θCDH[expname],DiffH)
            DiffR = vol_weight(-HBUDG["DiffR"]);push!(θCDR[expname],DiffR)
            avg1 = AdvH + AdvR + DiffH + DiffR + GTF
            push!(θC[expname], Float32.(avg1))
            level_timeseries!(θ_contribution_TS, HBUDG, crop_vols,
                                         lvls, msk, tt)
            
        end
        sθ = load_object(datadir(fdir1))
        θ̄ = volume_mean(sθ; weights = crop_vols)
        # σ = sqrt(volume_mean(variance; weights = crop_vols))
        # push!(θz_std[expname], Float32.(σ))
        push!(θz[expname], Float32.(θ̄))
        θ_depths[expname][:, tt] .= ma_horiz_avg(sθ, crop_vols)
        θ_zonal[expname][:, :, tt] .= ma_zonal_avg(sθ, crop_vols)
    end
    fdir = filedir * expname * "/" * fname1 * "_" * suffix * "_1_.jld2"
    θ1 = load_object(datadir(fdir))
    θ_budget_deltas[expname] = total_level_change(θ_contribution_TS, θ1, crop_vols, lvls)
    @time GC.gc() #garbage collecting 

end

for (key,values) in shortnames
    expname = key
    nt = length(tecco)
    filter_heat_budget_terms!(expname, lvls,
    crop_vols, nt, GTF,
    θz, θC, θCAH, θCAR, θCDH, θCDR, θ_depths, θ_zonal)
end

DFtoDict(x) = Dict(name => x[!, name] for name in names(x))
PDtoDict(x) = Dict(name => x[name].values for name in x.columns)
DFtoPD(x) = pd.DataFrame(Dict(name => x[!, name] for name in names(x)))
PDtoDF(x) = DataFrame(Dict(name => x[name].values for name in x.columns))

#=Plots heat budget terms (advection, diffusion) statistics=#
@time include("plot_divergence_hb_residual.jl")
#=Reconstructs temperature time series from the heat budget terms=#
@time include("plot_divergence_hb_reconstruction.jl")
#=Plots Δθ for deep (however defined) and surface. Also plots these values against each other
The question here is there a relationship between deep and surface values 
In future, surface should be taken as surface to 2km
=#
@time include("plot_divergence_hb_deltas.jl")

#=Depth plot for heat gain from individual heat budget terms 
=#
# @time include("plot_divergence_depth.jl")
# #=Same as above but looking at anomalies about iter0
# =#
@time include("plot_divergence_depth_deltas.jl")
@time include("plot_divergence_depth_deltas_anomaly.jl")

"""Show
"""