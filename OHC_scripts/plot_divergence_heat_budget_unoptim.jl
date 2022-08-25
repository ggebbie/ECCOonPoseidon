include("../src/intro.jl")
include("../src/OHC_helper.jl")
using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, JLD2,
NCDatasets, Printf, 
DataFrames, LaTeXStrings
import NaNMath as nm
using .OHC_helper
import PyPlot as plt
using PyCall
using Plots; gr()
using Plots.Measures

@pyimport seaborn as sns
@pyimport pandas as pd
# colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]
colors =  sns.color_palette("deep")[1:4]
sns.set_theme(context = "talk", palette = sns.color_palette("deep"));#sns.set_context("talk")
pygui(false)
cm = pyimport("cmocean.cm");colorway = cm.balance;
#fid = h5open(fname, "r")
include(srcdir("config_exp.jl"))
(ϕ,λ) = latlonC(γ)
area = readarea(γ)
const tecco = collect(Float64, 1992+1/24:1/12:2018)
runpath,diagpath = listexperiments(exprootdir())
shortnames = expnames()
ignore_list = ["noIA", "129ff"]
shortnames = reduce_dict(shortnames, ignore_list)

marks = expsymbols()
nexp = length(shortnames) # number of experiments

uplvl = -2000; botlvl = -3000
tt3km = findall( botlvl .<= z[:].<= uplvl)
region = "NPAC"; 
ocean_mask = wet_pts(Γ)
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; region)
const msk = PAC_msk;
cell_depths = get_cell_depths(msk, ΔzF, Γ.hFacC)
cell_volumes = get_cell_volumes(area, cell_depths);

""" get the bottom heating """
lvls = tt3km
GTF1 = get_geothermalheating(γ)
GTF = 1.82e-12 #derived by guess and check 
"""Compare vertical convergence of pot. temperature"""
filedir = "ECCO_vars/"
# fname1 = "Deep_THETA"

suffix = ""
fname1 = "Deep_STHETA"
fname2 = "THETA_BUDG"

θC, θCAR, θCAH, θCDH, θCDR, θz = Dict(), Dict(), Dict(), Dict(), Dict(), Dict()
θ_depths, θ_zonal = Dict(), Dict(), Dict()
budg_names = ["AdvR", "AdvH", "DiffH", "DiffR"]
lvls = tt3km
crop_vols = cell_volumes[:, lvls]
tot_vol = Float32(sum(cell_volumes[:, lvls]))
vol_weight(x) = Float32(sum(x .* msk) / tot_vol)
#x =  h5open(datadir(filedir*fname2*"_"*expname *".jld2"), "r")
E,F = trend_matrices(tecco)

function filter_heat_budget_terms!(expname::String, lvls::Vector{Int64},
                            crop_vols::MeshArrays.gcmarray{T, 2, Matrix{T}},
                            nt :: Int64, GTF :: Float64,
                            θz::Dict, θC::Dict, θCAH::Dict, θCAR::Dict, 
                            θCDH::Dict, θCDR::Dict, θ_depths::Dict, θ_zonal::Dict,
                            θ_budget_deltas::Dict) where T<:Real
    for var in [θz, θC, θCAH, θCAR, θCDH, θCDR]
        var[expname] = Float32[]
    end
    θ_depths[expname] = zeros(Float32, length(lvls), length(tecco))
    θ_zonal[expname] = zeros(length(lvls), 270, length(tecco))

    @time sθ = load_object(datadir(filedir*fname1*"_"*expname * ".jld2"))
    @time HBUDG = load_object(datadir(filedir*fname2*"_"*expname * ".jld2"))
    θ_contribution_TS = level_timeseries(HBUDG, sθ, crop_vols, msk, lvls)
    θ_budget_deltas[expname] = total_level_change(θ_contribution_TS, tecco, F)

    println("constructing time series")
    @time for tt in 1:nt
            AdvH  = vol_weight(HBUDG["AdvH"][tt]);push!(θCAH[expname],AdvH)
            AdvR  = vol_weight(-HBUDG["AdvR"][tt]);push!(θCAR[expname],AdvR)
            DiffH = vol_weight(HBUDG["DiffH"][tt]);push!(θCDH[expname],DiffH)
            DiffR = vol_weight(-HBUDG["DiffZ"][tt]);push!(θCDR[expname],DiffR)
            avg1 = (AdvH + AdvR + DiffH + DiffR + GTF) 
            avg2 = volume_mean(sθ[tt]; weights = crop_vols)
            push!(θC[expname], Float32.(avg1))
            push!(θz[expname], Float32.(avg2))
            θ_depths[expname][:, tt] .= ma_horiz_avg(sθ[tt], crop_vols)
            θ_zonal[expname][:, :, tt] .= ma_zonal_avg(sθ[tt], crop_vols)
    end
    sθ = nothing 
    HBUDG = nothing 
    @time GC.gc() #garbage collecting 
end
θ_budget_deltas = Dict()
for (key,values) in shortnames
    expname = key
    nt = length(tecco)
    filter_heat_budget_terms!(expname, lvls,
    crop_vols, nt, GTF,
    θz, θC, θCAH, θCAR, θCDH, θCDR, θ_depths, θ_zonal,
    θ_budget_deltas)
end
# for var in [:θC, :θCAR, :θCAH, :θCDH, :θCDR, :θz]
#     eval(:($var = pd.DataFrame($var, index = tecco).rename_axis("time") ))
# end
#convert DataFrame to 
xDFtoDict(x) = Dict(name => x[!, name] for name in names(x))
PDtoDict(x) = Dict(name => x[name].values for name in x.columns)
DFtoPD(x) = pd.DataFrame(Dict(name => x[!, name] for name in names(x)))
PDtoDF(x) = DataFrame(Dict(name => x[name].values for name in x.columns))

# for var in [:θC, :θCAR, :θCAH, :θCDH, :θCDR, :θz]
#     eval(:($var = PDtoDict($var)))
# end
#:θ_depths, :θ_zonal are 2-d arrays, so they will not work w/ DataFrames 
@time include("plot_divergence_hb_zonal.jl")
@time include("plot_divergence_hb_verticaltrends.jl")
@time include("plot_divergence_hb_verticaltrends_diff.jl")
@time include("plot_divergence_hb_residual.jl")
@time include("plot_divergence_hb_reconstruction.jl")
@time include("plot_divergence_hb_reconstruction_smoothed.jl")
@time include("plot_divergence_hb_deltas.jl")
@time include("plot_divergence_LinReg.jl")

@time include("plot_divergence_hb_clim.jl")
# @time include("plot_divergence_depth.jl")