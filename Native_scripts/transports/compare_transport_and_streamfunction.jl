#analysis should complete within 12-13 minutes 
#using 12 threads 
# julia --threads=6 --project=@. ./extract_heat_budget_native.jl

include("../../src/intro.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, 
    LaTeXStrings, PyCall
import NaNMath as nm
import PyPlot as plt

include(srcdir("config_exp.jl"))

cmo = pyimport("cmocean.cm");
@pyimport seaborn as sns;

colors =  sns.color_palette("deep")[1:4]
sns.set_theme(context = "talk", font_scale = 1.0,
              palette = sns.color_palette("deep"));
sns.set_style("white")

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

ocean_mask = wet_pts(Γ)
region = "PAC"; 
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "false", include_bering = true)

ϕ_avg = zonal_average(ϕ, area .* PAC_msk)
ϕ_avg = ϕ_avg[isfinite.(ϕ_avg)]

NPAC_boundidx = Base.findmin(abs.(ϕ_avg .- 23))[2]

tecco = 1992+1/24:1/12:2018

region = "NPAC"; 
uplvl = -2e3; botlvl = -3e3; suffix = "2to3"
lvls = findall(botlvl .< -z[:] .< uplvl)

expname = "iter129_bulkformula"

#load actual transports
fname = datadir("native/" * region * "_" * expname * "_TRSP" * ".jld2")
transports = load(fname)["transports"]
Vin = transports.Vin
Wtop = transports.Wtop
Wbot = transports.Wbot

fname = datadir("Ψ_EulBol_timeseries_PAC_" * expname *".jld2")
Ψ = load(fname)["Ψ_exp_timeseries"][:, NPAC_boundidx, :]
ΨVin = Ψ[lvls[end] + 1, :] .- Ψ[lvls[1], :]

fig, ax = plt.subplots()
ax.plot(tecco, Vin)
ax.plot(tecco, ΨVin)
fig

fname = datadir("Ψ_EulBol_timeseries_PAC_" * expname *".jld2")
Ψ = load(fname)["Ψ_exp_timeseries"][lvls[1], NPAC_boundidx, :]
ΨWtop =-Ψ 
# ΨWtop = -mean(Float32, Ψ, dims = 1)[1, :]
# Ψ = Ψ[2, :, :]
# Ψ = Ψ[1, :]
fig, ax = plt.subplots()
ax.plot(tecco, Wtop)
ax.plot(tecco, ΨWtop)
mean(Wtop .- ΨWtop)
fig

fname = datadir("Ψ_EulBol_timeseries_PAC_" * expname *".jld2")
Ψ = load(fname)["Ψ_exp_timeseries"][lvls[end]+1, NPAC_boundidx, :]
# ΨWbot = -mean(Float32, Ψ, dims = 1)[1, :, :]
ΨWbot = - Ψ
# Ψ = Ψ[2, :, :]
# Ψ = Ψ[1, :]
fig, ax = plt.subplots()
ax.plot(tecco, Wbot)
ax.plot(tecco, ΨWbot)
fig


fig, ax = plt.subplots()
ax.plot(tecco, 1e-6 .* (Wbot - Wtop + Vin), label = "actual")
ax.plot(tecco, 1e-6 .* (ΨWbot - ΨWtop + ΨVin), label = "approx")
ax.legend()
fig

obtain_Vin(Ψ, NPAC_boundidx) = Ψ[lvls[end] + 1, NPAC_boundidx, :] .- Ψ[lvls[1], NPAC_boundidx, :]
function obtain_Wtop(Ψ, NPAC_boundidx)
    Win = Ψ[lvls[1:2] .- 1, NPAC_boundidx, :]
    Win = mean(Win, dims = 1)[1, :]
    return -Win
end

function obtain_Wbot(Ψ, NPAC_boundidx)
    Wout = Ψ[lvls[end-1:end] .+ 1, NPAC_boundidx, :]
    Wout = mean(Win, dims = 1)[1, :]
    return -Wout
end
