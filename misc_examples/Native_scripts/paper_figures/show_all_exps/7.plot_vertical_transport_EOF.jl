include("../../../src/intro.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson,
    BenchmarkTools, LaTeXStrings, PyCall, DataFrames
import NaNMath as nm
import PyPlot as plt
import NumericalIntegration
using LinearAlgebra

cumul_integrate = NumericalIntegration.cumul_integrate
include(srcdir("config_exp.jl"))

cmo = pyimport("cmocean.cm");

#blue, red, #green, orangeDataFrames

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

ocean_mask = wet_pts(Γ)

region = "PAC"; 
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "false", include_bering = true)

ϕ_avg = zonal_average(ϕ, area .* PAC_msk)
ϕ_avg = ϕ_avg[isfinite.(ϕ_avg)]

open_streamfunction(expname, type) = jldopen(datadir("Ψ_" *type * "_timeseries_PAC_" * expname *".jld2"))["Ψ_exp_timeseries"]

Ψ_ϕ = jldopen(datadir("Ψ_EulBol_timeseries_PAC_iter0_bulkformula.jld2"))["ϕ_avg"]

lvls = findall( -3000 .<= -z[:].<= -2000)
uplvl = -2e3; botlvl = -3e3; suffix = "2to3"
lvls = findall( botlvl .<= z[:].<= uplvl)
tecco = 1992+1/24:1/12:2018

include(srcdir("plot_and_dir_config.jl"))
sns.set_theme(context = "paper", style = "ticks",
              palette = colors, rc = custom_params);

Ψ_init = open_streamfunction("only_init", "EulBol")
Ψ_0 = open_streamfunction("iter0_bulkformula", "EulBol")
ΔΨ = Ψ_init .- Ψ_0
Wres = -(ΔΨ[20:end, 1:end-2, :] .- ΔΨ[20:end, 3:end, :])
Wres = Wres[:, 47:end, :]
nZ, nY, nT = size(Wres) 
Y = reshape(Wres, (nZ * nY, nT))
Y = Y .- mean(Y, dims = 2)
Y[isnan.(Y)] .= 0.0

s = svd(Y; full = true)
iPC = 4
PC1 = (s.U[:, 1:iPC]' * Y)'
# PC1 = (PC1 .- mean(PC1, dims = 1)) ./ std(PC1, dims = 1)
cumsum(s.S.^2 ./ sum(s.S.^2))[1:7]

objective = 1 .* Wres
objective = objective .- mean(objective, dims = 3)

reconstruction = zeros(size(objective))
i = 1; j = 1

Rsq = zeros(size(objective)[1:2])
for i in 1:size(reconstruction, 1), j in 1:size(reconstruction, 2)

    df = DataFrame(PC1)
    df[!, "y"] = objective[i, j, :] .* 1e10
    ols = lm(@formula(y ~ x1 + x2 + x3 + x4 ), df)
    
    Rsq[i, j] = round(r2(ols); digits=5)

end

fig = plt.figure(constrained_layout=true, figsize = (8, 4))
gs = fig.add_gridspec(1, 6)
obj(i,j) = get(gs, (i,j))
slice(i,j) = pycall(pybuiltin("slice"), PyObject, i,j)
f3_ax1 = fig.add_subplot(obj(0, slice(0, 4)))
f3_ax2 = fig.add_subplot(obj(0, slice(4, 6)))

[f3_ax2.plot(tecco, PC1[:, i] ./ std(PC1[:, i]), label = "PC " * string(i)) for i = 1:iPC]
f3_ax2.legend(frameon=false)



CM = f3_ax1.pcolormesh(Ψ_ϕ[48:end-1], z[20:end],  100 .* Rsq, cmap = cmo.thermal, vmin = 0, vmax = 100)
f3_ax1.invert_yaxis()
f3_ax1.set_title("Variance Explained by First 5 PCs")
fig.colorbar(CM, orientation = "horizontal", fraction = 0.05)
f3_ax1.set_facecolor("black")
f3_ax1.set_xticks(0:20:60)
f3_ax1.set_xlim(0, 60)
lab = string.(abs.(collect(0:20:60)))
lab = lab .* ["", "°N", "°N", "°N"]
f3_ax1.set_xticklabels(lab)
fig

nm.mean(Rsq[24:28, :])