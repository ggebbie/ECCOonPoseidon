include("../../../src/intro.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, DSP, PyCall, 
    RollingFunctions, LinearAlgebra, Statistics, FFTW
import PyPlot as plt 
import NaNMath as nm
include(srcdir("config_exp.jl"))

include(srcdir("MeshArraysPlots.jl"))
cmo = pyimport("cmocean.cm");
tecco = 1992+1/24:1/12:2018
nz = 50

(ϕ,λ) = latlonC(γ)
area = readarea(γ)
λ_wrap = wrap_λ(λ)
ocean_mask = wet_pts(Γ)

region = "PAC"
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; region = region)

region = "NPAC"
NPAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; region = region)

include(srcdir("plot_and_dir_config.jl"))

normalize(x, dims) = (x .- mean(x, dims = dims)) ./ std(x, dims = dims)
normalize(x) = (x[:] .- mean(x[:])) ./ std(x[:])

corr(x, y) = cov(normalize(x), normalize(y)) / (var(normalize(x)) * var(normalize(y)))
corr(x, y, dims) = cov(normalize(x, dims), normalize(y, dims); dims = dims) / (var(normalize(x); dims = dims) * var(normalize(y); dims = dims))

vars =  ["only_kappa",  "iter0_bulkformula"]
ϕ_avg = jldopen(datadir("Ψ_Eul_timeseries_PAC_only_init" *".jld2"))["ϕ_avg"]
ϕ_ref = findall( 23 .<= ϕ_avg .<= 56)
z_ref = findall( -3300 .<= -z[:].<= -2000)

Ws = Dict(); ΨEul = Dict()
wθ_resid = Dict(); wθ_eul = Dict(); wθ_bol = Dict()

for (i, expname) in enumerate(vars)
    read_file = jldopen(datadir("Ψ_Eul_timeseries_PAC_" * expname *".jld2"))
    Ψ_exp = read_file["Ψ_exp_timeseries"]
    ΨEul[expname] = 1 .* Ψ_exp
    # Ws[expname] = -(Ψ_exp[[15, 40],  ϕ_ref[1], :] .- Ψ_exp[[15, 40],  ϕ_ref[end], :])
    # Ws[expname] = -(Ψ_exp[[15, 40],  ϕ_ref[1], :] .- Ψ_exp[[15, 40],  ϕ_ref[end], :])
    fname = datadir("native/" * expname * "NPAC_THETA_budget_Wθ_eul_bol_2to3.jld2")
    vars = jldopen(fname)["dθ"]
    wθ_resid[expname] =  vars["wθBot"] .- vars["wθTop"]
    wθ_bol[expname] = vars["wθBot_bol"] .- vars["wθTop_bol"]
    wθ_eul[expname] = wθ_resid[expname] .- wθ_bol[expname]
end

# Ws["only_kappa"] .-= Ws["iter0_bulkformula"]
wθ_resid["only_kappa"] .-= wθ_resid["iter0_bulkformula"]
ΨEul["only_kappa"] .-= ΨEul["iter0_bulkformula"]

VAR = -1e-6 .* ΨEul["only_kappa"][z_ref, ϕ_ref, :]
Y = reshape(1 .* VAR, :, 312)
Y = Y .- mean(Y, dims = 2)
Y[isnan.(Y)] .= 0.0

s = svd(Y; full = true)
PC1 = vec(s.U[:, 1]' * Y)
F = (PC1'*PC1)\PC1' # least squares estimator

objective = 1e-6 .* ΨEul["only_kappa"]
objective = objective .- mean(objective, dims = 3)
regress = reshape(objective , :, 312)* F' 
regress = reshape(regress, size(objective)[1:2]...)

fig,axs=plt.subplots(figsize = (12.5, 7))
Ψ_bounds = round(nm.maximum(abs.(regress)), digits = 2)
levels = -Ψ_bounds:0.03:Ψ_bounds
CM = axs.contourf(ϕ_avg, z,  regress, cmap=cmo.balance, levels = levels, 
vmin = -Ψ_bounds, vmax = Ψ_bounds, extend = "both")
axs.invert_yaxis()
fig.colorbar(CM, orientation = "horizontal", fraction = 0.05)
fig

fig,axs=plt.subplots(figsize = (12.5, 7))
axs.plot(tecco, PC1)
# axs.plot("First PC of ")
fig

reconstruction = zeros(size(objective))
for i in 1:size(reconstruction, 1), j in 1:size(reconstruction, 2)
    reconstruction[i, j, :] .= regress[i, j] .* PC1
end

denom = sum((objective .- mean(objective, dims = 3)).^2, dims = 3)[:, :, 1]
num = sum((reconstruction .- objective).^2, dims = 3)[:, :, 1]
Rsq = 100 .* (1 .- (num ./ denom))


fig,axs=plt.subplots(figsize = (12.5, 7))
levels = 0:10:100
CM = axs.contourf(ϕ_avg, z,  Rsq, levels = levels, cmap = cmo.thermal, 
vmin = 0, vmax = 100, extend = "both")
axs.invert_yaxis()
fig.colorbar(CM, orientation = "horizontal", fraction = 0.05)
fig