include("../../src/intro.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, DSP, PyCall, 
    RollingFunctions, LinearAlgebra, Statistics, FFTW, LaTeXStrings, Optim,
    GLM
import PyPlot as plt 
import NaNMath as nm
include(srcdir("config_exp.jl"))
@pyimport matplotlib.patches as patches

include(srcdir("MeshArraysPlots.jl"))

tecco = 1992+1/24:1/12:2018
nz = 50

(ϕ,λ) = latlonC(γ)
area = readarea(γ)
λ_wrap = wrap_λ(λ)
ocean_mask = wet_pts(Γ)
cmo = pyimport("cmocean.cm")

region = "PAC"
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; region = region)

region = "NPAC"
NPAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; region = region)

cell_depths = get_cell_thickness(PAC_msk, ΔzF, Γ.hFacC); 

include(srcdir("plot_and_dir_config.jl"))

ϕ_avg = zonal_average(ϕ, PAC_msk .* area); wherenan = (!isnan).(ϕ_avg .* 1)

normalize(x, dims) = (x .- mean(x, dims = dims)) ./ std(x, dims = dims)
normalize(x) = (x[:] .- mean(x[:])) ./ std(x[:])
square_error(x, y) = sum((x .- mean(x) .* (y .- mean(y))))
corr(x, y) = cov(x, y) / (std(x) * std(y))
var_explained(x, y) = 100 * (1 - (var(y - x) / var(y)) )
# corr(x, y, dims) = cov(x, , normalize(y, dims); dims = dims) / (var(normalize(x); dims = dims) * var(normalize(y); dims = dims))
# obj_func(t, x) =  (x[1] .* exp.(-t .* x[2])) .* (t .< x[6]) .+ (x[3] .* exp.(-t .* x[4])) .* (t .>= x[6]) .+ x[5] 
obj_func(t, x) =  (x[1] .* exp.(-t .* x[2])) .+ x[3] 
# obj_func(t, x) =  (x[1] .* exp.(-t .* x[2]) .* (1 .- exp.(-t .* x[3]))) .+ x[4] 

x_est = (4 .* collect(1:100)) .+ 1000*rand(100) .+ 1000 
corr(x_est, collect(1:100))^2
var_explained(4 .* collect(1:100) .+ 1000, x_est)

sum_squares(y, t, x) = sum((y .- obj_func(t, x)).^2) 

vars =  ["only_init",  "iter0_bulkformula"]
ϕ_avg = jldopen(datadir("Ψ_Eul_timeseries_PAC_only_init" *".jld2"))["ϕ_avg"]
Ws = Dict(); ΨEul = Dict(); ΨBol = Dict(); ΨEulBol = Dict()
# z_ref = findall( -3300 .<= -z[:].<= -2000)
wθ_resid = Dict(); wθ_eul = Dict(); wθ_bol = Dict()

for (i, expname) in enumerate(vars)
    read_file = jldopen(datadir("Ψ_Eul_timeseries_PAC_" * expname *".jld2"))
    Ψ_exp = read_file["Ψ_exp_timeseries"]
    ΨEul[expname] = 1 .* Ψ_exp

    read_file = jldopen(datadir("Ψ_EulBol_timeseries_PAC_" * expname *".jld2"))
    Ψ_exp = read_file["Ψ_exp_timeseries"]
    ΨBol[expname] = Ψ_exp .- ΨEul[expname]
    ΨEulBol[expname] = 1 .* Ψ_exp
end

ΨEulBol["only_init"] .-= ΨEulBol["iter0_bulkformula"]
ΨEul["only_init"] .-= ΨEul["iter0_bulkformula"]
ΨBol["only_init"] .-= ΨBol["iter0_bulkformula"]

Wres = -(ΨEulBol["only_init"][:, 1:end-2, :] .- ΨEulBol["only_init"][:, 3:end, :])
# WBol = -(ΨBol["only_wind"][:, 1:end-2, :] .- ΨBol["only_wind"][:, 3:end, :])

ϕ_avg = zonal_average(ϕ, PAC_msk .* area); 
area_avg = zonal_sum(PAC_msk .*area)
area_avg = area_avg[(!isnan).(ϕ_avg), :][:]
ϕ_avg = ϕ_avg[(!isnan).(ϕ_avg), :][:]
t = (collect(tecco) .- tecco[1]);

# W =  1 .* Wres[35:44, 70:end, :]
W =  1 .* Wres[:, 70:end, :]

nZ, nY, nT = size(W)  
Y = reshape(W, (nZ * nY, nT))
Y = Y .- mean(Y, dims = 2)
Y[isnan.(Y)] .= 0.0

s = svd(Y; full = true)
iPC = 5
PC1 = (s.U[:, 1:iPC]' * Y)'
var_exp_tot = 100 * round.(s.S.^2 ./ sum(s.S.^2), digits = 2)
var_exp_tot = Int.(floor.(var_exp_tot))
cumsum(var_exp_tot)[1:7]

ΨCorr = zeros(size(W)[1:2])

X = 1 .* PC1 #predictors
for i in 1:size(W, 1), j in 1:size(W, 2)
    y = Float64.(W[i, j, :][:])
    if (!isnan)(sum(y))
        # data = data .* 1e13
        ols = lm(X, y)
        ΨCorr[i, j] = round(r2(ols); digits=5)
    else
        ΨCorr[i, j] = NaN
    end
end


sns.set_theme(context = "paper", style = "ticks",
              palette = colors, rc = custom_params);

fig = plt.figure(constrained_layout=true, figsize = (10, 4))
gs = fig.add_gridspec(1, 6)
obj(i,j) = get(gs, (i,j))
slice(i,j) = pycall(pybuiltin("slice"), PyObject, i,j)
f3_ax1 = fig.add_subplot(obj(0, slice(0, 3)))
f3_ax2 = fig.add_subplot(obj(0, slice(3, 6)))

f3_ax2.set_title("First " * string(iPC) * " PCs (Total Variance of " * L"\Delta^{\mathbf{I}} W"* " Explained = " * string(cumsum(var_exp_tot)[iPC]) * "%)")
[f3_ax2.plot(tecco, PC1[:, i] ./ std(PC1[:, i]), linewidth = 1, label = "PC " * string(i) * "(" * string(var_exp_tot[i]) * "%)") for i = 1:iPC]
f3_ax2.legend(frameon=false, ncols = 3, loc = "lower center")

f3_ax1.invert_yaxis()
f3_ax1.set_title("Variance of " *L"\Delta^{\mathbf{I}} W" *" Explained by First " * string(iPC) * " PCs")
CM = f3_ax1.pcolormesh(ϕ_avg[70:end-2], z[30:45], 100 .* ΨCorr, cmap=cmo.thermal, 
vmin = 0, vmax = 100)
f3_ax1.set_facecolor("black")
f3_ax1.set_xticks(0:20:60)
lab = string.(abs.(collect(0:20:60)))
lab = lab .* ["", "°N", "°N", "°N"]
f3_ax1.set_xlim(20, 60)
f3_ax1.set_xticklabels(lab)
fig.colorbar(CM, ax = f3_ax1, orientation = "horizontal")
fig