include("../../../src/intro.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, DSP, PyCall, 
    RollingFunctions, LinearAlgebra, Statistics, FFTW, LaTeXStrings, Optim, Statistics
import PyPlot as plt 
import NaNMath as nm
using LsqFit
using ImageFiltering

include(srcdir("config_exp.jl"))
@pyimport matplotlib.patches as patches
@pyimport matplotlib.colors as mcolors

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

ϕ_avg = zonal_average(ϕ, PAC_msk .* area); 
area_avg = zonal_sum(PAC_msk .*area)
area_avg = area_avg[(!isnan).(ϕ_avg), :][:]
ϕ_avg = ϕ_avg[(!isnan).(ϕ_avg), :][:]
t = (collect(tecco) .- tecco[1]);

Wres = -(ΨEulBol["only_init"][:, 1:end-2, :] .- ΨEulBol["only_init"][:, 3:end, :])

W = 1 .* Wres[:, :, :]

@. m(t, p) = p[1] * exp.(-p[2] * t) + p[3]
function jacobian_model(x,p)
    J = Array{Float64}(undef, length(x), length(p))
    @. J[:,1] = exp(-x*p[2])     #dmodel/dp[1]
    @. @views J[:,2] = -x*p[1]*J[:,1] #dmodel/dp[2], thanks to @views we don't allocate memory for the J[:,1] slice
    @. J[:,3] = 1 #dmodel/dp[2], thanks to @views we don't allocate memory for the J[:,1] slice
    return J
end

istart = 12 * 1
y = W[3, 12, :]
p0 = [0.0, 1e-7, 0.0]

tstarts = collect(1:12*2:12*15)
var_exp = []
for istart in tstarts
    fit = curve_fit(m, jacobian_model, t[istart:end] .- t[istart], y[istart:end], p0; lower = [-Inf, 0.0, -Inf], )
    estimate3 = m(t[istart:end] .- t[istart], fit.param) 
    push!(var_exp, var_explained(estimate3, y[istart:end]))
end

fig, ax = plt.subplots()
ax.bar(tstarts ./ 12, var_exp)
fig

istart = tstarts[findmax(var_exp)[2]]
fig, ax = plt.subplots()
fit = curve_fit(m, jacobian_model, t[istart:end] .- t[istart], y[istart:end], p0; lower = [-Inf, 0.0, -Inf], )
estimate3 = m(t[istart:end] .- t[istart], fit.param) 
ax.plot(t[istart:end] .- t[istart], y[istart:end])
ax.plot(t[istart:end] .- t[istart], estimate3, label = "estimate 2")
ax.legend()
fig

estimates = []
tstarts = collect(1:12:12*25)
for istart in tstarts
    push!(estimates, zeros(size(W)[1:2]..., length(tecco) - istart + 1))
end

total_variance_explained = zeros(length(tstarts))
wherenans = isnan.(sum(W, dims = 3)[:, :, 1])
W_var = 1 .* W; W_var[isnan.(W_var)] .= 0.0

COEFS = zeros(size(W)[1:2]..., length(tstarts), 3)
converged_ = zeros(size(W)[1:2]..., length(tstarts))

for k in 1:length(tstarts)
    print(k)
    istart = tstarts[k]

    for i in 1:size(W)[1], j in 1:size(W)[2]
        y = 1 .* W[i, j, :][:]
        if (!wherenans[i, j])
            fit = curve_fit(m, jacobian_model, t[istart:end] .- t[istart], y[istart:end], p0; lower = [-Inf, 0.0, -Inf], )
            estimates[k][i, j, :] = m(t[istart:end] .- t[istart], fit.param) 
            COEFS[i, j, k, :] .= fit.param
            converged_[i, j, k] = fit.converged

        else 
            estimates[k][i, j, :] .= NaN
            COEFS[i, j, k, :] .= NaN
            converged_[i, j, k] = NaN
        end
    end
    estimates_var = 1 .* estimates[k]; estimates_var[isnan.(estimates_var)] .= 0.0
    total_variance_explained[k] = 1 .- (var(W_var[:, :, istart:end] .- estimates_var) ./ var(W_var[:, :, istart:end]))
end

var_explained_byt = 0 .* W[:, :, 1:length(tstarts)]
where_variance_explained =  0 .* W[:, :, 1:length(tstarts)]

for k in 1:length(tstarts)
    istart = tstarts[k]
    var_explained_byt[:, :, k] .= 1 .- (var(W[:, :, istart:end] .- estimates[k], dims = 3)[:, :, 1] ./ var(W[:, :, istart:end], dims = 3)[:, :, 1])
    where_variance_explained[:, :, k] .= tstarts[k]
end

where_max = findmax(var_explained_byt, dims = 3)[2]
converged_mask = iszero.(Float64.(converged_[where_max])[:, :, 1])
var_max = 100 .*  var_explained_byt[where_max][:, :, 1]

wherenan = isnan.(var_max)[:, :, 1]
where_max = findmax(var_explained_byt, dims = 3)[2]
time_max = (where_variance_explained[where_max] ./ 12) .+ 1992
time_max = time_max[:, :, 1]

decay_times = inv.(COEFS[:, :, :, 2]  ./ 12)[where_max][:, :, 1]
equil_val = (86400 * 365 * COEFS[:, :, :, 3][where_max][:, :, 1]) ./ area_avg[2:end-1]'
time_max[isnan.(equil_val)] .= NaN

fig,axs=plt.subplots(2, 2, figsize = (17, 17))

axs[1].set_title("Model Start Year, " * L"t_0")
CM = axs[1].contourf(ϕ_avg[2:end-1], z, time_max, cmap="Spectral_r", levels = 1992:1:2017, vmin = 1992, vmax = 2017, extend = "max")
fig.colorbar(CM, ax = axs[1], ticks = 1992:5:2017,
orientation = "horizontal", fraction  =0.04, label = "year")

axs[2].set_title("Equillibrium Rate, " * L"\lambda")
CM = axs[2].contourf(ϕ_avg[2:end-1], z, decay_times, cmap="Spectral_r", levels = 0:1:100, vmin = 0, vmax = 100, extend = "max")
fig.colorbar(CM, ax = axs[2], orientation = "horizontal", fraction  =0.04, 
ticks = 0:20:100, label = "years")

axs[3].set_title("Equillibrium Value " * L"\overline{W_\infty^{res}}")
CM = axs[3].contourf(ϕ_avg[2:end-1], z, equil_val, cmap=cmo.balance, levels = -9:1.5:9, vmin = -9, vmax = 9, extend = "both")
fig.colorbar(CM, ax = axs[3], orientation = "horizontal", fraction  =0.04,
ticks = -9:3:9, label = "m year" * L" ^{-1}")

axs[4].set_title("Variance explained by Equillibrium Model")
CM = axs[4].contourf(ϕ_avg[2:end-1], z, var_max, cmap=cmo.thermal, levels = 0:10:100, 
vmin = 0, vmax = 100, shading = "gourad")
fig.colorbar(CM, ax = axs[4], orientation = "horizontal", fraction  =0.04, label = "% Variance Explained")

# cmap_mask = mcolors.ListedColormap(["white", "white", "white"])
mask = 1.0 .* converged_mask; mask[iszero.(mask)] .= NaN
none_map = mcolors.ListedColormap(["none"])

for ax in axs
    ax.set_facecolor("black")
    ax.invert_yaxis()
    ax.set_xticks(-40:20:65)
    ax.set_xlim(-34, 65)
    lab = string.(abs.(collect(-40:20:60)))
    lab = lab .* ["°S", "°S", "", "°N", "°N", "°N"]
    ax.set_xticklabels(lab)
    ax.set_ylabel("Depth [m]", fontweight = "bold")
    ax.set_xlabel("Latitude", fontweight = "bold")
    rect = patches.Rectangle((23, 2000), 59, 1000, linewidth=3, edgecolor="black",facecolor="none", alpha = 0.7)
    ax.add_patch(rect)
    # cs = ax.contourf(ϕ_avg[2:end-1], z,mask , levels = 3, cmap = cmap_mask)
    # cs = ax.contourf(ϕ_avg[2:end-1], z,converged_mask, levels = 3, hatches=["", ".."],  alpha=0.0)
        cs = ax.pcolor(ϕ_avg[2:end-1], z,mask, cmap=none_map,
                   hatch="\\\\", edgecolor="black", lw=0, zorder=2)
end
cs = axs[3].pcolor(ϕ_avg[2:end-1], z,mask, cmap=none_map,
hatch="\\ \\", edgecolor="red", lw=0, zorder=2)
fig.tight_layout()

fig_labs = uppercase.(["a", "b", "c", "d", "e", "f"])
for (i, a) in enumerate(axs)
    a.annotate(fig_labs[i], (0.93, 0.02), fontsize = 30, color = "white", 
    xycoords="axes fraction", fontweight = "bold")
end
fig
fig.savefig(plotsdir("native/paper_figures/6.equillibrium_model.png"), dpi = 400)

close("all")