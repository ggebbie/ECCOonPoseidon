include("../../../src/intro.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, DSP, PyCall, 
    RollingFunctions, LinearAlgebra, Statistics, FFTW, LaTeXStrings, Optim, Statistics
import PyPlot as plt 
import NaNMath as nm
using LsqFit

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
WEul = -(ΨEul["only_init"][:, 1:end-2, :] .- ΨEul["only_init"][:, 3:end, :])
# WBol = -(ΨBol["only_wind"][:, 1:end-2, :] .- ΨBol["only_wind"][:, 3:end, :])

ϕ_avg = zonal_average(ϕ, PAC_msk .* area); 
area_avg = zonal_sum(PAC_msk .*area)
area_avg = area_avg[(!isnan).(ϕ_avg), :][:]
ϕ_avg = ϕ_avg[(!isnan).(ϕ_avg), :][:]
t = (collect(tecco) .- tecco[1]);

W = 1 .* Wres[30:43, 72:end, :]

# function optimize_func(y, t)
#     sol = Optim.optimize(x -> sum_squares(y, t, x), [0.0, 0.0, 0.0], Newton());
#     coefs = Optim.minimizer(sol);
#     return coefs 
# end
# coefs = optimize_func(y[istart:end], t[istart:end] .- t[istart])
# estimate = obj_func(t[istart:end] .- t[istart], coefs) 
# using BenchmarkTools

@. m(t, p) = p[1] * exp.(-p[2] * t) + p[3]
function jacobian_model(x,p)
    J = Array{Float64}(undef, length(x), length(p))
    @. J[:,1] = exp(-x*p[2])     #dmodel/dp[1]
    @. @views J[:,2] = -x*p[1]*J[:,1] #dmodel/dp[2], thanks to @views we don't allocate memory for the J[:,1] slice
    @. J[:,3] = 1 #dmodel/dp[2], thanks to @views we don't allocate memory for the J[:,1] slice
    return J
end

istart = 12 * 1
y = W[10, 12, :]
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
tstarts = collect(1:12:12*24)
for istart in tstarts
    push!(estimates, zeros(size(W)[1:2]..., length(tecco) - istart + 1))
end

total_variance_explained = zeros(length(tstarts))
wherenans = isnan.(sum(W, dims = 3)[:, :, 1])
W_var = 1 .* W; W_var[isnan.(W_var)] .= 0.0

for k in 1:length(tstarts)
    print(k)
    istart = tstarts[k]
    COEFS = zeros(size(W)[1:2]..., length(tstarts))

    for i in 1:size(W)[1], j in 1:size(W)[2]
        y = 1 .* W[i, j, :][:]
        if (!wherenans[i, j])
            fit = curve_fit(m, jacobian_model, t[istart:end] .- t[istart], y[istart:end], p0; lower = [-Inf, 0.0, -Inf], )
            estimates[k][i, j, :] = m(t[istart:end] .- t[istart], fit.param) 
            COEFS[i, j, k] = inv(fit.param[2])
        else 
            estimates[k][i, j, :] .= NaN
        end
    end
    estimates_var = 1 .* estimates[k]; estimates_var[isnan.(estimates_var)] .= 0.0
    total_variance_explained[k] = 1 .- (var(W_var[:, :, istart:end] .- estimates_var) ./ var(W_var[:, :, istart:end]))

end


mean_variance_explained = zeros(length(tstarts))

for k in 1:length(tstarts)
    istart = tstarts[k]
    var_exp = 1 .- (var(W[:, :, istart:end] .- estimates[k], dims = 3)[:, :, 1] ./ var(W[:, :, istart:end], dims = 3)[:, :, 1])
    mean_variance_explained[k] = nm.mean(var_exp)
end

fig, ax = plt.subplots(1, 2, sharey = true)
ax[1].plot(tstarts ./ 12, 100 .* total_variance_explained)
ax[2].plot(tstarts ./ 12, 100 .* mean_variance_explained)
fig

istart = tstarts[findmax(mean_variance_explained)[2]]
k = findmax(mean_variance_explained)[2]
var_exp = 1 .- (var(W[:, :, istart:end] .- estimates[k], dims = 3)[:, :, 1] ./ var(W[:, :, istart:end], dims = 3)[:, :, 1])

fig,axs=plt.subplots( figsize = (17, 6), sharey = true)
cms = []
axs.set_facecolor("black")
CM = axs.contourf(ϕ_avg[67:end-1], z, variances[:, :, k], cmap=cmo.thermal, levels = 10, 
vmin = 0, vmax = 100, shading = "gourad")
fig.colorbar(CM, ax = axs, orientation = "horizontal", fraction  =0.04, label = "correlation")
axs.invert_yaxis()
fig