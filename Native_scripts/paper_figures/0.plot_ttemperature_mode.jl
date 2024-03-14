include("../../src/intro.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, LaTeXStrings,
    PyCall, BenchmarkTools
import PyPlot as plt

using LinearAlgebra, PyPlot, Statistics

using SavitzkyGolay

include(srcdir("plot_and_dir_config.jl"))
sns.set_theme(context = "paper", style = "ticks",
              palette = colors, rc = custom_params);
@pyimport matplotlib.patches as patches

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

ocean_mask = wet_pts(Γ)
region = "NPAC"; 
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; region)

cell_depths = get_cell_thickness(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = get_cell_volumes(area, cell_depths)
ΔV = lateral_sum(cell_volumes)

lvls = findall( -3000 .<= -z[:].<= -2000)
mid_depths(x) = vec(sum(Float32, x[lvls, :] .* ΔV[lvls], dims = 1) / sum(Float32, ΔV[lvls]))

tecco = 1992+1/24:1/12:2018; nz = 50
get_datafiles(expname, key) = filter(x -> occursin("data",x),searchdir(diagpath[expname],key) )

adjust_exps =  jldopen(datadir(region * "_temperature_sens_exps.jld2"))["adjust_exps"]
adjust_exps["Difference"] = adjust_exps["iter129_bulkformula"] .- adjust_exps["iter0_bulkformula"]
lw = 2

only_init = mid_depths(adjust_exps["only_init"])
iter0 = mid_depths(adjust_exps["iter0_bulkformula"])


offset_difference(x) = (x[2:end] .- x[1:end-1])
offset_addition(x) = (x[2:end] .+ x[1:end-1])
derivative(x, t) = offset_difference(x) ./ offset_difference(t)
interp_2_center(x) =  offset_addition(x) ./ 2

#setup d²x/dt, dx/dt, x?
signal = Float32.(only_init .- iter0)
initial_sig = 1 .* signal[1]
# signal .-=initial_sig
# signal = signal
# signal = savitzky_golay(signal, 13, 7).y
t = collect(0:(length(signal) - 1)) ./ 12

fig, ax = plt.subplots(figsize = (6.0, 5))
ax.plot(signal , label = "true diff")
fig
dsignal = derivative(signal, t)
t_interp = interp_2_center(t)
# dt_interp = offset_difference(t_interp)
ddsignal = derivative(dsignal, t_interp)
dsignal = interp_2_center(dsignal)

signal_interp = interp_2_center(signal)
signal_interp = interp_2_center(signal_interp)

b = -savitzky_golay(ddsignal, 11, 1).y
ax.plot(t, signal, label = "true diff")

A = ones(length(dsignal), 2)
A[:,1] .= savitzky_golay(dsignal, 11, 1).y
A[:,2] .= savitzky_golay(signal_interp, 11, 1).y

# A[:,2] .= 1 .*signal_interp
U1, U2 = pinv(A) * b

tmp = sqrt((U1^2) - (4*U2))
a, b = 0.5 .* [-U1 + tmp, -U1 - tmp]

X = [exp.(a.*t) .- exp.(b.*t) exp.(a.*t)] 

# C = U3 / (-a*b)
B = signal

As = pinv(X) * B
As_iter129 = 1 .* As

reconstruct = X * As
# reconstruct .+= 0.001 .* exp.(a.*t)

t_extrap = collect(0:2500) ./ 12

X = [exp.(a.*t_extrap) .- exp.(b.*t_extrap) exp.(a.*t_extrap)] 
reconstruct = X * As

fig, ax = plt.subplots( figsize = (4, 3))
ax.plot(t_extrap, reconstruct, c = "k", linestyle = "--", label = "Bi-Exponential  Model"); 
ax.plot(t, signal, label = "True Model Difference"); 
ax.set_xlabel("years from initialization")
ax.legend()
ax.set_title(L"\theta'")
ax.set_ylabel("[°C]")
fig.tight_layout()
fig

ones_ = ones(length(t))

X = [exp.(a.*t) .- exp.(b.*t) exp.(a.*t) ones_] 
B = 1 * only_init

As = pinv(X) * B
As_it129 = 1 .* As
ones_ = ones(length(t_extrap))
X = [exp.(a.*t_extrap) .- exp.(b.*t_extrap) exp.(a.*t_extrap) ones_] 
reconstruct = X * As
only_init_ode = 1 .* reconstruct
# X = [exp.(a.*t) .- exp.(b.*t) exp.(a.*t)] 

# reconstruct = X * As
# reconstruct_init = 1 .* reconstruct

fig, ax = plt.subplots(figsize = (6.0, 5))
ax.plot(only_init, label = "true diff")
# F = mean(only_init .- reconstruct)
ax.plot(reconstruct, label = "approx diff")
ax.legend()
fig


# X = [exp.(a.*t) .- exp.(b.*t) exp.(a.*t)] 
# C = U3 / (-a*b)
B = iter0 

# As = pinv(X) * B
As_iter0 = As[1:2] .- As_iter129
As_iter0 = vcat(As_iter0, As[3])
ones_ = ones(length(t_extrap))
X = [exp.(a.*t_extrap) .- exp.(b.*t_extrap) exp.(a.*t_extrap) ones_] 
reconstruct = X * As_iter0
reconstruct_i0 = 1 .* reconstruct
fig, ax = plt.subplots(figsize = (6.0, 5))
ax.plot(iter0, label = "true diff")
ax.plot(reconstruct, label = "true diff")
fig

X = [exp.(a.*t_extrap) .- exp.(b.*t_extrap) exp.(a.*t_extrap) ones_] 
reconstruct = X * As

fig, ax = plt.subplots(figsize = (6.0, 5))
ax.plot(t, only_init, label = "true diff")
ax.plot(t_extrap, only_init_ode, label = "true diff")
ax.plot(t, iter0, label = "true diff")
ax.plot(t_extrap, reconstruct_i0, label = "true diff")
fig


As_iter0 .-As

As_iter129