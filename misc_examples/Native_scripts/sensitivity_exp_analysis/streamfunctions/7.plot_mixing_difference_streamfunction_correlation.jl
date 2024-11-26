include("../../../src/intro.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, DSP, PyCall, 
    RollingFunctions, LinearAlgebra, Statistics, FFTW, LaTeXStrings
import PyPlot as plt 
import NaNMath as nm
include(srcdir("config_exp.jl"))

include(srcdir("MeshArraysPlots.jl"))

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

ϕ_avg = zonal_average(ϕ, PAC_msk .* area); wherenan = (!isnan).(ϕ_avg .* 1)

normalize(x, dims) = (x .- mean(x, dims = dims)) ./ std(x, dims = dims)
normalize(x) = (x[:] .- mean(x[:])) ./ std(x[:])

corr(x, y) = cov(normalize(x), normalize(y)) / (var(normalize(x)) * var(normalize(y)))
corr(x, y, dims) = cov(normalize(x, dims), normalize(y, dims); dims = dims) / (var(normalize(x); dims = dims) * var(normalize(y); dims = dims))

function regress_matrix(x)

    nt = length(x)

    # would be neat to return whatever type goes in
    E = zeros(Float32,nt,1)
    E[:, 1] .= x
    
    F = (E'*E)\E' # least squares estimator

    return F
end

pad_arr(x) = vcat(reverse(x), x, reverse(x))
unpad_arr(x, arr_len) = x[arr_len+1:2*arr_len]
function low_pass(signal)
    nt = length(signal)
    signal_mean = mean(signal)
    ff = digitalfilter(Lowpass(1/24, fs = 1),Butterworth(2))
    filtered_signal = filtfilt(ff, pad_arr(signal .- signal_mean))
    return unpad_arr(filtered_signal, nt) .+ signal_mean
end

function get_transports(diagpath::Dict{String, String}, 
    expname::String, γ::gcmgrid)

    filelist = searchdir(diagpath[expname],"state_2d_set1") # first filter for state_3d_set1
    datafilelist_τ  = filter(x -> occursin("data",x),filelist) # second filter for "data"

    nt = length(datafilelist_τ);
    nϕ = length(ϕ_avg)
    curlτ = zeros(Float32, nϕ, nt);

    @time for tt = 1:nt
        println(tt)

        Tname = datafilelist_τ[tt]
        τx, τy = extract_ocnTAU(diagpath, expname , Tname, γ)
        τcurl = MeshArrays.curl(τx, τy, Γ)
        curlτ[:, tt] .= zonal_average(τcurl, PAC_msk .* area)[:]

    end

    return curlτ
end

τ_mean = Dict()
τ_mean["only_kappa"]   = get_transports(diagpath, "only_wind", γ)

vars =  ["only_kappa",  "iter0_bulkformula"]
ϕ_avg = jldopen(datadir("Ψ_Eul_timeseries_PAC_only_init" *".jld2"))["ϕ_avg"]
Ws = Dict(); ΨEul = Dict()
z_ref = findall( -3300 .<= -z[:].<= -2000)
wθ_resid = Dict(); wθ_eul = Dict(); wθ_bol = Dict()

for (i, expname) in enumerate(vars)
    read_file = jldopen(datadir("Ψ_Eul_timeseries_PAC_" * expname *".jld2"))
    Ψ_exp = read_file["Ψ_exp_timeseries"]
    ΨEul[expname] = 1 .* Ψ_exp
end

ΨEul["only_kappa"] .-= ΨEul["iter0_bulkformula"]

ϕ_avg = zonal_average(ϕ, PAC_msk .* area); 
τ = 1 .* τ_mean["only_kappa"][(!isnan).(ϕ_avg), :]
τ = τ[2:end-1, :]
predi = -(ΨEul["only_kappa"][:, 1:end-2, :] .- ΨEul["only_kappa"][:, 3:end, :])

ΨCorr = zeros(size(predi)[1:2])
for i in 1:size(predi, 1), j in 1:size(predi, 2)
    x = predi[i, j, :][:]
    y = τ[j, :]
    ΨCorr[i, j] = corr(x, y) 
end


ϕ_avg = zonal_average(ϕ, PAC_msk .* area); 
ϕ_avg = ϕ_avg[(!isnan).(ϕ_avg), :][:]

fig = plt.figure()
ax1 = fig.add_subplot(5, 1, 1)
ax2 = fig.add_subplot(5, 1, (2, 5), sharex=ax1)

mean_Ψ =  1e-6 .* mean(predi, dims = 3)[:,:,1]
Ψ_bounds = 0.5
levels = -0.5:0.1:0.5
CM = ax2.contourf(ϕ_avg[2:end-1], z, mean_Ψ,levels = levels, cmap=cmo.balance, 
vmin = -Ψ_bounds, vmax = Ψ_bounds, extend = "both")
ax2.invert_yaxis()
ax1.plot(ϕ_avg[2:end-1, :], mean(τ, dims = 2));
fig.colorbar(CM, orientation = "horizontal")
fig


fig,axs=plt.subplots(figsize = (15, 7.5), sharex = true)
levels = -1:0.1:1
CM = axs.contourf(ϕ_avg[2:end-1], z, ΨCorr,levels = levels, cmap=cmo.balance, 
vmin = -1, vmax = 1, extend = "both")
axs.invert_yaxis()
fig.colorbar(CM)
fig
