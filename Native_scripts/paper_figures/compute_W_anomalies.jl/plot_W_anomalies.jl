include("../../../src/intro.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, LaTeXStrings,
    PyCall, BenchmarkTools
using NaNStatistics
import PyPlot as plt

include(srcdir("plot_and_dir_config.jl"))

adjust_exps = jldopen(datadir(region * "W_interpolated_lev40.jld2"))["adjust_exps"]

lat, lon = adjust_exps["lat"], adjust_exps["lon"]
mask = adjust_exps["RegularMask"]
weights = cos.(deg2rad.(lat)) .* mask
weights= cos.(deg2rad.(abs.(lat))).*(111.1*111.1*1000*1000) .* mask
weighted_nanmean(x, weights; dims = 1) =  NaNStatistics.nansum(x .* weights, dims = dims) ./ NaNStatistics.nansum(weights, dims = dims)
weighted_nansum(x, weights; dims = 1) =  NaNStatistics.nansum(x .* weights, dims = dims) 
NaNStatistics.nansum
W129 = adjust_exps["mean_tau_adjusts"]
W0 = adjust_exps["mean_tau_noadjusts"]

ΔW129 = weighted_nansum(W129 .- 0, weights, dims = 2)[:, 1, :]

fig, ax = plt.subplots(1, 1, sharey = true)
lons = 1 .* lon[:, 1]
lons[lons .< 0] .+= 360
ax.pcolormesh(lons, 1:312, ΔW129', vmin = -1e5, vmax = 1e5, cmap = "bwr")
ax.set_ylim([0, 50])
fig

fig, ax = plt.subplots(1, 2, sharey = true)
lons = lon[:, 1]
lons[lons .< 0] .+= 360
ax[1].plot(1:312, NaNStatistics.nansum(ΔW129, dims = 1)[:])
fig

ΔW129 = weighted_nansum(W129 , weights, dims = 2)[:, 1, :]
ΔWSC = weighted_nansum(WSC, weights , dims = 2)[:, 1, :]
ΔW0 = weighted_nansum(W0 , weights, dims = 2)[:, 1, :]


fig, ax = plt.subplots(1, 3, sharey = true)
lons = lon[:, 1]
lons[lons .< 0] .+= 360
ax[1].plot(1:180, NaNStatistics.nansum(ΔW129, dims = 1)[:])
ax[2].plot(1:180, NaNStatistics.nansum(ΔWSC, dims = 1)[:])
ax[3].plot(1:180, NaNStatistics.nansum(ΔW0, dims = 1)[:])

fig

function mean_of_every_n_elements(arr, n)
    means = []
    for i in 1:n:length(arr)
        if i + n - 1 <= length(arr)
            push!(means, mean(arr[i:i+n-1]))
        else
            push!(means, mean(arr[i:end]))  # for remaining elements if they are less than n
        end
    end
    return means
end

fig, ax = plt.subplots(1, 3, sharey = true)
lons = lon[:, 1]
lons[lons .< 0] .+= 360
ax[1].plot(mean_of_every_n_elements(1:180, 6), mean_of_every_n_elements(NaNStatistics.nansum(ΔW129, dims = 1)[:], 6))
ax[2].plot(mean_of_every_n_elements(1:180, 6),  mean_of_every_n_elements(NaNStatistics.nansum(ΔWSC, dims = 1)[:], 6))
ax[3].plot(mean_of_every_n_elements(1:180, 6),  mean_of_every_n_elements(NaNStatistics.nansum(ΔW0, dims = 1)[:], 6))

fig
