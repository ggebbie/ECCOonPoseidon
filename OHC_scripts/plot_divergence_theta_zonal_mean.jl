using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, JLD2,
GoogleDrive,NCDatasets, NetCDF, Printf, RollingFunctions,
LaTeXStrings
using PyCall
using DataFrames, GLM
using Plots
import NaNMath as nm

using ColorSchemes
include(srcdir("config_exp.jl"))

import Base: extrema

function plot_zonal_contours(X, Y, zonal_var, clims, title)
    jcf = Plots.contourf(X, Y, zonal_var,
    xlabel = "latitude [º]",
    xticks = -70:20:70,
    ylabel = "depth [m]", 
    linewidth = 0.1,
    levels = 25,
    clim = clims,
    c = :delta, 
    colorbar_title= "Sv",
    title = title, 
    thickness_scaling = 1.5)
    return jcf
end

function plot_θ_mean(X, Y, θ_zonal, θ_labels, clims, θ_sym)
    ps = Vector{Any}(missing, length(keys(θ_zonal)))
    i = [0] 
    for ex in keys(θ_zonal)
        i .+= 1
        ps[i[1]] = plot_zonal_contours(X, Y, reverse(θ_zonal[ex], dims = 1),
        clims, θ_labels[i[1]])
    end
    l = @layout [a;b;c;d]
    p = Plots.plot(ps[1], ps[2], ps[3],ps[4],size=(1401,2101), link = :both,
    plot_title =  "Time mean "* θ_sym  * " "* region,
    plot_titlefontsize	= 17, plot_titlevspan = 0.025, layout = l)
    return p 
end

Y = reverse(z[lvls])
X = OHC_helper.ma_zonal_sum(ϕ .* area) ./ OHC_helper.ma_zonal_sum(area)
θ_zonal_mean = Dict(key => dropdims(mean(θ_zonal[key], dims = 3), dims = 3) for key in keys(θ_zonal))
θ_zonal_mean_anom = Dict(key => value .- θ_zonal_mean["iter0_bulkformula"] 
                        for (key, value) in θ_zonal_mean)
            
θ_labels = [L"\theta^{129}", L"\theta^{\Delta F}",L"\theta^{\Delta T}", L"\theta^{0}"]
θ_sym = "θ"
clims =  extrema(θ_zonal_mean);

p = plot_θ_mean(X, Y, θ_zonal_mean, θ_labels, clims, θ_sym)
savefig(p,plotsdir() * "/OHC_Divergence/Contours/" * "θmean_"* 
region * ".png")

clims = extrema(θ_zonal_mean_anom); clims = (-maximum(abs.(clims)), maximum(abs.(clims)))

p = plot_θ_mean(X, Y, θ_zonal_mean_anom, θ_labels, clims, θ_sym)
savefig(p,plotsdir() * "/OHC_Divergence/Contours/" * "θmean_iter0_anom_"* 
region * ".png")