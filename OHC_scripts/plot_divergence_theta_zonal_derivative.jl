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


function plot_zonal_contours(X, Y, zonal_var, clims, title, color = :amp)
    jcf = Plots.contourf(X, Y, 1e2 .* zonal_var,
    xlabel = "latitude [º]",
    xticks = -70:20:70,
    ylabel = "depth [m]", 
    linewidth = 0.5,
    levels = 10,
    clim = clims,
    c = color, 
    colorbar_title= "cK / m",
    title = title, 
    thickness_scaling = 1.5)
    return jcf
end

function plot_dθ_mean(X, Y, θ_zonal, θ_labels, clims, θ_sym, color = :amp)
    ps = Vector{Any}(missing, length(keys(θ_zonal)))
    i = [0] 
    for ex in keys(θ_zonal)
        i .+= 1
        ps[i[1]] = plot_zonal_contours(X, Y, reverse(θ_zonal[ex], dims = 1),
        clims, θ_labels[i[1]],color)
    end
    l = @layout [a;b;c;d]
    p = Plots.plot(ps[1], ps[2], ps[3],ps[4],size=(1401,2101), link = :both,
    plot_title =  "Time mean "* θ_sym  * " "* region,
    plot_titlefontsize	= 17, plot_titlevspan = 0.025, layout = l)
    return p 
end

function expdθdz(θ_zonal, z)
    dθdz_zonal = Dict()
    dz = z[1:end-1] .- z[2:end]
    for (key, value) in θ_zonal
        dθdz_zonal[key] = (value[1:end-1, :] .- value[2:end, :]) ./ dz
    end
    dθdz_zonal
end
Y = reverse((z[lvls][1:end-1] .+ z[lvls][2:end])/2)
X = OHC_helper.ma_zonal_sum(ϕ .* area) ./ OHC_helper.ma_zonal_sum(area)

dθdz_zonal_mean = expdθdz(θ_zonal_mean, z[lvls])
dθdz_zonal_mean_anom = Dict(key => value .- dθdz_zonal_mean["iter0_bulkformula"] 
                        for (key, value) in dθdz_zonal_mean)
            
dθdz_labels = "∂" .* θ_labels .* "/∂z"
dθdz_sym = "∂θ/∂z"
clims =  1e2 .* extrema(dθdz_zonal_mean); #cK

p = plot_dθ_mean(X, Y, dθdz_zonal_mean, dθdz_labels, clims, dθdz_sym)
savefig(p,plotsdir() * "/OHC_Divergence/Contours/" * "dθdzmean_"* 
region * "_" * suffix * ".png")

clims = extrema(dθdz_zonal_mean_anom); clims = 1e2 .* (-maximum(abs.(clims)), maximum(abs.(clims)))#cK
θ_sym = "∂θ/∂z Anomaly"
dθdz_labels_anom =  [sym * L" - \theta^{0}" for sym in θ_labels]
dθdz_labels_anom = "∂(" .*  dθdz_labels_anom .*  ")/∂z"
p = plot_dθ_mean(X, Y, dθdz_zonal_mean_anom, dθdz_labels_anom, clims, θ_sym, :delta)
savefig(p,plotsdir() * "/OHC_Divergence/Contours/" * "dθdzmean_iter0_anom_"* 
region * "_" * suffix * ".png")