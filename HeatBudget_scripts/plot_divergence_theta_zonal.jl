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
function extrema(d::Dict)
    extremas = map(x-> nm.extrema(x), values(d))
    clims = (minimum(minimum.(extremas)), maximum(maximum.(extremas)))
    return clims
end
function plot_zonal_contours(X, Y, zonal_var, clims, title)
    jcf = Plots.contourf(X, Y, reverse(zonal_var, dims = 1),
    xlabel = "latitude [º]",
    xticks = -70:20:70,
    ylabel = "depth [m]", 
    linewidth = 0.1,
    levels = 16,
    clim = clims,
    c = :balance, 
    colorbar_title= "°C",
    title = title, 
    thickness_scaling = 1.5)
    return jcf
end

function plot_all_zonal_θ(X, Y, θ_zonal, nt::Int64, tecco, clims; anom = "")
    exps = collect(keys(θ_zonal))

    a = Animation()

    ps = Vector{Any}(missing, length(exps))
    for tt in 1:nt
        for (i, ex) in enumerate(exps)
            ps[i] = plot_zonal_contours(X, Y, θ_zonal[ex][:, :, tt], clims, ex)
        end
        l = @layout [a;b;c]
        p = Plots.plot(ps[1], ps[2], ps[3], size=(1401,1401), link = :both,
        plot_title =  "Zonally averaged θ "*anom* region * "\n Time: " * 
        string(round(tecco[tt], digits = 6)),
        plot_titlefontsize	= 17, plot_titlevspan = 0.075, layout = l)
        frame(a, p)
    end
    return a
end


function make_zonal_pics(X, Y, θ_zonal, tt::Int64, tecco)
    exps = collect(keys(θ_zonal))
    extremas = map(x-> nm.extrema(x), values(θ_zonal))
    clims = (minimum(minimum.(extremas)), maximum(maximum.(extremas)))
    clims = (-maximum(abs.(clims)), maximum(abs.(clims)))

    ps = Vector{Any}(missing, length(exps))
    for (i, ex) in enumerate(exps)
        baseline = θ_zonal["iter0_bulkformula"][:, :, tt]
        cf = plot_zonal_contours(X, Y, θ_zonal[ex][:, :, tt], clims, ex)
        ps[i] = cf
    end
    l = @layout [a;b;c]
    p = Plots.plot(ps[1], ps[2], ps[3], size=(1401,1401), link = :both,
    plot_title =  "Zonally averaged θ "* region * "\n Time: " * 
    string(round(tecco[nt], digits = 6)),
    plot_titlefontsize	= 17, plot_titlevspan = 0.075, layout = l)
    return p 
end

function plot_zonal_θ(X, Y, θ_zonal, exp, nt::Int64, tecco, clims; anom = "")
    a = Animation()

    for tt in 1:nt
        ps = plot_zonal_contours(X, Y, θ_zonal[exp][:, :, tt], clims, exp)
        p = Plots.plot(ps, size=(1400, 700), link = :both,
        plot_title =  "Zonally averaged θ " * anom * region * "\n Time: " * 
        string(round(tecco[tt], digits = 6)) * " \n ",
        plot_titlefontsize = 12, plot_titlevspan = 0.1)
        frame(a, p)
    end
    return a
end

Y = reverse(z[lvls])
X = OHC_helper.ma_zonal_sum(ϕ .* area) ./ OHC_helper.ma_zonal_sum(area)
nt = length(tecco)
exps = collect(keys(shortnames))
clims = extrema(θ_zonal)
a = plot_all_zonal_θ(X, Y, θ_zonal, nt, tecco, clims)
mp4(a, plotsdir() * "/OHC_Divergence/Gifs/" * "θ_zonalavg_" * 
region * "_" * suffix * ".mp4", fps = 5)

exp = "iter129_bulkformula"
expname = "iter129"
a = plot_zonal_θ(X, Y, θ_zonal, exp, nt, tecco, clims)
mp4(a, plotsdir() * "/OHC_Divergence/Gifs/" * "θ_zonalavg_" *expname* 
        region * "_" * suffix * ".mp4", fps = 5)

#taking the anomaly with iter0 as the reference
θ_zonal_anom = Dict(key => θ_zonal[key] .- θ_zonal["iter0_bulkformula"] 
                    for key in keys(θ_zonal))
clims = extrema(θ_zonal_anom)
clims = (-maximum(abs.(clims)), maximum(abs.(clims)))
a = plot_all_zonal_θ(X, Y, θ_zonal_anom, nt, tecco, clims)
mp4(a, plotsdir() * "/OHC_Divergence/Gifs/" * "θ_zonalavg_iter0anom" * 
region * "_" * suffix * ".mp4", fps = 5)

exp = "iter129_bulkformula"
expname = "iter129"
a = plot_zonal_θ(X, Y, θ_zonal_anom, exp, nt, tecco, clims)
mp4(a, plotsdir() * "/OHC_Divergence/Gifs/" * "θ_zonalavg_iter0anom_" *expname* 
        region * "_" * suffix * ".mp4", fps = 5)

#taking the anomaly with the first time step as the reference 
θ_zonal_anom = Dict(key => θ_zonal[key] .- θ_zonal["iter129_bulkformula"][:, :, 1] 
                    for key in keys(θ_zonal))
exp = "iter129_bulkformula"
expname = "iter129"
clims = nm.extrema(θ_zonal_anom[exp])
clims = (-maximum(abs.(clims)), maximum(abs.(clims)))

a = plot_zonal_θ(X, Y, θ_zonal_anom, exp, nt, tecco, clims)
mp4(a, plotsdir() * "/OHC_Divergence/Gifs/" * "θ_zonalavg_initanom_" *expname* 
        region * "_" * suffix * ".mp4", fps = 5)
# gif(a, fps = 10)

# preview if needed gif(anim, fps = 5)


