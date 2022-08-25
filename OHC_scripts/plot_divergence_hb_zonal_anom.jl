using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, JLD2,
GoogleDrive,NCDatasets, NetCDF, Printf, RollingFunctions,
LaTeXStrings
using PyCall
using DataFrames, GLM
using Plots
# import Plots: Animation, frame, gif, cgrad
# import Plots: contourf as jcontourf
# import Plots: plot as jplot
using ColorSchemes
include(srcdir("config_exp.jl"))


Y = reverse(z[lvls])
X = OHC_helper.ma_zonal_sum(ϕ .* area) ./ OHC_helper.ma_zonal_sum(area)
ϕ_d = similar(crop_vols)
exps = collect(keys(shortnames))
c1 = minimum([nm.minimum(θ_zonal[ex] .- θ_zonal["iter0_bulkformula"]) for ex ∈ exps])
c2 = maximum([nm.maximum(θ_zonal[ex] .- θ_zonal["iter0_bulkformula"]) for ex ∈ exps])
c1=-maximum(abs.([c1, c2]))
c2=maximum(abs.([c1, c2]))

tt = 1
nt = length(tecco)
function make_zonal_gif(exps, nt::Int64)
    a = Animation()
    for tt in 1:nt
        ps = []
        for ex in exps
            var = θ_zonal[ex][:, :, tt] .- θ_zonal["iter0_bulkformula"][:, :, tt]
            horiz_avg = @views reverse(var, dims = 1)
            jcf = Plots.contourf(X, Y, horiz_avg,
            xlabel = "latitude [º]",
            xticks = -70:20:70,
            ylabel = "depth [m]", 
            linewidth = 0.1,
            clim = (c1, c2),
            c = :balance, 
            colorbar_title= "°C",
            title = ex, 
            titlefontsize = 12)
            push!(ps, jcf)
        end
        l = @layout [a;b;c]
        p = Plots.plot(ps[1], ps[2], ps[3], size=(1400,1400), 
        margin = 4mm, link = :both,
        plot_title =  "Zonally averaged θ Anomaly "* region * "\n Time: " * 
        string(round(tecco[tt], digits = 6)),
        plot_titlefontsize	= 13, plot_titlevspan = 0.075, layout = l)
        frame(a, p)
    end
    gif(a, plotsdir() * "/OHC_Divergence/Gifs/" * "θ_zonal_avg_anom_" * 
        region * "_" * suffix * "_.gif", fps = 10)
end

function make_zonal_pic(exps, nt::Int64)
    a = Animation()
    ps = []
    for ex in exps
        baseline = θ_zonal["iter0_bulkformula"][:, :, nt]
        horiz_avg = @views reverse(θ_zonal[ex][:, :, nt] .- baseline, dims = 1)
        jcf = Plots.contourf(X, Y, horiz_avg,
        xlabel = "latitude [º]",
        xticks = -70:20:70,
        ylabel = "depth [m]", 
        linewidth = 0.1,
        clim = (c1, c2),
        c = :balance, 
        colorbar_title= "°C",
        title = ex, 
        titlefontsize = 12)
        push!(ps, jcf)
    end
    l = @layout [a;b;c]
    p = Plots.plot(ps[1], ps[2], ps[3], size=(1400,2000), 
    margin = 4mm, link = :both,
    plot_title =  "Zonally averaged θ Anomaly "* region * "\n Time: " * 
    string(round(tecco[nt], digits = 6)),
    plot_titlefontsize	= 13, plot_titlevspan = 0.075,  thickness_scaling = 2, layout = l)
    savefig(p, plotsdir() * "/OHC_Divergence/Gifs/" * "θ_zonal_avg_anom_$nt" * 
        region * "_" * suffix * "_.pdf")
end

function make_zonal_gif_iter129(nt::Int64)
    a = Animation()
    for tt in 1:nt
        var = θ_zonal["iter129_bulkformula"][:, :, tt] .- θ_zonal["iter0_bulkformula"][:, :, tt]
        horiz_avg = @views reverse(var, dims = 1)
        jcf = Plots.contourf(X, Y, horiz_avg,
        xlabel = "latitude [º]",
        xticks = -70:20:70,
        ylabel = "depth [m]", 
        linewidth = 0.1,
        clim = (c1, c2),
        c = :balance, 
        colorbar_title= "°C",
        title = L"\theta^{129}", 
        titlefontsize = 12)

        p = Plots.plot(jcf, size=(1400,1400), 
        margin = 4mm, link = :both,
        plot_title =  "Zonal avg. Pacific θ Anomaly " * "\n Time: " * 
        string(round(tecco[tt], digits = 6)),
        plot_titlefontsize	= 13, plot_titlevspan = 0.075)
        frame(a, p)
    end
    gif(a, plotsdir() * "/OHC_Divergence/Gifs/" * "θ_zonal_avg_anom_iter129_" * 
        region * "_" * suffix * "_.gif", fps = 10)
end

make_zonal_gif_iter129(nt)
make_zonal_gif(collect(keys(shortnames)), nt)

make_zonal_pic(exps, 1)
make_zonal_pic(exps, nt)
# savefig(plotsdir() * "/OHC_Divergence/" * "θDepthTS_" * region * suffix * ".png")

