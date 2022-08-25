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
c1 = minimum([nm.minimum(θ_zonal[ex]) for ex ∈ exps])
c2 = maximum([nm.maximum(θ_zonal[ex]) for ex ∈ exps])
tt = 1
nt = length(tecco)
function make_zonal_gif(exps, nt::Int64)
    a = Animation()
    for tt in 1:nt
        ps = []
        for ex in exps
            baseline = 0
            horiz_avg = @views reverse(θ_zonal[ex][:, :, tt] .- baseline, dims = 1)
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
    gif(a, plotsdir() * "/OHC_Divergence/Gifs/" * "θ_zonal_avg_" * 
        region * "_" * suffix * "_.gif", fps = 10)
end

exps = collect(keys(shortnames))
make_zonal_gif(exps, nt)

# savefig(plotsdir() * "/OHC_Divergence/" * "θDepthTS_" * region * suffix * ".png")

