#analysis should complete within 12 minutes 
#using 12 threads 
include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, BenchmarkTools
using .OHC_helper
using Plots
using ColorSchemes
import NaNMath as nm
include(srcdir("config_exp.jl"))

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

runpath,diagpath = listexperiments(exprootdir());
ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)
marks = expsymbols()
 
ocean_mask = wet_pts(Γ)
region = "NPAC"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "full")
msk = PAC_msk;
(msk == ocean_mask) && (region = "GLOB")
tecco = 1992+1/24:1/12:2018
function plot_zonal_contours(X, Y, zonal_var, clims, title)
    jcf = Plots.contourf(X, Y, reverse(zonal_var, dims = 1),
    xlabel = "latitude [º]",
    xticks = -70:20:70,
    ylabel = "depth [m]", 
    linewidth = 0.1,
    levels = 16,
    clim = clims,
    c = :delta, 
    colorbar_title= "Sv",
    title = title, 
    thickness_scaling = 1.5)
    return jcf
end

function plot_all_zonal_Ψ(X, Y, Ψ_zonal, nt::Int64, tecco, clims)
    exps = collect(keys(Ψ_zonal))
    a = Animation()
    ps = Vector{Any}(missing, length(exps))
    for tt in 1:nt
        for (i, ex) in enumerate(exps)
            ps[i] = plot_zonal_contours(X, Y, 1e-6 .* Ψ_zonal[ex][tt, :, :]', 
            clims, ex)
        end
        l = @layout [a;b;c]
        p = Plots.plot(ps[1], ps[2], ps[3], size=(1401,1401), link = :both,
        plot_title =  "Zonally averaged Ψ "* region * "\n Time: " * 
        string(round(tecco[tt], digits = 6)),
        plot_titlefontsize	= 17, plot_titlevspan = 0.075, layout = l)
        frame(a, p)
    end
    return a
end

function plot_Ψ(X, Y, Ψ_zonal, expname, nt::Int64, tecco, clims)
    a = Animation()
    for tt in 1:nt
        ps = plot_zonal_contours(X, Y, 1e-6 .* Ψ_zonal[expname][tt, :, :]', 
        clims, expname)
        p = Plots.plot(ps, size=(1400, 700), link = :both,
        plot_title =  "Zonally averaged Ψ "* region * "\n Time: " * 
        string(round(tecco[tt], digits = 6)),
        plot_titlefontsize = 12, plot_titlevspan = 0.1)
        frame(a, p)
    end
    return a
end


Ψ_exp = Dict()
Ψ_exp_mean = Dict()
X=(collect(-89.0:89.0)); Y=reverse(z); #coordinate variables

@time for (key,values) in shortnames
    expname = key
    println(expname)
    Ψ_exp[expname] = extract_meridionalΨ(expname,diagpath, Γ, γ, msk)
    Ψ_exp_mean[expname] =  dropdims(mean(Ψ_exp[expname], dims = 1), dims =1)

end

clims = 1e-6 .* extrema(Ψ_exp_mean)
clims = (-maximum(abs.(clims)), maximum(abs.(clims)))

ps = Vector{Any}(missing, length(keys(shortnames)))
for (i, ex) in enumerate(keys(shortnames))
    ps[i] = plot_zonal_contours(X, Y, 1e-6 .* Ψ_exp_mean[ex]', 
    clims, ex)
end
l = @layout [a;b;c]
p = Plots.plot(ps[1], ps[2], ps[3], size=(1401,1401), link = :both,
plot_title =  "Zonally averaged Ψ̄ "* region ,
plot_titlefontsize	= 17, plot_titlevspan = 0.1, layout = l)
savefig(p,plotsdir() * "/OHC_Divergence/MassTransport/" * "Ψmean_zonal_"* 
region * ".png")

Ψ_exp_mean_anom = Dict(key => Ψ_exp_mean[key] .- Ψ_exp_mean["iter0_bulkformula"] 
                    for key in keys(Ψ_exp_mean))
clims = 1e-6 .* extrema(Ψ_exp_mean_anom)
clims = (-maximum(abs.(clims)), maximum(abs.(clims)))

ps = Vector{Any}(missing, length(keys(shortnames)))
for (i, ex) in enumerate(keys(shortnames))
    ps[i] = plot_zonal_contours(X, Y, 1e-6 .* Ψ_exp_mean_anom[ex]', 
    clims, ex)
end
l = @layout [a;b;c]
p = Plots.plot(ps[1], ps[2], ps[3], size=(1401,1401), link = :both,
plot_title =  "Zonally averaged Ψ̄ Anomaly "* region,
plot_titlefontsize	= 17, plot_titlevspan = 0.1, layout = l)
savefig(p,plotsdir() * "/OHC_Divergence/MassTransport/" * "Ψmean_zonal_initanom_"* 
region * ".png")

# clims = 1e-6 .* extrema(Ψ_exp)
# nt = length(tecco)
# a = plot_all_zonal_Ψ(X, Y, Ψ_exp, nt, tecco, clims)
# mp4(a, plotsdir() * "/OHC_Divergence/Gifs/" * "Ψ_zonalavg_" * 
# region * ".mp4", fps = 5)

# expname = "iter129_bulkformula"
# explabel = "iter129"
# a = plot_Ψ(X, Y, Ψ_exp, expname, nt, tecco, clims)
# mp4(a, plotsdir() * "/OHC_Divergence/Gifs/" * "Ψ_zonalavg_" * explabel * "_" *  
# region * ".mp4", fps = 5)

# Ψ_exp_anom = Dict(key => Ψ_exp[key] .- Ψ_exp["iter0_bulkformula"] 
#                     for key in keys(Ψ_exp))
# clims = 1e-6 .* extrema(Ψ_exp_anom)
# clims = (-maximum(abs.(clims)), maximum(abs.(clims)))

# a = plot_all_zonal_Ψ(X, Y, Ψ_exp_anom, nt, tecco, clims)
# mp4(a, plotsdir() * "/OHC_Divergence/Gifs/" * "Ψ_zonalavg_iter0anom_" * 
# region * ".mp4", fps = 5)

# expname = "iter129_bulkformula"
# explabel = "iter129"
# a = plot_Ψ(X, Y, Ψ_exp_anom, expname, nt, tecco, clims)
# mp4(a, plotsdir() * "/OHC_Divergence/Gifs/" * "Ψ_zonalavg_iter0anom_" * explabel * "_" *  
# region * ".mp4", fps = 5)