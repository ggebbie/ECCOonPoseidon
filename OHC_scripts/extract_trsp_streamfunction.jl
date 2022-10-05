include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, BenchmarkTools, LaTeXStrings
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
region = "PAC"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "full")
msk = PAC_msk;
(msk == ocean_mask) && (region = "GLOB")
tecco = 1992+1/24:1/12:2018

X=collect(-89.0:89.0); Y=reverse(z); #coordinate variables

suffix = "sfctobot"
uplvl = -2e3; botlvl = -3e3
lvls = findall( botlvl .<= Y[:].<= uplvl)

function plot_zonal_contours(X, Y, zonal_var, clims, title, 
    color = :delta, levels = 30, c_labels = false    )
    jcf = Plots.contourf(X, Y, zonal_var,
    xlabel = "latitude [º]",
    xticks = -70:20:70,
    ylabel = "depth [m]", 
    linewidth = 1,
    levels = levels,
    clim = clims,
    c = color, 
    colorbar_title= "Sv",
    title = title, 
    thickness_scaling = 1.5,
    contour_labels = c_labels)
    return jcf
end

function plot_Ψ_mean(X, Y, Ψ_zonal, Ψ_labels, clims, Ψ_sym, lvls)
    ps = Vector{Any}(missing, length(keys(Ψ_zonal)))
    i = [0] 
    for ex in keys(Ψ_zonal)
        i .+= 1
        ps[i[1]] = plot_zonal_contours(X, Y[lvls], 1e-6 .* Ψ_zonal[ex][:, lvls]', 
        clims, Ψ_labels[i[1]])
    end
    l = @layout [a;b;c;d]
    p = Plots.plot(ps[1], ps[2], ps[3],ps[4],size=(1401,2101), link = :both,
    plot_title =  "Time mean "* Ψ_sym  * " "* region,
    plot_titlefontsize	= 17, plot_titlevspan = 0.025, layout = l)
    return p 
end
             
function plot_Ψ_std(X, Y, Ψ_zonal, Ψ_labels, clims, Ψ_sym, lvls)
    ps = Vector{Any}(missing, length(keys(Ψ_zonal)))
    i = [0] 
    for ex in keys(Ψ_zonal)
        i .+= 1
        ps[i[1]] = plot_zonal_contours(X, Y[lvls], 1e-6 .* Ψ_zonal[ex][:, lvls]', 
        clims, Ψ_labels[i[1]], :amp, 6, true)
    end        
    l = @layout [a;b;c;d]
    p = Plots.plot(ps[1], ps[2], ps[3],ps[4],size=(1401,2101), link = :both,
    plot_title =  Ψ_sym  * " "* region,
    plot_titlefontsize	= 17, plot_titlevspan = 0.025, layout = l)
    return p 
end
     

Ψ_exp_mean = Dict(); Ψ_exp_std = Dict()
Ψ_exp_labels = [L"\Psi^{\Delta F, \Delta T }", L"\Psi^{\Delta F}",
                L"\Psi^{\Delta T}", L"\Psi^{0}"]

@time for (key,values) in shortnames
    expname = key
    println(expname)
    Ψ_exp_mean[expname], Ψ_exp_std[expname] = extract_meridionalΨ̄(expname,diagpath, Γ, γ, msk)
end

Ψ_sym = L"\Psi"
clims = 1e-6 .* extrema(Ψ_exp_mean); 
clims = (-maximum(abs.(clims)), maximum(abs.(clims)))
p = plot_Ψ_mean(X, Y, Ψ_exp_mean, Ψ_exp_labels, clims, Ψ_sym, lvls)
savefig(p,plotsdir() * "/OHC_Divergence/MassTransport/" * "Ψmean2_"* 
region * "_" *suffix * ".png")

Ψ_sym = L"\Psi" * " Standard Deviation"
Ψ_std_exp_labels = "SD(" .* Ψ_exp_labels .* ")"
clims = 1e-6 .* extrema(Ψ_exp_std); 
clims = (0, maximum(clims))
p = plot_Ψ_std(X, Y, Ψ_exp_std, Ψ_std_exp_labels, clims, Ψ_sym, lvls)
savefig(p,plotsdir() * "/OHC_Divergence/MassTransport/" * "Ψstd2_"* 
region * "_" *suffix * ".png")

Ψ_sym = L"\Psi" * " Anomaly"
Ψ_exp_mean_anom = Dict(key => Ψ_exp_mean[key] .- Ψ_exp_mean["iter0_bulkformula"] 
                    for key in keys(Ψ_exp_mean))

Ψ̄_anom_labels = [exp * L"- \Psi^{0}" for exp in Ψ_exp_labels]      

clims = 1e-6 .* extrema(Ψ_exp_mean_anom)
clims = (-maximum(abs.(clims)), maximum(abs.(clims)))
p = plot_Ψ_mean(X, Y, Ψ_exp_mean_anom, Ψ̄_anom_labels, clims, Ψ_sym, lvls)
savefig(p,plotsdir() * "/OHC_Divergence/MassTransport/" * "Ψmean2_iter0anom_"* 
region * suffix * ".png")

Ψ_exp_mean_anom = Dict(key => Ψ_exp_mean[key] .- Ψ_exp_mean["iter129_bulkformula"] 
                    for key in keys(Ψ_exp_mean))
Ψ̄_anom_labels = [exp * L"- \Psi^{\Delta F, \Delta T}" for exp in Ψ_exp_labels]      

clims = 1e-6 .* extrema(Ψ_exp_mean_anom)
clims = (-maximum(abs.(clims)), maximum(abs.(clims)))
p = plot_Ψ_mean(X, Y, Ψ_exp_mean_anom, Ψ̄_anom_labels, clims, Ψ_sym, lvls)
savefig(p,plotsdir() * "/OHC_Divergence/MassTransport/" * "Ψmean2_iter129anom_"* 
region * suffix * ".png")