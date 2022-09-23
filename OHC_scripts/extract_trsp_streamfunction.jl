#analysis should complete within 12 minutes 
#using 12 threads 
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
region = "NPAC"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "full")
msk = ocean_mask;
(msk == ocean_mask) && (region = "GLOB")
tecco = 1992+1/24:1/12:2018

suffix = "sfctobot"
uplvl = +Inf; botlvl = -Inf
lvls = findall( botlvl .<= z[:].<= uplvl)


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

function plot_Ψ_mean(X, Y, Ψ_zonal, Ψ_labels, clims, Ψ_sym)
    ps = Vector{Any}(missing, length(keys(Ψ_zonal)))
    i = [0] 
    for ex in keys(Ψ_zonal)
        i .+= 1
        ps[i[1]] = plot_zonal_contours(X, Y, 1e-6 .* Ψ_zonal[ex]', 
        clims, Ψ_labels[i[1]])
    end
    l = @layout [a;b;c;d]
    p = Plots.plot(ps[1], ps[2], ps[3],ps[4],size=(1401,2101), link = :both,
    plot_title =  "Time mean "* Ψ_sym  * " "* region,
    plot_titlefontsize	= 17, plot_titlevspan = 0.025, layout = l)
    return p 
end

Ψ_exp = Dict()
Ψ_exp_mean = Dict()
X=(collect(-89.0:89.0)); Y=reverse(z[lvls]); #coordinate variables

@time for (key,values) in shortnames
    expname = key
    println(expname)
    Ψ_exp[expname] = extract_meridionalΨ(expname,diagpath, Γ, γ, msk)
    Ψ_exp_mean[expname] =  dropdims(mean(Ψ_exp[expname], dims = 1), dims =1)
    Ψ_exp_mean[expname] = Ψ_exp_mean[expname][:, lvls]
end

Ψ_sym = L"\Psi"
Ψ_exp_labels = [L"\Psi^{\Delta F, \Delta T }", L"\Psi^{\Delta F}",
                L"\Psi^{\Delta T}", L"\Psi^{0}"]

clims = 1e-6 .* extrema(Ψ_exp_mean)
clims = (-maximum(abs.(clims)), maximum(abs.(clims)))
p = plot_Ψ_mean(X, Y, Ψ_exp_mean, Ψ_exp_labels, clims, Ψ_sym)
savefig(p,plotsdir() * "/OHC_Divergence/MassTransport/" * "Ψmean_"* 
region * ".png")

Ψ_exp_mean_anom = Dict(key => Ψ_exp_mean[key] .- Ψ_exp_mean["iter0_bulkformula"] 
                    for key in keys(Ψ_exp_mean))
Ψ̄_anom_labels = [exp * L"- \Psi^{0}" for exp in Ψ_exp_labels]      

clims = 1e-6 .* extrema(Ψ_exp_mean_anom)
clims = (-maximum(abs.(clims)), maximum(abs.(clims)))
p = plot_Ψ_mean(X, Y, Ψ_exp_mean_anom, Ψ̄_anom_labels, clims, Ψ_sym)
savefig(p,plotsdir() * "/OHC_Divergence/MassTransport/" * "Ψmean_iter0anom_"* 
region * suffix * ".png")

Ψ_exp_mean_anom = Dict(key => Ψ_exp_mean[key] .- Ψ_exp_mean["iter129_bulkformula"] 
                    for key in keys(Ψ_exp_mean))
Ψ̄_anom_labels = [exp * L"- \Psi^{\Delta F, \Delta T}" for exp in Ψ_exp_labels]      

clims = 1e-6 .* extrema(Ψ_exp_mean_anom)
clims = (-maximum(abs.(clims)), maximum(abs.(clims)))
p = plot_Ψ_mean(X, Y, Ψ_exp_mean_anom, Ψ̄_anom_labels, clims, Ψ_sym)
savefig(p,plotsdir() * "/OHC_Divergence/MassTransport/" * "Ψmean_iter129anom_"* 
region * suffix * ".png")