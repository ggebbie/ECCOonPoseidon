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
msk = ocean_mask;
(msk == ocean_mask) && (region = "GLOB")
tecco = 1992+1/24:1/12:2018
"""
    function extract_sst34
    extract by reading multiple files
"""
function getUeVe!(U::MeshArrays.gcmarray{T, 2, Matrix{T}}, V::MeshArrays.gcmarray{T, 2, Matrix{T}}, 
    τx::MeshArrays.gcmarray{T, 1, Matrix{T}}, τy::MeshArrays.gcmarray{T, 1, Matrix{T}}, 
    finv::MeshArrays.gcmarray{T, 1, Matrix{T}}, 
    ρ0::MeshArrays.gcmarray{T, 1, Matrix{T}}, Hinv::MeshArrays.gcmarray{T, 1, Matrix{T}}, 
    hinv::MeshArrays.gcmarray{T, 1, Matrix{T}}) where {T<:Real}
    for ff in 1:5
        τxff = τx.f[ff]; τyff = τy.f[ff]
        fρ0inv = finv[ff] .* inv(1030)
        Hinvff = Hinv.f[ff]; hinvff = hinv.f[ff]
        for k in 1:50
            V.f[ff, k] .= ((hinvff .* (k == 1)) .- Hinvff ) .* τxff
            V.f[ff, k] .=  fρ0inv .* V.f[ff, k]

            U.f[ff, k] .= ((hinvff .* (k == 1)) .- Hinvff ).* τyff
            U.f[ff, k] .=  -fρ0inv .* U.f[ff, k]
        end
    end
end
function extract_meridionalΨEK(expname,diagpath,fileroot, Γ, γ, Hinv, finv, mask)
    #exf_zflux_set1
    filelist = searchdir(diagpath[expname],fileroot) # first filter for state_3d_set1
    datafilelist  = filter(x -> occursin("data",x),filelist)[24:288] # second filter for "data"
    LC=LatitudeCircles(-89.0:89.0,Γ)
    nz=size(Γ.hFacC,2); nl=length(LC); nt = length(datafilelist)

    Ψs = zeros(nt, nl, nz)

    ρ0 = MeshArray(γ,Float32,1)[:, 1]; fill!(ρ0, 1027)
    hinv = MeshArray(γ,Float32,1)[:, 1]; fill!(hinv, inv(10))
    
    Threads.@threads for tt=1:nt
        fname = datafilelist[tt]
        EXF = γ.read(diagpath[expname]*fname,MeshArray(γ,Float32,15))
        EXF = EXF .* mask; τx = EXF[:, 14]; τy = EXF[:, 15]
        U = MeshArray(γ,Float32,50); V = MeshArray(γ,Float32,50);
        getUeVe!(U, V, τx, τy, finv, ρ0, Hinv, hinv);
        hFacC =  Γ.hFacC
        U = U .* hFacC; V = V .* hFacC 
        (Utr,Vtr) = UVtoTransport(U,V,Γ)
        ov=Array{Float64,2}(undef,nl,nz)
        for z=1:nz
            @inbounds Uz = Utr[:,z] 
            @inbounds Vz = Vtr[:,z]
            for l=1:nl
                @inbounds ov[l,z] = ThroughFlowDim(Uz,Vz, LC[l])            
            end
        end
        ov=reverse(cumsum(reverse(ov,dims=2),dims=2),dims=2)
        ψ =ov ; ψ[ψ.==0.0].=NaN
        Ψs[tt, :, :] .= ψ
        GC.safepoint()
        println(tt)
    end
    return Ψs
end



function plot_zonal_contours(X, Y, zonal_var, clims, title)
    jcf = Plots.contourf(X, Y, reverse(zonal_var, dims = 1),
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

function plot_all_zonal_Ψ(X, Y, Ψ_zonal, nt::Int64, tecco, clims)
    Ψ_sym = L"\Psi_{Ek}" * " "
    exps = collect(keys(Ψ_zonal))
    a = Animation()
    ps = Vector{Any}(missing, length(exps))
    for tt in 1:nt
        for (i, ex) in enumerate(exps)
            ps[i] = plot_zonal_contours(X, Y, 1e-6 .* Ψ_zonal[ex][tt, :, :]', 
            clims, ex)
        end
        l = @layout [a;b;c; d]
        p = Plots.plot(ps[1], ps[2], ps[3], ps[4], size=(1401, 1867), link = :both,
        plot_title =  "Zonally averaged " * Ψ_sym * region * "\n Time: " * 
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


fileroot = "state_2d_set1"
Ψ_exp = Dict()
Ψ_exp_mean = Dict()
X=(collect(-89.0:89.0)); Y=reverse(z); #coordinate variables
Hinv = smush(get_cell_depths(msk, ΔzF, Γ.hFacC)); Hinv[findall(Hinv .== 0)] = Inf
Hinv = 1 ./ Hinv; Hinv = Float32.(Hinv)

f(ϕ) = abs(ϕ) > 2 ? (2 * (2π / 86400) * sin(deg2rad(ϕ))) : Inf32
fs = MeshArray(γ,Float32,1)[:, 1]; fs.f .= map.(x -> f(x), ϕ)
finv = 1 ./ fs; finv = Float32.(finv) 

@time for (key,values) in shortnames
    expname = key
    println(expname)
    Ψ_exp[expname] = extract_meridionalΨEK(expname,diagpath,fileroot, Γ, γ,
                                           Hinv, finv, msk)
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
savefig(p,plotsdir() * "/OHC_Divergence/MassTransport/" * "ΨEKmean_zonal_"* 
region * ".png")

Ψ_exp_mean_anom = Dict(key => Ψ_exp_mean[key] .- Ψ_exp_mean["iter0_bulkformula"] 
                    for key in keys(Ψ_exp_mean))
clims = 1e-6 .* extrema(Ψ_exp_mean_anom)
clims = (-maximum(abs.(clims)), maximum(abs.(clims)))

for (i, ex) in enumerate(keys(shortnames))
    ps[i] = plot_zonal_contours(X, Y, 1e-6 .* Ψ_exp_mean_anom[ex]', 
    clims, ex)
end
l = @layout [a;b;c]
p = Plots.plot(ps[1], ps[2], ps[3], size=(1401,1401), link = :both,
plot_title =  "Zonally averaged Ψ̄ Anomaly "* region,
plot_titlefontsize	= 17, plot_titlevspan = 0.1, layout = l)
savefig(p,plotsdir() * "/OHC_Divergence/MassTransport/" * "ΨEKmean_zonal_initanom_"* 
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