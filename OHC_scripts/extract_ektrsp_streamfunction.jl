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
msk = PAC_msk;
(msk == ocean_mask) && (region = "GLOB")
tecco = 1992+1/24:1/12:2018

function extract_meridionalΨEK(expname,diagpath,fileroot, Γ, γ, Hinv, finv, mask)
    #exf_zflux_set1
    filelist = searchdir(diagpath[expname],fileroot) # first filter for state_3d_set1
    datafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    LC=LatitudeCircles(-89.0:89.0,Γ)
    nz=size(Γ.hFacC,2); nl=length(LC); nt = length(datafilelist)

    Ψs = zeros(nt, nl, nz)

    ρ0 = MeshArray(γ,Float32,1)[:, 1]; fill!(ρ0, 1027)
    hinv = MeshArray(γ,Float32,1)[:, 1]; fill!(hinv, inv(10))
    hFacC =  Γ.hFacC
    surfthick = hFacC[:,1]; surfthick[findall(surfthick .== 0)] = Inf32
    surfcell_inv = 1 ./ surfthick ; surfcell_inv = Float32.(surfcell_inv)
    hinv = hinv .* surfcell_inv
    Threads.@threads for tt=1:nt
        fname = datafilelist[tt]
        EXF = γ.read(diagpath[expname]*fname,MeshArray(γ,Float32,15))
        EXF = EXF .* mask; τx = EXF[:, 14]; τy = EXF[:, 15]
        U = MeshArray(γ,Float32,50); V = MeshArray(γ,Float32,50);
        getUVek!(U, V, τx, τy, finv, Hinv, hinv);
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
        Ψs[tt, :, :] .= -ψ
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

fileroot = "state_2d_set1"
Ψ_exp = Dict()
Ψ_exp_mean = Dict()
X=(collect(-89.0:89.0)); Y=reverse(z); #coordinate variables
Hinv = smush(get_cell_depths(msk, ΔzF, Γ.hFacC)); Hinv[findall(Hinv .== 0)] = Inf
Hinv = 1 ./ Hinv; Hinv = Float32.(Hinv)

f0(ϕ) = abs(ϕ) > 2 ? (2 * (2π / 86400) * sind(ϕ)) : Inf32
fs = MeshArray(γ,Float32,1)[:, 1]; fs.f .= map.(x -> f0.(x), ϕ.f)
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
Ψ_sym = L"\Psi^{Ek}"
Ψ_exp_labels = [L"\Psi^{\Delta F, \Delta T }", L"\Psi^{\Delta F}",
                L"\Psi^{\Delta T}", L"\Psi^{0}"]

p = plot_Ψ_mean(X, Y, Ψ_exp_mean, Ψ_exp_labels, clims, Ψ_sym)
savefig(p,plotsdir() * "/OHC_Divergence/MassTransport/" * "ΨEKmean2_"* 
region * ".png")

Ψ_exp_mean_anom = Dict(key => Ψ_exp_mean[key] .- Ψ_exp_mean["iter0_bulkformula"] 
                    for key in keys(Ψ_exp_mean))
Ψ̄_anom_labels = [exp * L"- \Psi^{0}" for exp in Ψ_exp_labels]      
clims = 1e-6 .* extrema(Ψ_exp_mean_anom)
clims = (-maximum(abs.(clims)), maximum(abs.(clims)))
Ψ_sym = L"\Psi^{Ek}" * "Anomalies"
p = plot_Ψ_mean(X, Y, Ψ_exp_mean_anom, Ψ̄_anom_labels, clims, Ψ_sym)
savefig(p,plotsdir() * "/OHC_Divergence/MassTransport/" * "ΨEKmean2_0anom_"* 
region * ".png")

Ψ_exp_mean_anom = Dict(key => Ψ_exp_mean[key] .- Ψ_exp_mean["iter129_bulkformula"] 
                    for key in keys(Ψ_exp_mean))
Ψ̄_anom_labels = [exp * L"- \Psi^{\Delta F, \Delta T}" for exp in Ψ_exp_labels]      
clims = 1e-6 .* extrema(Ψ_exp_mean_anom)
clims = (-maximum(abs.(clims)), maximum(abs.(clims)))
Ψ_sym = L"\Psi_{Ek}" * "Anomalies"
p = plot_Ψ_mean(X, Y, Ψ_exp_mean_anom, Ψ̄_anom_labels, clims, Ψ_sym)
savefig(p,plotsdir() * "/OHC_Divergence/MassTransport/" * "ΨEKmean2_129anom_"* 
region * ".png")
