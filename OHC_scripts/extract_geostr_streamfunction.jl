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

"""
    function getUeVe!
    computes velocities resulting Geoman Transport and compensating return flow
    U, V : empty arrays to be filled
    η : (unrotated) windstress on sea surface 
    finv : 1/f
    Hinv : 1 / (total depth of column at (x, y))
    hinv : 1 / (depth of surface cell at (x, y))

"""
function getUVgeo!(U::MeshArrays.gcmarray{T, 2, Matrix{T}}, V::MeshArrays.gcmarray{T, 2, Matrix{T}}, 
    η::MeshArrays.gcmarray{T, 1, Matrix{T}}, 
    Φ_hyd::MeshArrays.gcmarray{T, 2, Matrix{T}},
    iDXC::MeshArrays.gcmarray{T, 1, Matrix{T}}, 
    iDYC::MeshArrays.gcmarray{T, 1, Matrix{T}},
    finv::MeshArrays.gcmarray{T, 1, Matrix{T}}) where {T<:Real}
    g = 0.0; ρ0 = 1
    P_anom = Φ_hyd
    ∂xP = MeshArray(γ,Float32,50); ∂yP = MeshArray(γ,Float32,50); 
    for z=1:50
        P_anom_z = P_anom[:, z]
        ∂xP_temp, ∂yP_temp = gradient(P_anom_z, iDXC, iDYC) 
        ∂xP.f[:, z] .= ∂xP_temp.f; ∂yP.f[:, z] .= ∂yP_temp.f
    end    
    ∂xη, ∂yη = gradient(η, iDXC, iDYC) 
    for ff in 1:5
        ∂xη_ff = ∂xη.f[ff]; ∂yη_ff = ∂yη.f[ff]
        ∂xP_ff = ∂xP.f[ff]; ∂yP_ff = ∂yP.f[ff]
        finv_ff = finv[ff]
        for k in 1:50
            V.f[ff, k] .=  finv_ff .* (g.*∂xη_ff .+ ∂xP_ff)
            U.f[ff, k] .= -finv_ff .* (g.*∂yη_ff .+ ∂yP_ff)
        end
    end
end


"""
    function getUeVe!
    computes velocities resulting Geoman Transport and compensating return flow
    U, V : empty arrays to be filled
    η : (unrotated) windstress on sea surface 
    finv : 1/f
    Hinv : 1 / (total depth of column at (x, y))
    hinv : 1 / (depth of surface cell at (x, y))

"""
function getUVgeo_P0!(U::MeshArrays.gcmarray{T, 2, Matrix{T}}, V::MeshArrays.gcmarray{T, 2, Matrix{T}}, 
    η::MeshArrays.gcmarray{T, 1, Matrix{T}}, 
    H::MeshArrays.gcmarray{T, 1, Matrix{T}}, Hinv::MeshArrays.gcmarray{T, 1, Matrix{T}}, 
    z::Vector{Float64},
    finv::MeshArrays.gcmarray{T, 1, Matrix{T}},
    iDXC::MeshArrays.gcmarray{T, 1, Matrix{T}}, iDYC::MeshArrays.gcmarray{T, 1, Matrix{T}}) where {T<:Real}
    
    g = 9.81
    HηH = H .+ η; HηH = Hinv .*  HηH
    ∂xHηH, ∂yHηH = gradient(HηH, iDXC, iDYC) 

    zstar = MeshArray(γ,Float32,50)
    for (ij, k) in eachindex(zstar)
        @inbounds zstar.f[ij, k] .= H.f[ij] .* (z[k] .- η.f[ij]) .* inv.(H.f[ij] .+  η.f[ij])        
    end
    
    for ff in 1:5
        ∂xHηH_ff = ∂xHηH.f[ff]; ∂yHηH_ff = ∂yHηH.f[ff]; finv_ff = finv[ff]
        for k in 1:50
            V.f[ff, k] .= zstar[ff] .* finv_ff .* g .* ∂xHηH_ff
            U.f[ff, k] .= -zstar[ff] .* finv_ff .* g .* ∂yHηH_ff
        end
    end
end

function extract_meridionalΨGeo(expname,diagpath, Γ, γ, H, Hinv, finv, z,mask)
    #exf_zflux_set1
    # fileroot2 = "state_3d_set2"
    # filelist = searchdir(diagpath[expname],fileroot2) # first filter for state_3d_set1
    # datafilelist2  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    fileroot1 = "state_2d_set1"
    filelist = searchdir(diagpath[expname],fileroot1) # first filter for state_3d_set1
    datafilelist1  = filter(x -> occursin("data",x),filelist) # second filter for "data"

    LC=LatitudeCircles(-89.0:89.0,Γ)
    nz=size(Γ.hFacC,2); nl=length(LC); nt = length(datafilelist1)
    Ψs = zeros(nt, nl, nz)
    iDXC=MeshArray(γ,Float32,1)[:, 1]; iDYC=MeshArray(γ,Float32,1)[:, 1]
    for a=1:5
        iDXC.f[a].= inv.(Γ.DXC.f[a])
        iDYC.f[a].= inv.(Γ.DYC.f[a])
    end
    hFacC =  Γ.hFacC .* mask
    Threads.@threads for tt=1:5
        #fname2 = datafilelist2[tt]
        # Φ_hyd= γ.read(diagpath[expname]*fname2,MeshArray(γ,Float32,100))[:, 51:end] #hydrostatic pressure
        fname1 = datafilelist1[tt]
        η = γ.read(diagpath[expname]*fname1,MeshArray(γ,Float32,1))[:, 1] #hydrostatic pressure
        U = MeshArray(γ,Float32,50); V = MeshArray(γ,Float32,50);
        getUVgeo_P0!(U, V, η, H, Hinv, z, finv, iDXC, iDYC)
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
        ψ = ov ; ψ[ψ.==0.0].=NaN
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
X=(collect(-89.0:89.0)); Y=reverse(z); #coordinate variables
H = smush(get_cell_depths(msk, ΔzF, Γ.hFacC)); H = Float32.(H)
Hinv = deepcopy(H); Hinv[findall(Hinv .== 0)] = Inf
Hinv = 1 ./ Hinv; Hinv = Float32.(Hinv)

f0(ϕ) = abs(ϕ) > 2 ? (2 * (2π / 86400) * sind(ϕ)) : Inf32
fs = MeshArray(γ,Float32,1)[:, 1]; fs.f .= map.(x -> f0.(x), ϕ.f)
finv = 1 ./ fs; finv = Float32.(finv) 

@time for (key,values) in shortnames
    expname = key
    println(expname)
    Ψ_exp[expname] = extract_meridionalΨGeo(expname, diagpath, Γ, γ, H, Hinv, finv, z, msk)
    Ψ_exp_mean[expname] =  dropdims(mean(Ψ_exp[expname], dims = 1), dims =1)
end

clims = 1e-6 .* extrema(Ψ_exp_mean)
clims = (-maximum(abs.(clims)), maximum(abs.(clims)))
Ψ_sym = L"\Psi^{Geo}"
Ψ_exp_labels = [L"\Psi_{129}", L"\Psi_{\Delta F}",
                L"\Psi_{\Delta T}", L"\Psi_{0}"]
p = plot_Ψ_mean(X, Y, Ψ_exp_mean, Ψ_exp_labels, clims, Ψ_sym)
savefig(p,plotsdir() * "/OHC_Divergence/MassTransport/" * "ΨGeomean_"* 
region * ".png")

Ψ_exp_mean_anom = Dict(key => Ψ_exp_mean[key] .- Ψ_exp_mean["iter0_bulkformula"] 
                    for key in keys(Ψ_exp_mean))
Ψ̄_anom_labels = [exp * L"- \Psi_{0}" for exp in Ψ_exp_labels]      
clims = 1e-6 .* extrema(Ψ_exp_mean_anom)
clims = (-maximum(abs.(clims)), maximum(abs.(clims)))
Ψ_sym = L"\Psi^{Geo}" * "Anomalies"
p = plot_Ψ_mean(X, Y, Ψ_exp_mean_anom, Ψ̄_anom_labels, clims, Ψ_sym)
savefig(p,plotsdir() * "/OHC_Divergence/MassTransport/" * "ΨGeomean_0anom_"* 
region * ".png")

Ψ_exp_mean_anom = Dict(key => Ψ_exp_mean[key] .- Ψ_exp_mean["iter129_bulkformula"] 
                    for key in keys(Ψ_exp_mean))
Ψ̄_anom_labels = [exp * L"- \Psi_{129}" for exp in Ψ_exp_labels]      
clims = 1e-6 .* extrema(Ψ_exp_mean_anom)
clims = (-maximum(abs.(clims)), maximum(abs.(clims)))
Ψ_sym = L"\Psi_{Geo}" * "Anomalies"
p = plot_Ψ_mean(X, Y, Ψ_exp_mean_anom, Ψ̄_anom_labels, clims, Ψ_sym)
savefig(p,plotsdir() * "/OHC_Divergence/MassTransport/" * "ΨGeomean_129anom_"* 
region * ".png")
