#analysis should complete within 12-13 minutes 
#using 12 threads 
# julia --threads=6 --project=@. ./extract_heat_budget_native.jl

include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, 
    BenchmarkTools, LaTeXStrings, PyCall
using .OHC_helper
import NaNMath as nm
import PyPlot as plt
include(srcdir("config_exp.jl"))
include("./IsopycnalHelpers.jl")

cmo = pyimport("cmocean.cm");

#using . is faster 
(ϕ,λ) = latlonC(γ)
area = readarea(γ)

runpath,diagpath = listexperiments(exprootdir());

# abbreviations for each experiment for labels, etc.
ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)
marks = expsymbols()
nexp = length(shortnames) # number of experiments

function calc_UV_conv3D_Isopycnal!(CONV::MeshArrays.gcmarray{T, 2, Matrix{T}}, uFLD::MeshArrays.gcmarray{T, 2, Matrix{T}}, 
    vFLD::MeshArrays.gcmarray{T, 2, Matrix{T}},
    inv_DXG::MeshArrays.gcmarray{Float64, 1, Matrix{Float64}}, 
    inv_DYG::MeshArrays.gcmarray{Float64, 1, Matrix{Float64}}) where T<:Real
        tmpU, tmpV = exch_UV_cs3D(uFLD,vFLD)
        for a in eachindex(uFLD)
            ff = a[1]
            (s1,s2)=size(uFLD.f[a])
            @inbounds tmpU1=view(tmpU.f[a],1:s1,1:s2)
            @inbounds tmpU2=view(tmpU.f[a],2:s1+1,1:s2)
            @inbounds tmpV1=view(tmpV.f[a],1:s1,1:s2)
            @inbounds tmpV2=view(tmpV.f[a],1:s1,2:s2+1)
            @inbounds CONV.f[a] .= ((tmpU1.-tmpU2) .* inv_DXG.f[ff]) .+ ((tmpV1.-tmpV2) .* inv_DYG.f[ff])
        end
end

region = "PAC56"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not")
cell_depths = OHC_helper.get_cell_depths(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = Float32.(OHC_helper.get_cell_volumes(area, cell_depths));

expname = "iter129_bulkformula"

tecco = 1992+1/24:1/12:2018
nz = 50; tstart = 12*3; tstop = 312; nt = tstop - tstart + 1;
E,F = trend_matrices(tecco[tstart:tstop])

sig2grid = sigma2grid()
nσ = length(sig2grid)
_, σidx = findmin(abs.(sig2grid .- 36.92));
ptidx = 100

get_datafilelist(x, expname) = filter(x -> occursin("data",x), 
                             searchdir(runpath[expname]*"sigma2/", x * "_on_sigma2"))

function filter_heat_budget(inv_DXG::MeshArrays.gcmarray{T, 1, Matrix{T}}, 
                            inv_DYG::MeshArrays.gcmarray{T, 1, Matrix{T}},
                            expname::String, γ::gcmgrid) where T<:Real
    
    dflist_TH        = get_datafilelist("THETA", expname)
    dflist_THADVx    = get_datafilelist("ADVx_TH", expname)
    dflist_THADVy    = get_datafilelist("ADVy_TH", expname)
    dflist_THADVr    = get_datafilelist("ADVr_TH", expname)
    dflist_THDFx     = get_datafilelist("DFxE_TH", expname)
    dflist_THDFy     = get_datafilelist("DFyE_TH", expname)
    dflist_THDFr     = get_datafilelist("DFr_TH", expname)
    dflist_P          = get_datafilelist("p", expname)

    var_names = ["θ", "κxyθ", "uvθ", "wθ", "κzθ"]
    vars = Dict(varname => zeros(10) for varname in var_names)
    κθ_conv3D = MeshArray(γ,Float32,nσ); uθ_conv3D = MeshArray(γ,Float32,nσ);
    wθ_conv = MeshArray(γ,Float32,nσ); kθr_conv = MeshArray(γ,Float32,nσ);

    @time for tt = 1:10
        println(tt)
        θ_gcm = γ.read(runpath[expname]*"sigma2/"*dflist_TH[tt],MeshArray(γ,Float32,nσ)); 
        κθx = γ.read(runpath[expname]*"sigma2/"*dflist_THDFx[tt],MeshArray(γ,Float32,nσ)); 
        κθy = γ.read(runpath[expname]*"sigma2/"*dflist_THDFy[tt],MeshArray(γ,Float32,nσ)); 
        uθ  = γ.read(runpath[expname]*"sigma2/"*dflist_THADVx[tt],MeshArray(γ,Float32,nσ)); 
        vθ  = γ.read(runpath[expname]*"sigma2/"*dflist_THADVy[tt],MeshArray(γ,Float32,nσ)); 
        kθr = γ.read(runpath[expname]*"sigma2/"*dflist_THDFr[tt],MeshArray(γ,Float32,nσ)); 
        wθ  = γ.read(runpath[expname]*"sigma2/"*dflist_THADVr[tt],MeshArray(γ,Float32,nσ)); 
        Pσ  = γ.read(runpath[expname]*"sigma2/"*dflist_P[tt],MeshArray(γ,Float32,nσ)); 
        # for a in eachindex(κθx)
        #     κθx.f[a][isnan.(κθx.f[a])] .= 0.0 
        #     κθy.f[a][isnan.(κθy.f[a])] .= 0.0 
        #     uθ.f[a][isnan.(uθ.f[a])] .= 0.0 
        #     vθ.f[a][isnan.(vθ.f[a])] .= 0.0 
        # end
        calc_UV_conv3D_Isopycnal!(uθ_conv3D, uθ, vθ, inv_DXG, inv_DYG)
        calc_UV_conv3D_Isopycnal!(κθ_conv3D, κθx, κθy, inv_DXG, inv_DYG)

        for ff=1:5, k=2:nσ-1
            wθ_conv.f[ff, k]  .= (wθ.f[ff, k+1] .- wθ.f[ff, k-1]) ./ (Pσ.f[ff, k+1] .- Pσ.f[ff, k-1])
            kθr_conv.f[ff, k] .= (kθr.f[ff, k+1] .- kθr.f[ff, k-1]) ./ (Pσ.f[ff, k+1] .- Pσ.f[ff, k-1])
        end

        vars["θ"][tt]   =  θ_gcm[4, σidx][ptidx] #C / year
        println(vars["θ"][tt])

        vars["κxyθ"][tt] = κθ_conv3D[4, σidx][ptidx] 
        vars["uvθ"][tt]  = uθ_conv3D[4, σidx][ptidx] 
        vars["wθ"][tt]   = wθ_conv[4, σidx][ptidx] 
        vars["κzθ"][tt]  = kθr_conv[4, σidx][ptidx] 

    end

    return vars
end

sum_fluxes(d) = sum(d[varn][:] for varn in ["κxyθ", "uvθ", "wθ", "κzθ"])

expname = "iter129_bulkformula"
inv_DXG = inv_ma(deepcopy(Γ.DXG)); inv_DYG = inv_ma(deepcopy(Γ.DYG))
dθ = filter_heat_budget(inv_DXG, inv_DYG, expname, γ)

dθ["dθ"] = (dθ["θ"][2:end] .- dθ["θ"][1:end-1]) ./ 2.628e+6
dθ["sum"] = (sum_fluxes(dθ))
tecco_int = (tecco[2:nt] .+ tecco[1:nt-1]) ./ 2

fig, axes = plt.subplots(figsize=(10,5))
axes.plot(tecco[1:10], dθ["sum"], label = "reconstruction")
axes.plot(tecco_int[1:9], dθ["dθ"], label = "true")
axes.legend(); axes.set_title("dθ/dt")
fig
println(mean(dθ["dθ"])) #off by an order of magnitude again...? 
println(mean(dθ["sum"]))