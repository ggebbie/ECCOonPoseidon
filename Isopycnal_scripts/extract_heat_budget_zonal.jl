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
 
sig2grid = sigma2grid()
nσ = length(sig2grid)

region = "PAC56"; 
ocean_mask = OHC_helper.wet_pts(Γ)
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not")
areas3D = MeshArray(γ,Float32,nσ)
for ff=1:5, k=1:nσ
    areas3D.f[ff, k] .= area.f[ff] .* PAC_msk.f[ff]
end
ΔV = zonal_sum_3D(areas3D)
ΔV[ΔV .== 0] .= NaN

expname = "iter129_bulkformula"

tecco = 1992+1/24:1/12:2018
nz = 50; tstart = 12*3; tstop = 312; nt = tstop - tstart + 1;
E,F = trend_matrices(tecco[tstart:tstop])

get_datafilelist(x) = filter(x -> occursin("data",x), 
                             searchdir(runpath[expname]*"sigma2/", x * "_on_sigma2"))

function filter_heat_budget(savename, 
    diagpath::Dict{String, String}, 
    expname::String, 
    γ::gcmgrid)

    # Get list of files for salinity on sigma1

    dflist_TH        = get_datafilelist("THETA")
    dflist_THADVx    = get_datafilelist("ADVx_TH")
    dflist_THADVy    = get_datafilelist("ADVy_TH")
    dflist_THADVr    = get_datafilelist("ADVr_TH")
    dflist_THDFx     = get_datafilelist("DFyE_TH")
    dflist_THDFy     = get_datafilelist("DFxE_TH")
    dflist_THDFrE    = get_datafilelist("DFrE_TH")
    dflist_THDFrI    = get_datafilelist("DFrI_TH")

    var_names = ["dθ", "κxyθ", "uvθ", "wθ", "κzθ"]
    vars = Dict(varname => zeros(size(ΔV)) for varname in var_names)
    κθ_conv3D = MeshArray(γ,Float32,nσ); uθ_conv3D = MeshArray(γ,Float32,nσ);
    i = 0 
    @time for tt = tstart:tstop
        i+=1 
        println(tt)
        
        #horizontal convergences
        θ_gcm = γ.read(runpath[expname]*"sigma2/"*dflist_TH[tt],MeshArray(γ,Float32,nσ)); 
        κθx = γ.read(runpath[expname]*"sigma2/"*dflist_THDFx[tt],MeshArray(γ,Float32,nσ)); 
        κθy = γ.read(runpath[expname]*"sigma2/"*dflist_THDFy[tt],MeshArray(γ,Float32,nσ)); 
        uθx = γ.read(runpath[expname]*"sigma2/"*dflist_THADVx[tt],MeshArray(γ,Float32,nσ)); 
        uθy = γ.read(runpath[expname]*"sigma2/"*dflist_THADVy[tt],MeshArray(γ,Float32,nσ)); 

        for ijk in eachindex(κθx)
            κθx.f[ijk][isnan.(κθx.f[ijk])] .= 0.0
            κθy.f[ijk][isnan.(κθy.f[ijk])] .= 0.0
            uθx.f[ijk][isnan.(uθx.f[ijk])] .= 0.0
            uθy.f[ijk][isnan.(uθy.f[ijk])] .= 0.0
        end
        OHC_helper.calc_UV_conv3D!(κθx, κθy, κθ_conv3D); 
        OHC_helper.calc_UV_conv3D!(uθx, uθy, uθ_conv3D);

        #vertical convergences
        wθ_conv = γ.read(runpath[expname]*"sigma2/"*dflist_THADVr[tt],MeshArray(γ,Float32,nσ)); 
        κzθE_conv = γ.read(runpath[expname]*"sigma2/"*dflist_THDFrE[tt],MeshArray(γ,Float32,nσ)); 
        κzθI_conv = γ.read(runpath[expname]*"sigma2/"*dflist_THDFrI[tt],MeshArray(γ,Float32,nσ)); 

        for ijk in eachindex(wθ_conv)
            wθ_conv.f[ijk][isnan.(wθ_conv.f[ijk])] .= 0.0
            κzθE_conv.f[ijk][isnan.(κzθE_conv.f[ijk])] .= 0.0
            κzθI_conv.f[ijk][isnan.(κzθI_conv.f[ijk])] .= 0.0

        end

        κzθ_conv = κzθE_conv .+ κzθI_conv; 

        for ff=1:5, k=1:nz-1
            wθ_conv.f[ff, k] .= wθ_conv.f[ff, k+1] .- wθ_conv.f[ff, k] #in - out 
            κzθ_conv.f[ff, k] .= κzθ_conv.f[ff, k+1] .- κzθ_conv.f[ff, k] #does not do correct heatflux for z = 50 
        end

        θ_gcm_zonal = zonal_avg_3D(θ_gcm)
        wθ_conv = wθ_conv .* PAC_msk
        κzθ_conv = κzθ_conv .* PAC_msk
        κθ_conv3D = κθ_conv3D .* PAC_msk; 
        uθ_conv3D = uθ_conv3D .* PAC_msk
        
        #
        vars["dθ"]   .+= F[2,i] .* θ_gcm_zonal ./ 3.156e+7 #C / year -> C / s
        vars["κxyθ"] .+= (zonal_sum_3D(κθ_conv3D) ./ ΔV) ./nt #are the units here correct?.. probably not.. 
        vars["uvθ"]  .+= (zonal_sum_3D(uθ_conv3D) ./ ΔV) ./nt
        vars["wθ"]   .+= (zonal_sum_3D(wθ_conv) ./ ΔV) ./nt
        vars["κzθ"]  .+= (zonal_sum_3D(κzθ_conv) ./ ΔV) ./nt

    end

    return vars
end

sum_fluxes(d) = sum(d[varn] for varn in ["κxyθ", "uvθ", "wθ", "κzθ"])
expname = "iter129_bulkformula"
svename = datadir("sigma2/THETA_sigma2_budget_zonal_" * region * "_" * expname *".jld2")
dθ = filter_heat_budget(svename, diagpath, expname, γ)
jldsave(svename, dθ = dθ)
