#analysis should complete within 12-13 minutes 
#using 12 threads 
# julia --threads=6 --project=@. ./extract_heat_budget_native.jl

include("../../src/intro.jl")
include("../../src/OHC_helper.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, 
    BenchmarkTools, LaTeXStrings, PyCall
using .OHC_helper
import NaNMath as nm
import PyPlot as plt

include(srcdir("config_exp.jl"))

cmo = pyimport("cmocean.cm");
@pyimport seaborn as sns;
sns.set_theme(context = "poster", style = "ticks", font_scale = 1.0,
              palette = sns.color_palette("colorblind"));
colors =  sns.color_palette("deep")[1:4]

(ϕ,λ) = latlonC(γ)
area = readarea(γ)
ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)
tecco = 1992+1/24:1/12:2018

ocean_mask = OHC_helper.wet_pts(Γ)
region = "NPAC"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; region, extent = "not", include_bering = false)
cell_depths = OHC_helper.get_cell_depths(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = Float32.(OHC_helper.get_cell_volumes(area, cell_depths));
area_mask = area .* PAC_msk

function vertical_flux_profile(ds::MeshArray)
    nz = size(ds, 2)
    vol_avg = zeros(Float32, nz)

    for ff=1:5, k=1:nz
        vol_avg[k] += Float32(sum(ds[ff, k]))
    end
    return vol_avg[:]
end

W_dict = Dict(); N_dict = Dict()
runpath,diagpath = listexperiments(exprootdir());
mskC, mskW, mskS = OHC_helper.get_msk(Γ)
cs, sn = OHC_helper.get_cs_and_sn(γ)
ϕ_min_mask, ϕ_max_mask = OHC_helper.get_ϕ_max_min_mask(region, Γ, λ, ϕ, basins, basin_list)

function filter_volume_budget(diagpath::Dict{String, String}, expname::String, γ::gcmgrid) 
    filelist = searchdir(diagpath[expname],"trsp_3d_set1") # first filter for state_3d_set1
    datafilelist_uvw  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    nz = 50; nt = 312
    We = zeros(Float32, nz, nt); Wbol = zeros(Float32, nz, nt)
    Wresid = zeros(Float32, nz, nt)

    Ne = zeros(Float32, nz, nt); Nbol = zeros(Float32, nz, nt)
    Nresid = zeros(Float32, nz, nt)

    Wdict = Dict(); Ndict = Dict()
    @time for tt = 1:nt
        println(tt)
        # u, v, w = OHC_helper.extract_velocities(diagpath, expname , datafilelist_uvw[tt], γ)
        u, v, w, Ub, Vb, Wb = OHC_helper.extract_velocities_and_bolus(diagpath, expname , 
        datafilelist_uvw[tt], γ, Γ, mskC, mskW, mskS)

        w = w .* area_mask; Wb = Wb .* area_mask
        u, v = OHC_helper.UVtoTrsp(u, v, Γ)
        Ub, Vb = OHC_helper.UVtoTrsp(Ub, Vb, Γ)

        We[:, tt] .= vertical_flux_profile(w)
        Wbol[:, tt] .= vertical_flux_profile(Wb)
        Wresid[:, tt] .= vertical_flux_profile(w .+ Wb)

        _, N = OHC_helper.rotate_UV_native(u, v, cs, sn); N = N .* ϕ_min_mask
        _, Nb = OHC_helper.rotate_UV_native(Ub, Vb, cs, sn); Nb = Nb .* ϕ_min_mask

        Ne[:, tt] .= vertical_flux_profile(N)
        Nbol[:, tt] .= vertical_flux_profile(Nb)
        Nresid[:, tt] .= vertical_flux_profile(N .+ Nb)

    end
    Wdict["We"] = We
    Wdict["Wb"] = Wbol
    Wdict["Wr"] = Wresid

    Ndict["Ne"] = Ne
    Ndict["Nb"] = Nbol
    Ndict["Nr"] = Nresid
    return Wdict, Ndict
end

for expname in ["iter129_bulkformula", "iter0_bulkformula", "seasonalclimatology"]
    W, N = filter_volume_budget(diagpath, expname, γ)
    N_dict[expname] = N
    W_dict[expname] = W
end

savename = datadir("native/" * region * "_WTRSP_hovmoller.jld2")
jldsave(savename, W_dict = W_dict)

savename = datadir("native/" * region * "_NTRSP_hovmoller.jld2")
jldsave(savename, N_dict = N_dict)