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

colors =  sns.color_palette("deep")[1:4]

(ϕ,λ) = latlonC(γ)
area = readarea(γ)
ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)
tecco = 1992+1/24:1/12:2018

ocean_mask = OHC_helper.wet_pts(Γ)
region = "NPAC"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; region)
cell_depths = OHC_helper.get_cell_depths(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = Float32.(OHC_helper.get_cell_volumes(area, cell_depths));

cell_depths = OHC_helper.get_cell_depths(ocean_mask, ΔzF, Γ.hFacC); 
cell_volumes_full = Float32.(OHC_helper.get_cell_volumes(area, cell_depths));
uplvl = -2e3; botlvl = -3e3; suffix = "2to3"
lvls = findall( botlvl .<= z[:].<= uplvl)

bathy_mask = Γ.Depth .> 2000; bathy_mask = Float32.(bathy_mask)
abs_dist(x, r) = abs(x) < r
ϕ_min_mask, ϕ_max_mask = OHC_helper.get_ϕ_max_min_mask(region, Γ, λ, ϕ, basins, basin_list)

cs, sn = OHC_helper.get_cs_and_sn(γ)

transports_dict = Dict()
runpath,diagpath = listexperiments(exprootdir());
mskC, mskW, mskS = OHC_helper.get_msk(Γ)

function filter_volume_budget(diagpath::Dict{String, String}, expname::String, γ::gcmgrid, 
    cs::MeshArrays.gcmarray{T,1,Matrix{T}}, sn::MeshArrays.gcmarray{T,1,Matrix{T}}) where T<:AbstractFloat
    filelist = searchdir(diagpath[expname],"trsp_3d_set1") # first filter for state_3d_set1
    datafilelist_uvw  = filter(x -> occursin("data",x),filelist) # second filter for "data"

    Wtops = Float32[]; Wbots = Float32[]; Vsouths = Float32[]; Vnorths = Float32[]
    nt = 312;
    @time for tt = 1:nt
        println(tt)
        u, v, w, Ub, Vb, Wb = OHC_helper.extract_velocities_and_bolus(diagpath, expname , 
        datafilelist_uvw[tt], γ, Γ, mskC, mskW, mskS)
        for i in eachindex(u)
            u.f[i] .= u.f[i] .+ Ub.f[i]
            v.f[i] .= v.f[i] .+ Vb.f[i]
            w.f[i] .= w.f[i] .+ Wb.f[i]
        end
  
        Utrsp, Vtrsp = OHC_helper.UVtoTrsp(u, v, Γ)
        # OHC_helper.calc_UV_conv3D!(Utrsp[:, lvls], Vtrsp[:, lvls], uv_conv)
        Etrsp, Ntrsp = OHC_helper.rotate_UV_native(Utrsp[:, lvls], Vtrsp[:, lvls], cs, sn) 
        
        Wtrsp = w .* area
        Wtop = Wtrsp[:, lvls[1]]
        Wbot = Wtrsp[:, lvls[end] + 1]

        push!(Wtops,  sum(Wtop .* PAC_msk ))
        push!(Wbots,  sum(Wbot .* PAC_msk ))
        push!(Vsouths,  sum(Ntrsp .* ϕ_min_mask )) #according to the rotations, 
        push!(Vnorths,  sum(Ntrsp .* ϕ_max_mask )) #according to the rotations, 
    end
    return (Wtop = deepcopy(Wtops), Wbot = deepcopy(Wbots), 
            Vin = deepcopy(Vsouths), Vout = deepcopy(Vnorths))
end
vars = ["seasonalclimatology", "seasonalclimatology_iter0"]
for expname in vars
    transports_dict[expname] = filter_volume_budget(diagpath, expname, γ, cs, sn)
end

svename = datadir("native/" * region * "_TRSP_" * ".jld2")
jldsave(svename, transports_dict = transports_dict)


sum(cell_volumes[:, 38:end]) / sum(cell_volumes_full[:, 38:end])
lvls