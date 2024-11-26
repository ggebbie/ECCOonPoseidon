#analysis should complete within 4 minutes 
#using 1 threads 
# julia --threads=4 --project=@. ./extract_theta_native.jl

include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, 
Printf, LinearAlgebra, LaTeXStrings, PyCall, GibbsSeaWater
import PyPlot as plt
import NaNMath as nm
@pyimport matplotlib.animation as anim

cm = pyimport("cmocean.cm");
@pyimport seaborn as sns;

sns.set_theme(context = "poster", style = "ticks", font_scale = 1.0,
              palette = sns.color_palette("colorblind"));
              
include(srcdir("config_exp.jl"))

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

region = "PAC56"; 
ocean_mask = OHC_helper.wet_pts(Γ)
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not")
sns.set_theme(context = "poster", style = "white", font_scale = 1.0,
              palette = sns.color_palette("deep"));

#computed on standard pressures, 
#may have a different value if using pressures diagnosed from ECCO
#https://ecco-v4-python-tutorial.readthedocs.io/VectorCalculus_ECCO_barotropicVorticity.html
runpath,diagpath = listexperiments(exprootdir());
ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)

tecco= 1992+1/24:1/12:2018 # ecco years
nt = length(tecco); nz = length(z)

reg_λ = LLCcropC(λ,γ); reg_ϕ = LLCcropC(ϕ,γ)
reg_mask = LLCcropC(ocean_mask, γ)
reg_PAC = LLCcropC(PAC_msk, γ)
dxG = Γ.DXG .* Γ.hFacS[:, 1]; dyG  = Γ.DYG .* Γ.hFacW[:, 1]; 

V_ek_PAC = Dict()
mask_array!(x, mask) = (x[mask .== 0] .= NaN)
similar_zeros(x) = zeros(size(x)); 

for expname in keys(shortnames)
    filelist = searchdir(diagpath[expname],"state_2d_set1") # first filter for state_3d_set1
    datafilelist_τ  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    filelist = searchdir(diagpath[expname],"state_3d_set1") # first filter for state_3d_set1
    datafilelist_θ  = filter(x -> occursin("data",x),filelist) # second filter for "data"

    τxs = similar_zeros(reg_ϕ); 
    τys = similar_zeros(reg_ϕ); 
    σ0  =similar_zeros(reg_ϕ); 
    for tt = 1:nt
        println("year ",Int(floor(tecco[tt]))," month ",((tt-1)%12)+1)
        τname = datafilelist_τ[tt]
        θname = datafilelist_θ[tt]

        @time EXF = γ.read(diagpath[expname]*τname,MeshArray(γ,Float32,15))
        @time θSz = γ.read(diagpath[expname]*θname,MeshArray(γ,Float32,2*nz))
        σz = similar_zeros(reg_ϕ);
        ekz = 20 
        for i = 1:ekz #depth avg T & S first 300 meters
            θz = LLCcropC(θSz[:, 1+i],γ) 
            Sz = LLCcropC(θSz[:, nz+i],γ)

            SA = gsw_sa_from_sp.(Sz,0,0,30)
            CT = gsw_ct_from_pt.(SA,θz)
            σ0 .+= gsw_rho.(SA,CT,0) ./ ekz #reference at the surface
        end
        τx = EXF[:, 14] .* dxG; τy = EXF[:, 15] .* dyG; 
        τE, τN = rotate_uv(τx, τy, Γ)
        τx_reg = LLCcropC(τE,γ); τy_reg = LLCcropC(τN,γ); 
        τxs .+= τx_reg ./ nt
        τys .+= τy_reg ./ nt
        σ0 .+= σz ./ nt
    end

    σ0_Inf = copy(σ0); σ0_Inf[σ0_Inf.== 0] .= Inf
    Ω = 2π/86400
    f = 2Ω .* sind.(reg_ϕ)
    V_ek = -τxs ./ (f .* σ0_Inf)
    V_ek[abs.(reg_ϕ) .< 2] .= 0.0
    V_ek_PAC[expname] = 1e-6 .*  sum(V_ek.* reg_PAC, dims = 1)[:]
end

fig, ax = plt.subplots( figsize=(15,7.5))
ax.set_title("Integrated Surface Ekman Transport in Pacific Ocean")
ax.set_xlabel("Longitude")
ax.set_ylabel("[Sv]")
for expname in keys(shortnames)
    ax.plot(vec(reg_ϕ[1, :]), V_ek_PAC[expname], label = expname)
end
ax.set_xlim(20, 55)
ax.set_ylim(-7.5, 2)
ax.legend()
fig.tight_layout()
fig
