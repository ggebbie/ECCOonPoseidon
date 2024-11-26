#analysis should complete within 4 minutes 
#using 1 threads 
# julia --threads=4 --project=@. ./extract_theta_native.jl

include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, 
Printf, LinearAlgebra, LaTeXStrings, PyCall
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

region = "NPAC"; 
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
nt = length(tecco)

reg_λ = LLCcropC(λ,γ); reg_ϕ = LLCcropC(ϕ,γ)

reg_mask = LLCcropC(ocean_mask, γ)
reg_PAC = LLCcropC(PAC_msk, γ)
reg_PAC[reg_PAC .== 0.0] .= NaN
τxs = Dict(key => zeros(nt) for key in keys(shortnames)); 
τys = Dict(key => zeros(nt) for key in keys(shortnames)); 

for expname in keys(shortnames)
        filelist = searchdir(diagpath[expname],"state_2d_set1") # first filter for state_3d_set1
        τdatafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    for tt = 1:nt

        println("year ",Int(floor(tecco[tt]))," month ",((tt-1)%12)+1)
        Tname = τdatafilelist[tt]
        
        @time EXF = γ.read(diagpath[expname]*τdatafilelist[tt],MeshArray(γ,Float32,15))
        τx = EXF[:, 14]; τy = EXF[:, 15]; 
        # curlτ = curl(τx,τy,Γ)
        τE, τN = rotate_uv(τx, τy, Γ)

        τx_reg = LLCcropC(τE,γ); τy_reg = LLCcropC(τN,γ); 
        # curlτ_reg = LLCcropC(curlτ,γ);
        
        [expname][tt] = nm.mean(τx_reg .* reg_PAC)
        τys[expname][tt] = nm.mean(τy_reg .* reg_PAC)
        # curlτs .+= curlτ_reg ./ nt
    end
end

ΔF = τxs["iter129_bulkformula"] .- τxs["iter0_bulkformula"]

fig, ax = plt.subplots( figsize=(15,10))
ax.plot(tecco,  τxs["iter129_bulkformula"], label = "τx adjusted")
ax.plot(tecco,  τxs["iter0_bulkformula"], label = "τx unadjusted")
ax.set_xlabel("Years")
ax.set_ylabel("Pascals")
ax.set_title("Adjustment to Forcing")
ax.legend()
fig

fig, ax = plt.subplots( figsize=(15,10))
ax.plot(tecco,  ΔF)
ax.set_xlabel("Years")
ax.set_ylabel("Pascals")
ax.set_title("Adjustment to Forcing")
ax.legend()
fig

fig, ax = plt.subplots( figsize=(15,5))
integrate(x) = cumsum(cat([0], x[:], dims = 1))
int_ΔF= integrate(ΔF)
ax.plot(tecco, int_ΔF[1:end-1])
ax.set_xlabel("Years")
ax.set_ylabel("Pascals")
ax.set_title("Integrated Adjustment to Forcing")
ax.legend()
fig
