include("../../src/intro.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, PyCall
import PyPlot as plt 

@pyimport seaborn as sns
sns.set_theme(context = "talk", style = "ticks");

include(srcdir("config_exp.jl"))

tecco = 1992+1/24:1/12:2018
nz = 50

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

include("sections_helper.jl")


runpath,diagpath = listexperiments(exprootdir());
 
ocean_mask = wet_pts(Γ)
abs_dist(x, r) = abs(x) < r

region2 = "PAC56"
PAC56_mask = PAC_mask(Γ, basins, basin_list, ϕ, λ; region = region2, extent = "full")
sections = Dict()

sections["P01"] = latitude_section_mask(ϕ, 45.6, PAC56_mask)
sections["P02"] = latitude_section_mask(ϕ, 30.3, PAC56_mask)
sections["P03"] = latitude_section_mask(ϕ, 25.3, PAC56_mask)
sections["P06"] = latitude_section_mask(ϕ, -31.9, PAC56_mask)
sections["P21"] = latitude_section_mask(ϕ, -17.8, PAC56_mask)

sections["P09"] = longitude_section_mask(λ, 138.2, PAC56_mask)
sections["P10"] = longitude_section_mask(λ, 147.4, PAC56_mask)
sections["P14"] = longitude_section_mask(λ, 178.3, PAC56_mask)
sections["P16"] = longitude_section_mask(λ, -151.1, PAC56_mask)
sections["P17"] = longitude_section_mask(λ, -139.9, PAC56_mask)
sections["P18"] = longitude_section_mask(λ, -106.2, PAC56_mask)

function filter_heat_budget(diagpath::Dict{String, String}, 
    expname::String, γ::gcmgrid)
    nz = 50

    t_rand = zeros(3)

    t_rand[1] = rand(tecco[findall( 1993 .<= collect(tecco).<= 2002)])
    t_rand[2:3] =  t_rand[1] .+ [5 + rand(0:0.1:1), 10 + rand(0:0.1:1)]
    [t_rand[i] = Base.findmin(abs.(collect(tecco) .- t_rand[i]))[2] for i = 1:3]
    t_rand = Int.(t_rand)

    tecco_sample = tecco[t_rand]
    E,F_samp = trend_matrices(tecco_sample)
    println(F_samp)

    filelist = searchdir(diagpath[expname],"state_3d_set1") 
    datafilelist_θ  = filter(x -> occursin("data",x),filelist)[t_rand]

    nt = length(datafilelist_θ); 
    θ_trends = MeshArray(γ,Float32,nz); fill!(θ_trends, 0.0)

    @time for tt = 1:nt
        println(tt)
        θname = datafilelist_θ[tt]
        # get S on sigma1. Way to read a slice? (Didn't figure it out yet)
        θ = γ.read(diagpath[expname]*θname,MeshArray(γ,Float32,nz))
        for a in eachindex(θ)
            θ_trends.f[a] .+= F_samp[2, tt] .* θ.f[a]
        end
    end

    return θ_trends
end


WOCE_mask = MeshArray(γ,Float32); fill!(WOCE_mask, 0.0)
[WOCE_mask .+= mask for mask in values(sections)]
WOCE_mask = Γ.hFacC .* WOCE_mask
for a in eachindex(WOCE_mask)
    WOCE_mask.f[a][WOCE_mask.f[a] .> 0] .= 1
end

proj0 = ECCOonPoseidon.cartopy.crs.Robinson(central_longitude=-150)
projPC = ECCOonPoseidon.cartopy.crs.PlateCarree()

fig, ax = plt.subplots( figsize=(15,10), subplot_kw=Dict("projection"=> proj0))
ax.coastlines(resolution="110m")
gl = ax.gridlines(crs=projPC, draw_labels=true,
                    linewidth=2, color="gray", alpha=0, linestyle="--")
gl.top_labels = false; gl.right_labels = false
plot_basin_mask!(λ, ϕ, WOCE_mask, projPC)    
ax.set_extent((120, 295, -70, 70),crs=projPC)
fig.savefig(plotsdir("native/sections/WOCE_Sections.png"), bbox_inches = "tight")

fig

expname = "iter129_bulkformula"

n_iterations = Int(100)
trends = zeros(n_iterations)
θ_trends_ens = zeros(Float32, 50, n_iterations)
θ_std_ens = zeros(Float32, 50, n_iterations)

for iN in 1:n_iterations
    θ_trends = filter_heat_budget(diagpath, expname, γ)

    mean_trends = lateral_sum(θ_trends .* WOCE_mask) ./ lateral_sum(WOCE_mask)

    std_trends = MeshArray(γ,Float32, 50)
    for a in eachindex(std_trends)
        std_trends.f[a] .= (θ_trends.f[a] .- mean_trends[a[2]]).^2
    end

    std_trends = lateral_sum(std_trends .* WOCE_mask) ./ lateral_sum(WOCE_mask)

    θ_trends_ens[:, iN] .= mean_trends
    θ_std_ens[:, iN] .= sqrt.(std_trends)

end

mean_mean_trends = mean(θ_trends_ens, dims = 2)[:]

mean_std_trends = mean(θ_std_ens, dims = 2)[:]
std_std_trends = std(θ_std_ens, dims = 2)[:]

μ = mean_mean_trends .* 1000
σ = (mean_std_trends .+ std_std_trends) .* 1000 

fig, ax = plt.subplots( figsize=(5,10))
ax.plot(μ, z, color = "k", linestyle = "dashed", label = "WOCE Sampling")
ax.fill_betweenx(z, μ .- (2 .*σ), μ .+ (2 .*σ), alpha=0.3, color="k", label = "2σ uncertainty estimate")
ax.set_ylim(2000, 5000); ax.set_xlim(-10, 10)
ax.set_title("Pacific Ocean Temperature Trends")
ax.invert_yaxis()
ax.legend()
ax.set_xlabel("m°C per year")
fig.savefig(plotsdir("native/sections/WOCE_Trends_Trend_Dist.png"), bbox_inches = "tight")
fig