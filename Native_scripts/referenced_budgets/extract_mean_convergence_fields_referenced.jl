#analysis should complete within 12-13 minutes 
#using 12 threads 
# julia --threads=6 --project=@. ./extract_heat_budget_native.jl

include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, 
    BenchmarkTools, LaTeXStrings, PyCall
using .OHC_helper
import NaNMath as nm
include(srcdir("config_exp.jl"))

cmo = pyimport("cmocean.cm");
@pyimport seaborn as sns;
sns.set_theme(context = "notebook", style = "ticks", font_scale = 1.0,
              palette = sns.color_palette("colorblind"));
(ϕ,λ) = latlonC(γ)
area = readarea(γ)
ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)

ocean_mask = OHC_helper.wet_pts(Γ)
cell_depths = OHC_helper.get_cell_depths(ocean_mask, ΔzF, Γ.hFacC); 

H = OHC_helper.smush(cell_depths, γ); H[findall(H .== 0.0)] = Inf
inv_H = 1 ./ H
runpath,diagpath = listexperiments(exprootdir());
for expname in keys(shortnames)
    filelist = searchdir(diagpath[expname], "state_2d_set1") # first filter for state_3d_set1
    datafilelist_S  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    filelist = searchdir(diagpath[expname],"state_3d_set1") # first filter for state_3d_set1
    datafilelist_θ  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    filelist = searchdir(diagpath[expname],"trsp_3d_set1") # first filter for state_3d_set1
    datafilelist_uvw  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    filelist = searchdir(diagpath[expname],"trsp_3d_set2")
    datafilelist_H  = filter(x -> occursin("data",x),filelist) 
    filelist = searchdir(diagpath[expname],"trsp_3d_set3")
    datafilelist_R  = filter(x -> occursin("data",x),filelist)

    sθbar = MeshArray(γ,Float32,50); fill!(sθbar, 0.0)
    Ubar = MeshArray(γ,Float32,50); fill!(Ubar, 0.0)
    Vbar = MeshArray(γ,Float32,50); fill!(Vbar, 0.0)
    Wbar = MeshArray(γ,Float32,50); fill!(Wbar, 0.0)
    nt = 312
    for tt = 1:nt
        println(tt)
        sθ = OHC_helper.extract_sθ(expname,diagpath, γ, datafilelist_S[tt], datafilelist_θ[tt], inv_H)
        u, v, w = OHC_helper.extract_velocities(diagpath, expname, datafilelist_uvw[tt], γ)

        #take the time average of velocity and referenced faces
        for a in eachindex(sθbar)
            sθbar.f[a] .+= (sθ.f[a]) ./ nt
            Ubar.f[a] .+= (u.f[a]) ./ nt
            Vbar.f[a] .+= (v.f[a]) ./ nt
            Wbar.f[a] .+= (w.f[a]) ./ nt
        end
    end

    #interpolate θ to the faces after time average, since these averages are interchangeable 
    θUbar, θVbar, θWbar = MeshArray(γ,Float32,50), MeshArray(γ,Float32,50), MeshArray(γ,Float32,50)
    for k = 1:49 #dont do this for the bottom level
        tmp = sθbar[:, k]
        tmp1, tmp2 = OHC_helper.tofaces(tmp, tmp)
        θUbar.f[:, k] .= tmp1.f; θVbar.f[:, k] .= tmp2.f
    end

    dθref = Dict("sθbar" => sθbar, "sθUbar" => θUbar, 
                "sθVbar" => θVbar, "sθWbar" => θWbar,
                "ubar" => Ubar, "vbar" => Vbar, "wbar" => Wbar)

    savename = datadir("native/" * expname * "_THETA_budget_meanUθ" * ".jld2")
    println("saving file " * savename)
    jldsave(savename, dθref = dθref)
end 