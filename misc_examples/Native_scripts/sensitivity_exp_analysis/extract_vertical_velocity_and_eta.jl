include("../../src/intro.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, DSP, PyCall
import PyPlot as plt 
include(srcdir("config_exp.jl"))

include(srcdir("MeshArraysPlots.jl"))

tecco = 1992+1/24:1/12:2018
nz = 50

(ϕ,λ) = latlonC(γ)
area = readarea(γ)
λ_wrap = wrap_λ(λ)
ocean_mask = wet_pts(Γ)

region = "NPAC"
PAC_msk = PAC_mask(Γ, basins, basin_list, ϕ, λ; region = region)

include(srcdir("plot_and_dir_config.jl"))
@pyimport cmocean.cm as cmo

NW_PAC = region_mask(ocean_mask, λ_wrap, ϕ, (0, 90), (0, 190))
NE_PAC = region_mask(ocean_mask, λ_wrap, ϕ, (0, 90), (190, 360))

cs, sn = get_cs_and_sn(γ)

reg_mask = LLCcropC(PAC_msk,γ)

function get_transports(diagpath::Dict{String, String}, 
    expname::String, γ::gcmgrid)

    filelist = searchdir(diagpath[expname],"trsp_3d_set1") # first filter for state_3d_set1
    datafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"

    filelist = searchdir(diagpath[expname],"state_2d_set1") # first filter for state_3d_set1
    datafilelist_η = filter(x -> occursin("data",x),filelist) # second filter for "data"

    nt = length(datafilelist); 

    ηreg = zeros(size(reg_mask)..., nt)
    ws = zeros(size(reg_mask)..., nt)
    wtop = zeros(size(reg_mask)..., nt)
    ma_template = MeshArray(γ,Float32)

    @time for tt = 1:nt
        println(tt)
        fname = datafilelist[tt]
        # get S on sigma1. Way to read a slice? (Didn't figure it out yet)
        _, _, w = extract_eulerian_velocities(diagpath, expname, fname, γ)

        wtop[:, :, tt] .= LLCcropC(w[:, 1], γ)

        ws[:, :, tt] .= LLCcropC(w[:, 43], γ)

        fname_η = datafilelist_η[tt]
        η = γ.read(diagpath[expname]*fname_η,ma_template)

        ηreg[:, :, tt] .= LLCcropC(η, γ)


    end
    var_dict = Dict()
    var_dict["Wbot"] = ws
    var_dict["Wsfc"] = wtop
    var_dict["η"] = ηreg

    return var_dict
end
        
expname = "only_kappa"
vars = get_transports(diagpath, expname, γ); 
fname = datadir("native/w_oceτ_timeseries_regular_grid_" * expname * ".jld2")
jldsave(fname, vars = vars)

expname = "iter0_bulkformula"
vars_iter0 = get_transports(diagpath, expname, γ); 
fname = datadir("native/w_oceτ_timeseries_regular_grid_" * expname * ".jld2")
jldsave(fname, vars = vars_iter0)