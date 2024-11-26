include("../../src/intro.jl")

using Revise, ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, DrWatson, 
     LaTeXStrings, PyCall
import NaNMath as nm
import PyPlot as plt

include(srcdir("config_exp.jl"))

cmo = pyimport("cmocean.cm");
@pyimport seaborn as sns;

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

tecco = 1992+1/24:1/12:2018

ocean_mask = wet_pts(Γ)
runpath,diagpath = listexperiments(exprootdir());


function filter_volume_budget(diagpath::Dict{String, String}, expname::String, γ::gcmgrid, 
    Γ, mask,
    ϕ_min_mask, ϕ_max_mask, lvls)
    filelist = searchdir(diagpath[expname],"trsp_3d_set1") # first filter for state_3d_set1
    datafilelist_uvw  = filter(x -> occursin("data",x),filelist) # second filter for "data"

    cs, sn = get_cs_and_sn(γ)
    mskC, mskW, mskS = get_msk(Γ)

    area = readarea(γ) .* mask

    Wtops = Float32[]; Wbots = Float32[]; Vsouths = Float32[]; Vnorths = Float32[]
    nt = 312;
    @time for tt = 1:nt
        println(tt)
        u, v, w, Ub, Vb, Wb = extract_eulerian_and_bolus_velocities(diagpath, expname , 
        datafilelist_uvw[tt], γ, Γ, mskC, mskW, mskS)

        u .+= Ub
        v .+= Vb
        w .+= Wb
  
        Utrsp, Vtrsp = UVtoTrsp(u, v, Γ, γ)
        Etrsp, Ntrsp = rotate_UV_native(Utrsp[:, lvls], Vtrsp[:, lvls], cs, sn) 
        
        Wtrsp = w .* area
        Wtop = Wtrsp[:, lvls[1]]
        Wbot = Wtrsp[:, lvls[end] + 1]

        push!(Wtops,  sum(Wtop))
        push!(Wbots,  sum(Wbot))
        push!(Vsouths,  sum(Ntrsp .* ϕ_min_mask )) #according to the rotations, 
        push!(Vnorths,  sum(Ntrsp .* ϕ_max_mask )) #according to the rotations, 
    end
    return (Wtop = deepcopy(Wtops), Wbot = deepcopy(Wbots), 
            Vin = deepcopy(Vsouths), Vout = deepcopy(Vnorths))
end

vars = ["iter0_bulkformula", "iter129_bulkformula"]
for expname in vars
    transports = filter_volume_budget(diagpath, expname, γ, Γ, 
    PAC_msk, ϕ_min_mask, ϕ_max_mask, lvls)
    svename = datadir("native/" * region * "_" * expname * "_TRSP" * ".jld2")
    jldsave(svename, transports = transports)
end





sum(cell_volumes[:, 38:end]) / sum(cell_volumes_full[:, 38:end])
lvls