#analysis should complete within 12 minutes 
#using 12 threads 
include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise
using ECCOonPoseidon, ECCOtour,
    MeshArrays, MITgcmTools, JLD2, 
    DrWatson, BenchmarkTools, LaTeXStrings
using .OHC_helper
using Plots
using ColorSchemes
import NaNMath as nm
include(srcdir("config_exp.jl"))

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

runpath,diagpath = listexperiments(exprootdir());
ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)
marks = expsymbols()
 
ocean_mask = wet_pts(Γ)
region = "PAC"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "full")
msk = PAC_msk;
(msk == ocean_mask) && (region = "GLOB")
tecco = 1992+1/24:1/12:2018

Ψ_exp = Dict()
Ψ_exp_mean = Dict()
X=(collect(-89.0:89.0)); Y=reverse(z); #coordinate variables
Hinv = smush(get_cell_depths(msk, ΔzF, Γ.hFacC)); Hinv[findall(Hinv .== 0)] = Inf
Hinv = 1 ./ Hinv; Hinv = Float32.(Hinv)

f0(ϕ) = abs(ϕ) > 2 ? (2 * (2π / 86400) * sind(ϕ)) : Inf32
fs = MeshArray(γ,Float32,1)[:, 1]; fs.f .= map.(x -> f0.(x), ϕ.f)
finv = 1 ./ fs; finv = Float32.(finv) 

@time for (key,values) in shortnames
    expname = key
    println(expname)
    Ψ_exp_mean[expname] = extract_meridionalΨEK(expname,diagpath, Γ, γ,
                                           Hinv, finv, msk)
end

Ψ_sym = L"\Psi^{Ek}"

