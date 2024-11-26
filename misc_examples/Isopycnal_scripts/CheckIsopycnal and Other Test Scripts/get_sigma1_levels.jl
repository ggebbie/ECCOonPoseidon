println(expname)
filelist = searchdir(diagpath[expname], "state_2d_set1") # first filter for state_3d_set1
datafilelist_S  = filter(x -> occursin("data",x),filelist) # second filter for "data"
filelist = searchdir(diagpath[expname],"state_3d_set1") # first filter for state_3d_set1
datafilelist_θ  = filter(x -> occursin("data",x),filelist) # second filter for "data"

fnameS = datafilelist_S[1]
fnameθS = datafilelist_θ[1] #potential temp and abs. salinity

θS = γ.read(diagpath[expname]*fnameθS,MeshArray(γ,Float32,100))
ETAN = γ.read(diagpath[expname]*fnameS,MeshArray(γ,Float32,1))[:, 1]

θ = θS[:, 1:50]
S = θS[:, 51:end]

lon=0.;lat=30.

SA = similar(sθ); Θ = similar(sθ);σ = similar(sθ)
p₀ = 1000

for ff = 1:5, k = 1:50
    SA.f[ff, k] .= gsw_sa_from_sp.(S.f[ff, k],-pstdz[k],lon,lat)
    Θ.f[ff, k] .= gsw_ct_from_pt.(SA.f[ff, k],θ.f[ff, k])
    σ.f[ff, k] .= gsw_rho.(SA.f[ff, k],Θ.f[ff, k],p₀) .- 1000.
end

σavg = volume_average_by_depth(σ, cell_volumes, γ)

(z[37], z[43])
σtop, σbot = (σavg[37], σavg[43])



