"""
    extract_meridionalΨ̄wBolus(expname, diagpath, Γ, γ, mask)


This function reads multiple data files and calculates the 
time-mean stream function by integrating the meridional transport 
(euleraina + bolus) through latitude circles.

# Arguments
- `expname`: The name of the experiment.
- `diagpath`: The path to the diagnostic files.
- `Γ`: The grid data structure containing necessary fields.
- `γ`: An object representing the data format for reading files.
- `mask`: A mask representing the region of interest.

# Returns
- `ψ̄`: The total stream function as a 2D array with NaNs 
representing non-ocean points.

"""
function extract_meridional_Ψ_Mean(expname::String,diagpath::Dict, 
    Γ::NamedTuple, γ::gcmgrid, mask::MeshArray)
    fileroot = "trsp_3d_set1"
    filelist = searchdir(diagpath[expname],fileroot) # first filter for state_3d_set1
    datafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"

    (ϕ,λ) = latlonC(γ); area = Float32.(readarea(γ))
    ϕ = Float32.(ϕ); ϕ_avg = zonal_average(ϕ, area .* mask)
    ϕ_avg = ϕ_avg[isfinite.(ϕ_avg)]
    
    LC=LatitudeCircles(ϕ_avg,Γ)

    nz=size(Γ.hFacC,2); nl=length(LC); nt = length(datafilelist)
    U=0.0*Γ.hFacW #varmeta not quite correct
    V=0.0*Γ.hFacS #varmeta not quite correct
    fill!(U, 0); fill!(V, 0)
    mskC, mskW, mskS = get_msk(Γ)

    for tt=1:nt
        println(tt)
        fname = datafilelist[tt]
        u, v, w, Ub, Vb, Wb = extract_eulerian_and_bolus_velocities(diagpath, expname, 
                                fname, γ, Γ, mskC, mskW, mskS)
        for i in eachindex(U)
            U.f[i] .+= ((u.f[i] .+ Ub.f[i])./nt)
            V.f[i] .+= ((v.f[i] .+ Vb.f[i])./nt)
        end
    end

    (Utr,Vtr)=UVtoTransport(U,V,Γ)
    Utr = Utr .* mask; Vtr = Vtr .* mask;
    ov=Array{Float64,2}(undef,nl,nz)
    for z=1:nz
        [ov[l,z]=ThroughFlow(Utr[:,z], Vtr[:,z],LC[l]) for l=1:nl]
    end
    ov= reverse(cumsum(reverse(ov,dims=2),dims=2),dims=2)
    ψ̄ = -reverse(ov,dims=2); ψ̄[ψ̄.==0.0].=NaN
    ψ̄ = reverse(ψ̄', dims = 1)

    return ψ̄, ϕ_avg
end


"""
    extract_meridional_Ψ(expname, diagpath, Γ, γ, mask)

This function reads multiple data files and calculates the 
time-series of the stream function by integrating the meridional transport 
(eulerian and bolus) through latitude circles. This is the streamfunction that 
advects tracers.

# Arguments
- `expname`: The name of the experiment.
- `diagpath`: The path to the diagnostic files.
- `Γ`: The grid data structure containing necessary fields.
- `γ`: An object representing the data format for reading files.
- `mask`: A mask representing the region of interest.

# Returns
- `ψ̄`: The total stream function as a 3D (z, x, t) array with NaNs 
representing non-ocean points.

"""
function extract_meridional_Ψ(expname::String,diagpath::Dict, 
    Γ::NamedTuple, γ::gcmgrid, mask::MeshArray)
    fileroot = "trsp_3d_set1"
    filelist = searchdir(diagpath[expname],fileroot) # first filter for state_3d_set1
    datafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    
    (ϕ,λ) = latlonC(γ); area = Float32.(readarea(γ))
    ϕ = Float32.(ϕ); ϕ_avg = zonal_average(ϕ, area .* mask)
    ϕ_avg = ϕ_avg[isfinite.(ϕ_avg)]
    
    LC=LatitudeCircles(ϕ_avg,Γ)
    nz=size(Γ.hFacC,2); nl=length(LC); nt = length(datafilelist)
    ψ = zeros(nl,nz, nt)
    mskC, mskW, mskS = get_msk(Γ)
    #use this for more efficiency
    for tt=1:nt
        println(tt)
        fname = datafilelist[tt]
        u, v, w, Ub, Vb, Wb = extract_eulerian_and_bolus_velocities(diagpath, expname, 
                                fname, γ, Γ, mskC, mskW, mskS)
        u .+= Ub
        v .+= Vb

        u = u .* mask; v = v .* mask; 
        (Utr,Vtr)=UVtoTrsp(u,v,Γ, γ)
        ov=Array{Float64,2}(undef,nl,nz)
        for z=1:nz
            [ov[l,z]=ThroughFlow(Utr[:,z], Vtr[:,z],LC[l]) for l=1:nl]
        end
        ov= reverse(cumsum(reverse(ov,dims=2),dims=2),dims=2)
        ψ̄ = -reverse(ov,dims=2); ψ̄[ψ̄.==0.0].=NaN
        ψ[:, :, tt] .= ψ̄
    end

    ψ = reverse(permutedims(ψ, (2, 1, 3)), dims = 1)

    return ψ, ϕ_avg

end

#inverse of TabtoMsk
function LatitudeCirclesMasks(LC, γ)
    #low memory footprint MeshArray(γ, Int8, 180)
    nl=length(LC);

    Umask = MeshArray(γ, Int8, nl)
    Vmask = MeshArray(γ, Int8, nl)

    fill!(Umask, 0.0)
    fill!(Vmask, 0.0)

    for l = 1:nl
        IntegralPath = LC[l]
        tabW=IntegralPath.tabW
        tabS=IntegralPath.tabS
        nW = size(tabW,1); nS = size(tabS,1)
        for k=1:nW
            (a,i1,i2,w)=tabW[k,:]
            Umask.f[a, l][i1, i2] = w
        end

        for k=1:nS
            (a,i1,i2,w)=tabS[k,:]
            Vmask.f[a, l][i1, i2] = w
        end
    end

return Umask, Vmask

end
"""
    ThroughFlow(Utrsp::MeshArray, Vtrsp::MeshArray, IntegralPath)

Compute transport through an integration path.

This function calculates the transport through an integration path defined by latitude circles
using the provided velocity components Utrsp and Vtrsp.

# Arguments
- `Utrsp::MeshArray`: Velocity component in the eastward direction (u) as a MeshArray.
- `Vtrsp::MeshArray`: Velocity component in the northward direction (v) as a MeshArray.
- `IntegralPath`: An object representing the integration path, which contains required data.

# Returns
- `Float64`: The computed transport value through the integration path.

# Examples
```julia
result = ThroughFlow(Utrsp, Vtrsp, IntegralPath)
println("Transport value: ", result)
"""
function ThroughFlow(Utrsp::MeshArrays.gcmarray{T, 1, Matrix{T}}, 
    Vtrsp::MeshArrays.gcmarray{T, 1, Matrix{T}},IntegralPath) where {T<:Real}

    trsp=Array{Float32}(undef,1)

    #these will return the indexs and coefficients required for Transport Through a particular Latitude Circle
    
    tabW=IntegralPath.tabW
    tabS=IntegralPath.tabS
  
    nW = size(tabW,1); nS = size(tabS,1)
    trsp[1]=0.0
    for k=1:nW
        (a,i1,i2,w)=tabW[k,:]
        u=Utrsp.f[a][i1,i2]
        trsp[1]=trsp[1]+w*u
    end
    for k=1:nS
        (a,i1,i2,w)=tabS[k,:]
        v=Vtrsp.f[a][i1,i2]
        trsp[1]=trsp[1]+w*v
    end
  
    return trsp[1]
end


function UVtoTrsp(U::MeshArrays.gcmarray{T, 2, Matrix{T}},V::MeshArrays.gcmarray{T, 2, Matrix{T}},
    G::NamedTuple, γ::gcmgrid) where {T<:Real}
    nz = size(U, 2)
    Utrsp = MeshArray(γ, T, nz)
    Vtrsp = MeshArray(γ, T, nz)

    for a in eachindex(U)
        Utrsp.f[a] .= G.DRF[a[2]]*U.f[a].*G.DYG.f[a[1]]
        Vtrsp.f[a] .= G.DRF[a[2]]*V.f[a].*G.DXG.f[a[1]]
    end
    return Utrsp, Vtrsp
end
