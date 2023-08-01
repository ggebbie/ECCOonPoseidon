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
function ThroughFlow(Utrsp::MeshArray, Vtrsp::MeshArray,IntegralPath)

    trsp=Array{Float64}(undef,1)

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
function extract_meridional_Ψ_Mean_EulBolus(expname,diagpath, Γ, γ, mask)
    fileroot = "trsp_3d_set1"
    filelist = searchdir(diagpath[expname],fileroot) # first filter for state_3d_set1
    datafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    #[1:36, 312-36+1:312] comparing first and last years
    tecco = 1992+1/24:1/12:2018
    LC=LatitudeCircles(-89.0:89.0,Γ)
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
    return ψ̄
end

"""
    extract_meridionalΨ̄(expname, diagpath, Γ, γ, mask)


This function reads multiple data files and calculates the 
time-mean stream function by integrating the meridional transport 
(eulerian only) through latitude circles.

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

function extract_meridionalΨ̄(expname,diagpath, Γ, γ, mask)
    fileroot = "trsp_3d_set1"
    filelist = searchdir(diagpath[expname],fileroot) # first filter for state_3d_set1
    datafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    #[1:36, 312-36+1:312] comparing first and last years
    LC=LatitudeCircles(-89.0:89.0,Γ)
    nz=size(Γ.hFacC,2); nl=length(LC); nt = length(datafilelist)
    U=0.0*Γ.hFacW #varmeta not quite correct
    V=0.0*Γ.hFacS #varmeta not quite correct
    fill!(U, 0); fill!(V, 0)
    for tt=1:nt
        fname = datafilelist[tt]
        UV = γ.read(diagpath[expname]*fname,MeshArray(γ,Float32,100))
        u = UV[:, 1:50]
        v = UV[:, 51:100]
        for i in eachindex(U)
            U[i]=U[i]+(u[i]/nt)
            V[i]=V[i]+(v[i]/nt)
        end
    end

    (Utr,Vtr)=UVtoTransport(U,V,Γ)
    Utr = Utr .* mask; Vtr = Vtr .* mask;
    ov=Array{Float64,2}(undef,nl,nz)
    for z=1:nz
        UV=Dict("U"=>Utr[:,z],"V"=>Vtr[:,z],"dimensions"=>["x","y"])
        [ov[l,z]=ThroughFlow(UV,LC[l],Γ) for l=1:nl]
    end
    ov= reverse(cumsum(reverse(ov,dims=2),dims=2),dims=2)
    ψ̄ = -reverse(ov,dims=2); ψ̄[ψ̄.==0.0].=NaN
    return ψ̄
end

"""
    extract_meridional_ΨwBolus(expname, diagpath, Γ, γ, mask)


This function reads multiple data files and calculates the 
time-series of the stream function by integrating the meridional transport 
(eulerian and bolus) through latitude circles.

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
function extract_meridional_Ψ_EulBolus(expname,diagpath, Γ, γ, mask)
    fileroot = "trsp_3d_set1"
    filelist = searchdir(diagpath[expname],fileroot) # first filter for state_3d_set1
    datafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    #[1:36, 312-36+1:312] comparing first and last years
    LC=LatitudeCircles(-89.0:89.0,Γ)
    nz=size(Γ.hFacC,2); nl=length(LC); nt = length(datafilelist)

    ψ = zeros(nl,nz, nt)
    mskC, mskW, mskS = get_msk(Γ)

    for tt=1:nt
        println(tt)
        fname = datafilelist[tt]
        u, v, w, Ub, Vb, Wb = extract_eulerian_and_bolus_velocities(diagpath, expname, 
                                fname, γ, Γ, mskC, mskW, mskS)
        u .+= Ub
        v .+= Vb
        (Utr,Vtr)=UVtoTransport(u,v,Γ)
        Utr = Utr .* mask; Vtr = Vtr .* mask;
        ov=Array{Float64,2}(undef,nl,nz)
        for z=1:nz
            [ov[l,z]=ThroughFlow(Utr[:,z],Vtr[:,z],LC[l]) for l=1:nl]
        end
        ov= reverse(cumsum(reverse(ov,dims=2),dims=2),dims=2)
        ψ̄ = -reverse(ov,dims=2); ψ̄[ψ̄.==0.0].=NaN
        ψ[:, :, tt] .= ψ̄
    end

    return ψ

end

"""
extract_meridionalΨ̄timeseries(expname, diagpath, Γ, γ, mask)


This function reads multiple data files and calculates the 
time-series of the stream function by integrating the meridional transport 
(eulerian only) through latitude circles.

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
function extract_meridionalΨ̄timeseries(expname,diagpath, Γ, γ, mask)
    fileroot = "trsp_3d_set1"
    filelist = searchdir(diagpath[expname],fileroot) # first filter for state_3d_set1
    datafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    #[1:36, 312-36+1:312] comparing first and last years
    LC=LatitudeCircles(-89.0:89.0,Γ)
    nz=size(Γ.hFacC,2); nl=length(LC); nt = length(datafilelist)
    U=0.0*Γ.hFacW #varmeta not quite correct
    V=0.0*Γ.hFacS #varmeta not quite correct
    fill!(U, 0); fill!(V, 0)
    ψ = zeros(nl,nz, nt)
    for tt=1:nt
        println(tt)
        fname = datafilelist[tt]
        UV = γ.read(diagpath[expname]*fname,MeshArray(γ,Float32,100))
        u = UV[:, 1:50]
        v = UV[:, 51:100]
        (Utr,Vtr)=UVtoTransport(u,v,Γ)
        Utr = Utr .* mask; Vtr = Vtr .* mask;
        ov=Array{Float64,2}(undef,nl,nz)
        for z=1:nz
            UV=Dict("U"=>Utr[:,z],"V"=>Vtr[:,z],"dimensions"=>["x","y"])
            [ov[l,z]=ThroughFlow(UV,LC[l],Γ) for l=1:nl]
        end
        ov= reverse(cumsum(reverse(ov,dims=2),dims=2),dims=2)
        ψ̄ = -reverse(ov,dims=2); ψ̄[ψ̄.==0.0].=NaN
        ψ[:, :, tt] .= ψ̄
    end

    return ψ
end
