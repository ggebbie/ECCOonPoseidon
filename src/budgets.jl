function extract_ocnTAU(diagpath, expname, fname, γ)
    flux_forcing_exps = ["seasonalclimatology", "climatological_tau", "clim_tau_iter0"]
    if expname ∈ flux_forcing
        @time EXF = γ.read(diagpath[expname]*fname,MeshArray(γ,Float32,10))
        τx = EXF[:, 9]; τy = EXF[:, 10]; 
        return τx, τy
    else
        @time EXF = γ.read(diagpath[expname]*fname,MeshArray(γ,Float32,15))
        τx = EXF[:, 14]; τy = EXF[:, 15]; 
        return τx, τy
    end
end

function extract_eulerian_velocities(diagpath, expname , fnameuvw, γ)
    UV = γ.read(diagpath[expname]*fnameuvw,MeshArray(γ,Float32,150))
    u = UV[:, 1:50]; v = UV[:, 51:100]; w = UV[:, 101:150]
    return u, v, w
end


function extract_eulerian_and_bolus_velocities(
    diagpath, expname , fnameuvw, γ, Γ, 
    mskC::MeshArrays.gcmarray{T, 2, Matrix{T}}, 
    mskW::MeshArrays.gcmarray{T, 2, Matrix{T}}, 
    mskS::MeshArrays.gcmarray{T, 2, Matrix{T}})  where T<:Real
    UV = γ.read(diagpath[expname]*fnameuvw,MeshArray(γ,Float32,250))
    u = UV[:, 1:50]; v = UV[:, 51:100]; w = UV[:, 101:150]
    GM_PsiX = UV[:, 151:200]; GM_PsiY = UV[:, 201:250]
    Ub, Vb, Wb = calc_bolus(GM_PsiX,GM_PsiY, Γ, mskC, mskW, mskS)
    return u, v, w, Ub, Vb, Wb
end

function extract_lateral_heatbudget(diagpath, expname, fnameHθ, γ)
    dθλ = γ.read(diagpath[expname]*fnameHθ,MeshArray(γ,Float32,200))

    @inbounds κUθ = dθλ[:, 1:50]; 
    @inbounds κVθ = dθλ[:, 51:100]; 
    @inbounds Uθ = dθλ[:, 101:150]; 
    @inbounds Vθ = dθλ[:, 151:200]; 

    return κUθ, κVθ, Uθ, Vθ
end

function extract_verical_heatbudget(diagpath, expname, fnameRθ, γ)
    #vertical convergences
    dθr = γ.read(diagpath[expname]*fnameRθ,MeshArray(γ,Float32,150))
    @inbounds wθ = dθr[:, 1:50]; 
    @inbounds κzθ = dθr[ :, 51:100] .+ dθr[ :, 101:150];

    return κzθ, wθ
end

function calc_bolus(GM_PsiX::MeshArrays.gcmarray{T, 2, Matrix{T}},
    GM_PsiY::MeshArrays.gcmarray{T, 2, Matrix{T}}, Γ, 
    mskC::MeshArrays.gcmarray{T, 2, Matrix{T}}, 
    mskW::MeshArrays.gcmarray{T, 2, Matrix{T}}, 
    mskS::MeshArrays.gcmarray{T, 2, Matrix{T}})  where {T<:Real}
    nr=length(Γ.RC);

    bolusU = T.(similar(Γ.hFacW));
    bolusV = T.(similar(Γ.hFacS));

    @inbounds for k=1:nr-1;
        bolusU[:,k] =(GM_PsiX[:,k+1].-GM_PsiX[:,k]) ./ Γ.DRF[k];
        bolusV[:,k] =(GM_PsiY[:,k+1].-GM_PsiY[:,k]) ./ Γ.DRF[k];
    end
    @inbounds bolusU.f[:, nr] = -GM_PsiX.f[:,nr]./Γ.DRF[nr];
    @inbounds bolusV.f[:, nr] = -GM_PsiY.f[:,nr]./Γ.DRF[nr];

    #and its vertical part
    #   (seems correct, leading to 0 divergence)
    #tmpx and tmpy are the BOLUS transports
    tmp_x = T.(similar(GM_PsiX));
    tmp_y = T.(similar(GM_PsiX));
    bolusW = T.(similar(GM_PsiX));

    @inbounds for a in eachindex(tmp_x)
        tmp_x.f[a] .= GM_PsiX.f[a] .* Γ.DYG.f[a[1]]
        tmp_y.f[a] .= GM_PsiY.f[a] .* Γ.DXG.f[a[1]]
    end

    calc_UV_conv3D!(tmp_x, tmp_y, bolusW)
    @inbounds for a in eachindex(tmp_x)
        #negative for divergence instead of convergence
        #need to rescale bolus W
        bolusW.f[a] .= -1 .* bolusW.f[a] ./ Γ.RAC.f[a[1]] 
        bolusU.f[a] .= bolusU.f[a] .* mskW.f[a];
        bolusV.f[a] .= bolusV.f[a] .* mskS.f[a];
        bolusW.f[a] .= bolusW.f[a] .* mskC.f[a];
    end

    return bolusU, bolusV, bolusW
end

function exch_UV_cs3D(fldU::MeshArrays.gcmarray{T, 2, Matrix{T}},
    fldV::MeshArrays.gcmarray{T, 2, Matrix{T}}) where T<:Real
    fillval=0f0
    #step 1

    s=size.(fldU[:, 1].f)
    nz = size(fldU, 2)
    nf=fldU.grid.nFaces
    s=vcat(s,s[3]) #always has 5 faces in LLC90
    tp=fldU.grid.class

    FLDU=similar(fldU)
    FLDV=similar(fldV)
    (ovfW,ovfE,ovfS,ovfN,evfW,evfE,evfS,evfN)=MeshArrays.exch_cs_viewfunctions();
    for lvl=1:nz, a=1:nf
        FLDU.f[a, lvl] = fill(fillval,s[a][1]+1,s[a][2]);
        FLDV.f[a, lvl] = fill(fillval,s[a][1],s[a][2]+1);
        FLDU.f[a, lvl][1:s[a][1],1:s[a][2]] .= fldU.f[a, lvl];
        FLDV.f[a, lvl][1:s[a][1],1:s[a][2]] .= fldV.f[a, lvl];

        (jW, jE, jS, jN)=MeshArrays.exch_cs_target(s[a],1)
        (aW,aE,aS,aN,iW,iE,iS,iN)=MeshArrays.exch_cs_sources(a,s,1)

        if (!iseven)(a)
            (aE <= nf) && (FLDU.f[a, lvl][jE[1].-1,jE[2].-1].=ovfE(fldU.f[aE, lvl],iE[1],iE[2]))
            (aN <= nf) && (FLDV.f[a, lvl][jN[1].-1,jN[2].-1].=ovfN(fldU.f[aN, lvl],iN[1],iN[2]))
        else
            (aE <= nf) && (FLDU.f[a, lvl][jE[1].-1,jE[2].-1].=evfE(fldV.f[aE, lvl],iE[1],iE[2]))
            (aN <= nf) && (FLDV.f[a, lvl][jN[1].-1,jN[2].-1].=evfN(fldV.f[aN, lvl],iN[1],iN[2]))
        end
    end
    return FLDU,FLDV

end

function tofaces(uFLD::MeshArrays.gcmarray{T, 1, Matrix{T}},
    vFLD::MeshArrays.gcmarray{T, 1, Matrix{T}}, Γ) where T<:Real
    Ugrid = T.(similar(uFLD))
    Vgrid = T.(similar(uFLD))
    (tmpU,tmpV)=exch_UV(uFLD,vFLD)

    for a in 1:5
        (s1,s2)=size(uFLD.f[a])
        tmpU1=view(tmpU.f[a],1:s1,1:s2)
        tmpU2=view(tmpU.f[a],2:s1+1,1:s2)
        tmpV1=view(tmpV.f[a],1:s1,1:s2)
        tmpV2=view(tmpV.f[a],1:s1,2:s2+1)
        Ugrid.f[a] .=(tmpU1+tmpU2) ./ 2
        Vgrid.f[a] .=(tmpV1+tmpV2) ./ 2
    end

return Ugrid, Vgrid
end

function tofaces(uFLD::MeshArrays.gcmarray{T, 2, Matrix{T}},
    vFLD::MeshArrays.gcmarray{T, 2, Matrix{T}}, Γ) where T<:Real
    Ugrid = T.(similar(uFLD))
    Vgrid = T.(similar(uFLD))
    DXG = T.(Γ.DXG); DYG = T.(Γ.DXG); 
    (DXG,DYG)=exch_UV(DXG,DYG)

    (tmpU,tmpV)=exch_UV_cs3D(uFLD,vFLD)

    for a in eachindex(uFLD)
        (s1,s2)=size(uFLD.f[a])
        tmpU1=view(tmpU.f[a],1:s1,1:s2)
        tmpU2=view(tmpU.f[a],2:s1+1,1:s2)
        tmpV1=view(tmpV.f[a],1:s1,1:s2)
        tmpV2=view(tmpV.f[a],1:s1,2:s2+1)
        #setup weights
        tmpDXG1 =view(DXG.f[a[1]],1:s1,1:s2)
        tmpDXG2 =view(DXG.f[a[1]],2:s1+1,1:s2)
        tmpDYG1 =view(DYG.f[a[1]],1:s1,1:s2)
        tmpDYG2 =view(DYG.f[a[1]],1:s1,2:s2+1)

        Ugrid.f[a] .= ((tmpU1 .* tmpDXG1) + (tmpU2 .* tmpDXG2)) ./ (tmpDXG1 .+ tmpDXG2)
        Vgrid.f[a] .= ((tmpV1 .* tmpDYG1) + (tmpV2 .* tmpDYG2)) ./ (tmpDYG1 .+ tmpDYG2)
    end

return Ugrid, Vgrid
end