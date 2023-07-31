"""
    ThroughFlow(Utrsp::MeshArray, Vtrsp::MeshArray,IntegralPath)

Compute transport through an integration path
"""
function ThroughFlow(Utrsp::MeshArray, Vtrsp::MeshArray,IntegralPath)


    trsp=Array{Float64}(undef,1)

    #these will return the indexs and coefficients required for Transport Through a particular Latitude Circle
    
    tabW=IntegralPath.tabW
    tabS=IntegralPath.tabS
  
    trsp[1]=0.0
    for k=1:size(tabW,1)
        (a,i1,i2,w)=tabW[k,:]
        u=Utrsp.f[a][i1,i2]
        trsp[1]=trsp[1,i3,i4]+w*u
    end
    for k=1:size(tabS,1)
        (a,i1,i2,w)=tabS[k,:]

        v=Vtrsp.f[a][i1,i2]
        trsp[1]=trsp[1,i3,i4]+w*v
    end
  

    return trsp[1]
end