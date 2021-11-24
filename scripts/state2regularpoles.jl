# transfer output from native to regularpoles grid
# 1. read monthly-average fields
# 2. interpolate to Cartesian grid
# 3. save to self-describing file.
# 4. repeat with all fields.

include("../src/intro.jl")

using Revise
using ECCOtour, ECCOonPoseidon
using MeshArrays, MITgcmTools

# add transport fields and accommodate for staggered grid
Frootlist = ("trsp_3d_set1","state_3d_set1","state_3d_set2","state_2d_set1","state_2d_set2")
###################################

include(srcdir("config_exp.jl"))

include(srcdir("config_regularpoles.jl"))

# Froot= Frootlist[1] # for interactive use
for Froot in Frootlist

    filelist = searchdir(diagpath,Froot)
    filelist = searchdir(diagpath,Froot) 
    datafilelist  = filter(x -> occursin("data",x),filelist)

    tt = 0

    for Fname in datafilelist
        global tt += 1
        println("filename ",Fname)

        year,month = timestamp_monthly_v4r4(tt)

        fileoutput = diagpath*Fname

        if month < 10
            filesuffix = "_"*string(year)*"_0"*string(month)*".nc"
        else 
            filesuffix = "_"*string(year)*"_"*string(month)*".nc"
        end

        filein = Fname[1:end-5]
        pathin = diagpath

        @time varsregpoles =  mdsio2regularpoles(pathin,filein,γ,nx,ny,nyarc,λarc,nyantarc,λantarc)

        @time writeregularpoles(varsregpoles,γ,pathout,filesuffix,filelog,λC,lonatts,ϕC,latatts,z,depthatts)

    end
end
