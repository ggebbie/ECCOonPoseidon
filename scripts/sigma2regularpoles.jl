# transfer output from native to regularpoles grid
# 1. read monthly-average fields
# 2. interpolate to Cartesian grid
# 3. save to self-describing file.
# 4. repeat with all fields.

include("intro.jl")

using Revise
using ECCOtour, ECCOonPoseidon
using MeshArrays, MITgcmTools

include("config_exp.jl")
include("config_regularpoles.jl")

# DEFINE THE LIST OF SIGMA1 VALUES.
sig1grid = sigma1grid()

# file root names
Frootlist = ("state_3d_set1","state_3d_set2","state_2d_set1","state_2d_set2")
#Froot = "state_3d_set2"
#Froot = "state_2d_set1"
#Froot = "state_2d_set2"
###################################

# Froot= Frootlist[1] # for interactive use
for Froot in Frootlist
    pathin = sig1dir(expt)
    filelist = searchdir(pathin,"") # anything in directory
    datafilelist  = filter(x -> occursin("data",x),filelist)

    global tt = 0

    for Fname in datafilelist
        tt += 1
        println("filename ",Fname)

        year,month = timestamp_monthly_v4r4(tt)
        fileoutput = diagpath*Fname

        if month < 10
            filesuffix = "_"*string(year)*"_0"*string(month)*".nc"
        else 
            filesuffix = "_"*string(year)*"_"*string(month)*".nc"
        end

        filein = Fname[1:end-5]

        @time varsregpoles =  mdsio2regularpoles(pathin,filein,γ,nx,ny,nyarc,farc,iarc,jarc,warc,nyantarc,fantarc,iantarc,jantarc,wantarc)

        @time writeregularpoles(varsregpoles,γ,pathout,filesuffix,filelog,λC,lonatts,ϕC,latatts,z,depthatts)

    end
end
