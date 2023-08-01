println("Configure the experiment variables and directories")

# Check that the output has been copied to batou from poseidon
!isdir(exprootdir(expt)) && error("Experiment missing")

diagpath = diagdir(expt)

vastexprootdir() = "/vast/ECCOv4r4/exps"

update_paths_dict!(main_dict, add_dict) = [main_dict[key] = add_dict[key] for key ∈ keys(add_dict) if key ∉ keys(main_dict)]

# print output here
path_out = sig1dir(expt)
!isdir(path_out) ? mkdir(path_out) : nothing;

pathout = regpolesdir(expt)
!isdir(pathout) ? mkdir(pathout) : nothing;

################################################################
# get MITgcm / ECCOv4r4 LLC grid and depth information. Store in γ.
pth = MeshArrays.GRID_LLC90
γ = GridSpec("LatLonCap",pth)
Γ = GridLoad(γ;option="full")

#fix floating point errors for Float32 vars
check_zero(x) = abs(x) < 1e-4 #check if within float32 zero 
check_notzero(x) = abs(x) > 1e-7 #check if within float32 zero 
withineps(x) = convert(Float32,x)* check_notzero(x) #replace if close to zero 

for ff = 1:5
    Γ.AngleCS.f[ff] = withineps.(Γ.AngleCS.f[ff])
    Γ.AngleSN.f[ff] = withineps.(Γ.AngleSN.f[ff])
    isoneSN = check_zero.(1 .- Γ.AngleSN.f[ff])
    isoneCS = check_zero.(1 .- Γ.AngleCS.f[ff])
    
    Γ.AngleCS.f[ff][isoneSN] .= sign.(Γ.AngleCS.f[ff][isoneSN])
    Γ.AngleSN.f[ff][isoneCS] .= sign.(Γ.AngleSN.f[ff][isoneCS])

    Γ.AngleCS.f[ff] .= round.(Γ.AngleCS.f[ff], digits = 4)
    Γ.AngleSN.f[ff] .= round.(Γ.AngleSN.f[ff], digits = 4)
end

# no longer needed?
#γ = setupLLCgrid(datadir("grid/"))

nf = length(γ.fSize)

# get standard levels of MITgcm
z = depthlevels(γ)
z = -vec(z)

pstdz = pressurelevels(z)
p₀ = 1000.0; # dbar

Δz = read_mdsio(γ.path,"DRC")
ΔzF = vec(read_mdsio(γ.path,"DRF"))
Δz = vec(Δz)

#retrieves the array containing all basin IDs 
basins=read(joinpath(pth,"v4_basin.bin"),MeshArray(γ,Float32))

basin_list=["Pacific","Atlantic","indian","Arctic","Bering Sea",
            "South China Sea","Gulf of Mexico","Okhotsk Sea",
            "Hudson Bay","Mediterranean Sea","Java Sea","North Sea",
            "Japan Sea", "Timor Sea","East China Sea","Red Sea",
            "Gulf","Baffin Bay","GIN Seas","Barents Sea"];


            