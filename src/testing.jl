using ECCOtour
"""
    matrixfilter(F,froot,years,γ)
    writing it in a funny way to save computation
    issue with timeseries being read in different files
# Arguments
- `w`: weight of timemean
- `froot`: filename root
- `years`: iterator for multiple files
- `γ`: GCM grid (meshArray type)
"""
function matrixmean(w,froot,years,γ)
    nyr = length(years)

    fname = froot*string(years[1])
    println("initialize θout")
    field = read_bin(fname,Float32,γ)

    θout = MeshArray(γ,Float32) # some nans here
    fill!(θout,0.0)
    println("initialization finished")

    istart = 1
    for tt in 1:nyr
        fname = froot*string(years[tt])
        println(fname)
        field = read_bin(fname,Float32,γ)
        inx, int = size(field)
        for i in 1:int
            θout.f[:] .+= field.f[:, i] .* w
        end
    end

    return θout
end
