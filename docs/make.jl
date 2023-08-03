using ECCOonPoseidon
using Documenter

DocMeta.setdocmeta!(ECCOonPoseidon, :DocTestSetup, :(using ECCOonPoseidon); recursive=true)

makedocs(;
    modules=[ECCOonPoseidon],
    authors="G Jake Gebbie",
    repo="https://github.com/ggebbie/ECCOonPoseidon.jl/blob/{commit}{path}#{line}",
    sitename="ECCOonPoseidon.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://ggebbie.github.io/ECCOonPoseidon.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ggebbie/ECCOonPoseidon.jl",
    devbranch="main",
)
