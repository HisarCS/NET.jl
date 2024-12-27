using NET
using Documenter

DocMeta.setdocmeta!(NET, :DocTestSetup, :(using NET); recursive=true)

makedocs(;
    modules=[NET],
    authors="Emre Dayangaç",
    sitename="NET.jl",
    format=Documenter.HTML(;
        canonical="https://HisarCS.github.io/NET.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/EmreDay1/NET.jl",
    devbranch="main",
)
