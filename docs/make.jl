using MyFirstPackage
using Documenter

DocMeta.setdocmeta!(MyFirstPackage, :DocTestSetup, :(using MyFirstPackage); recursive=true)

makedocs(;
    modules=[MyFirstPackage],
    authors="JChen901",
    sitename="MyFirstPackage.jl",
    format=Documenter.HTML(;
        canonical="https://JChen901.github.io/MyFirstPackage.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JChen901/MyFirstPackage.jl",
    devbranch="main",
)
