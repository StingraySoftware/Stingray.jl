using Stingray
using Documenter

DocMeta.setdocmeta!(Stingray, :DocTestSetup, :(using Stingray); recursive=true)

makedocs(;
    modules=[Stingray],
    authors="Aman Pandey",
    repo="https://github.com/matteobachetti/Stingray.jl/blob/{commit}{path}#{line}",
    sitename="Stingray.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://github.com/matteobachetti/Stingray.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="https://github.com/matteobachetti/Stingray.jl",
    devbranch="main",
)
