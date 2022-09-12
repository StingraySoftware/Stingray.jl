using Stingray
using Documenter

DocMeta.setdocmeta!(Stingray, :DocTestSetup, :(using Stingray); recursive=true)

makedocs(;
    modules=[Stingray],
    sitename="Stingray.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://github.com/matteobachetti/Stingray.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "API Reference" => "api.md",
    ],
    strict=true,
    checkdocs=:export
)

deploydocs(;
    repo="https://github.com/StingraySoftware/Stingray.jl",
)
