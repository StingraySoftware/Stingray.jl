using Stingray
using Documenter

DocMeta.setdocmeta!(Stingray, :DocTestSetup, :(using Stingray); recursive=true)

makedocs(;
    modules=[Stingray],
    authors="Aman Pandey",
    repo=Remotes.GitHub("StingraySoftware", "Stingray.jl"),
    sitename="Stingray.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://github.com/StingraySoftware/Stingray.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="https://github.com/StingraySoftware/Stingray.jl",
    devbranch="main",
    target="build",
    push_preview=true
)
