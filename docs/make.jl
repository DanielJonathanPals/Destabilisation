using Destabilisation
using Documenter

DocMeta.setdocmeta!(Destabilisation, :DocTestSetup, :(using Destabilisation); recursive=true)

makedocs(;
    modules=[Destabilisation],
    authors="Daniel Pals <Daniel.Pals@tum.de>",
    repo="https://github.com/DanielJonathanPals/Destabilisation.jl/blob/{commit}{path}#{line}",
    sitename="Destabilisation.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "home.md",
        "Index" => "index.md",
        "Format Test" => "formatTests.md",
        "Matrices" => "matrices.md",
        "Coupling to dynamical system" => "coupling.md",
        "Integrate trajectories" => "integrate.md",
        "Fit VAR model" => "fitVARmodel.md",
    ],
)

deploydocs(
    repo = "github.com/DanielJonathanPals/Destabilisation",
)
