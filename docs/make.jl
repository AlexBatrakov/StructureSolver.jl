using Documenter
using StructureSolver

DocMeta.setdocmeta!(StructureSolver, :DocTestSetup, :(using StructureSolver); recursive=true)

makedocs(
    sitename = "StructureSolver.jl",
    modules = [StructureSolver],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://AlexBatrakov.github.io/StructureSolver.jl",
        edit_link = "main",
    ),
    pages = [
        "Home" => "index.md",
        "Quickstart" => "quickstart.md",
        "Tutorial: M–R curve" => "tutorial_mr.md",
        "EoS" => "eos.md",
        "Simulations" => "simulations.md",
        "Examples" => "examples.md",
        "Troubleshooting" => "troubleshooting.md",
        "API" => "api.md",
    ],
)

if get(ENV, "GITHUB_ACTIONS", "false") == "true"
    deploydocs(
        repo = "github.com/AlexBatrakov/StructureSolver.jl",
        devbranch = "main",
        push_preview = true,
    )
end
