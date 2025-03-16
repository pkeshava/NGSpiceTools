push!(LOAD_PATH, "../src/")

using Documenter

using NGSpiceTools  # Now bring it into scope


makedocs(
    sitename = "NGSpiceTools.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
    ),
    modules = [NGSpiceTools],  # Make sure this is included
    authors = "Pouyan Keshavarzian",
    pages = [
        "Home" => "index.md",
        "User Guide" => "guide.md",
        "API Reference" => "reference.md"
    ],
    repo = "",
    remotes = nothing  # Disable remote checks if not using Git
)

# Documenter can also automatically deploy docs to GitHub pages.
# See the Documenter.jl README for more information.
# deploydocs(
#     repo = "github.com/pkeshava/NGSpiceTools.jl.git",
#     devbranch = "main"
# )