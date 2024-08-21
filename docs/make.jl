using Documenter
# Using a placeholder module 

# a dummy or minimal package to test documentation


# Include the documentation files
makedocs(
    sitename = "PDE Library Documentation",
    format = Documenter.HTML(),  
    pages = [
        "Home" => "index.md",
        "Manual" => "manual.md",
        "Reference" => "reference.md"
    ]
)
