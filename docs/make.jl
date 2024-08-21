using Documenter
# Using a placeholder module 

# a dummy or minimal package to test documentation
module MyPDELib
    export hello
    
    function hello()
        println("Hello from MyPDELib!")
    end
end

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

# Generate the documentation files
Documenter.write_directory()
