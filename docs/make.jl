using Documenter, Fermi

Basics = "Basics" => "index.md"
Methods = "Methods" => ["hartreefock.md", "mp.md", "cc.md"]
Index = "Index" => "indice.md"

PAGES = [
    Basics,
    Methods,
    Index
]

makedocs(
    sitename="Fermi.jl",
    authors = "Gustavo Aroeira",
    format = Documenter.HTML(
        sidebar_sitename = false
    ),
    pages = PAGES
)

deploydocs(
    repo = "github.com/FermiQC/Fermi.jl.git",
)
