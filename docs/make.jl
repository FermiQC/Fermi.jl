using Documenter, Fermi

makedocs(sitename="Fermi Documentation",format = Documenter.HTML(prettyurls=false))
deploydocs(
           repo = "github.com/FermiQC/Fermi.jl.git",
          )
