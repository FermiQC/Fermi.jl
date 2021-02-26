module Output
import Fermi
using Formatting
#using Fermi.Options

#printstyle = Fermi.Options.printstyle
if !isdefined(Fermi.Output,:printstyle)
    if isinteractive()
        printstyle = ["stdout"]
    else
        printstyle = ["file"]
    end
end

export output
export @output

"""
    set_print(pstyle)

Used to set the print style for @output and output calls. 
Options are "none","file","stdout","both"
"""
function set_print(pstyle)
    if pstyle in ["none","file","stdout","both"]
        Fermi.Output.printstyle[1] = pstyle
        #revise(Fermi.Output)
    else
        return false
    end
end

function output(str,x...)
    if Fermi.Output.printstyle[1] == "stdout"
        f = format(str,x...)
        print(f)
    elseif Fermi.Output.printstyle[1] == "file"
        f = format(str,x...)
        open("output.dat","a") do file
            write(file,f)
            flush(file)
        end
    elseif Fermi.Output.printstyle[1] == "both"
        f = format(str,x...)
        print(f)
        open("output.dat","a") do file
            write(file,f)
            flush(file)
        end
    elseif Fermi.Output.printstyle[1] == "none"
    end
end

macro output(str,x...)
    return quote
        if Fermi.Output.printstyle[1] == "stdout"
                local f = format($str, $([esc(i) for i in x]...))
                print(f)
        elseif Fermi.Output.printstyle[1] == "file"
                local f = format($str, $([esc(i) for i in x]...))
                open("output.dat","a") do file
                    write(file,f)
                    flush(file)
                end
        elseif Fermi.Output.printstyle[1] == "both"
            if length(x) >= 1
                local f = format($str, $([esc(i) for i in x]...))
                open("output.dat","a") do file
                    write(file,f)
                    flush(file)
                end
                print(f)
            else
                local f = format($str)
                open("output.dat","a") do file
                    write(file,f)
                    flush(file)
                end
            end
        elseif Fermi.Output.printstyle[1] == "none"
        end
    end
end
end #module
