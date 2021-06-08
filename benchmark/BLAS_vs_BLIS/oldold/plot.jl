function get_tblis_data(N, C)

    open(joinpath(@__DIR__, "N$N", "tblis_c$(C)_t$(N).out"), "r") do lines

        for line = eachline(lines)
            m = match(r"Average time per iteration\s+?(\d+?\.\d+)", line)
            if m !== nothing
                return parse(Float64, m.captures[1])
            end
        end
    end
end

