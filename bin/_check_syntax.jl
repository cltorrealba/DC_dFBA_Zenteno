for f in ["main.jl", "pFBA_KKT_flux_Zenteno_vargam_simultaneous_v2.jl"]
    fp = joinpath(@__DIR__, f)
    try
        code = read(fp, String)
        ex = Meta.parseall(code)
        println("OK  $f  ($(length(code)) chars)")
    catch e
        println("ERR $f")
        showerror(stdout, e)
        println()
    end
end
