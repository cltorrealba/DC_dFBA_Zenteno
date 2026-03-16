using Libdl
bin = raw"C:\Users\ctorrealba\CoinHSL.v2023\bin"
libs = [joinpath(bin, "libhsl.dll"), joinpath(bin, "libcoinhsl.dll")]
println("bin exists: ", isdir(bin))
for lib in libs
    println("---")
    println("lib=", lib)
    if !isfile(lib)
        println("exists=false")
        continue
    end
    try
        h = Libdl.dlopen(lib)
        sym = Libdl.dlsym(h, :LIBHSL_isfunctional)
        ok = ccall(sym, Cint, ()) == 1
        println("functional=", ok)
        Libdl.dlclose(h)
    catch err
        println("functional=false (", err, ")")
    end
end
