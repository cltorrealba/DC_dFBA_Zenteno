#!/usr/bin/env julia

using Pkg

const HSL_IMPORT_ERROR = Ref{Union{Nothing,String}}(nothing)
const HAS_HSL_JLL = let
    try
        @eval import HSL_jll
        true
    catch err
        HSL_IMPORT_ERROR[] = sprint(showerror, err)
        false
    end
end

function banner(msg)
    println("\n", "="^80)
    println(msg)
    println("="^80)
end

function list_existing(paths::Vector{String})
    return [p for p in paths if ispath(p)]
end

function find_hsl_package_roots()
    roots = String[]

    for depot in DEPOT_PATH
        pkg_dir = joinpath(depot, "packages", "HSL_jll")
        if isdir(pkg_dir)
            for entry in readdir(pkg_dir; join=true)
                if isdir(entry)
                    push!(roots, entry)
                end
            end
        end
    end

    unique!(roots)
    return sort(roots)
end

function find_hsl_artifacts()
    dlls = String[]
    patterns = ("libhsl.dll", "libhsl_subset.dll", "libhsl_subset_64.dll", "coinhsl.dll")

    for depot in DEPOT_PATH
        artifacts_dir = joinpath(depot, "artifacts")
        if !isdir(artifacts_dir)
            continue
        end
        for (root, _, files) in walkdir(artifacts_dir)
            for f in files
                if f in patterns
                    push!(dlls, joinpath(root, f))
                end
            end
        end
    end

    unique!(dlls)
    return sort(dlls)
end

function package_root_from_pathof(pathof_hsl::AbstractString)
    # pathof(HSL_jll) -> .../HSL_jll/<hash>/src/HSL_jll.jl
    return normpath(joinpath(dirname(pathof_hsl), ".."))
end

function try_load_hsl_jll()
    try
        if !HAS_HSL_JLL
            return (
                loaded = false,
                pathof = nothing,
                package_root = nothing,
                libhsl_path = nothing,
                functional = false,
                error = HSL_IMPORT_ERROR[],
            )
        end
        functional = ccall((:LIBHSL_isfunctional, HSL_jll.libhsl), Cint, ()) != 0
        return (
            loaded = true,
            pathof = pathof(HSL_jll),
            package_root = package_root_from_pathof(pathof(HSL_jll)),
            libhsl_path = String(HSL_jll.libhsl_path),
            functional = functional,
            error = nothing,
        )
    catch err
        return (
            loaded = false,
            pathof = nothing,
            package_root = nothing,
            libhsl_path = nothing,
            functional = false,
            error = sprint(showerror, err),
        )
    end
end

function main()
    banner("Julia / Environment")
    println("VERSION              = ", VERSION)
    println("ACTIVE_PROJECT       = ", Base.active_project())
    println("DEPOT_PATH (existing)")
    for d in list_existing(String.(DEPOT_PATH))
        println("  - ", d)
    end

    banner("Search: HSL_jll package cache")
    pkg_roots = find_hsl_package_roots()
    if isempty(pkg_roots)
        println("No package roots found under DEPOT_PATH")
    else
        for p in pkg_roots
            println("  - ", p)
        end
    end

    banner("Search: HSL DLL artifacts")
    dlls = find_hsl_artifacts()
    if isempty(dlls)
        println("No libhsl*.dll / coinhsl.dll artifacts found under DEPOT_PATH")
    else
        for d in dlls
            println("  - ", d)
        end
    end

    banner("Load test: import HSL_jll")
    res = try_load_hsl_jll()
    if res.loaded
        println("loaded               = true")
        println("pathof(HSL_jll)      = ", res.pathof)
        println("package_root         = ", res.package_root)
        println("libhsl_path          = ", res.libhsl_path)
        println("LIBHSL_isfunctional  = ", res.functional)
    else
        println("loaded               = false")
        println("error                = ", res.error)
    end

    banner("Suggested notebook bootstrap lines")
    println("import Pkg")
    println("using LinearAlgebra")
    println("Pkg.activate(pwd())")
    if res.package_root !== nothing
        println("ENV[\"HSL_JLL_PATH\"] = raw\"", res.package_root, "\"")
        println("ENV[\"IPOPT_LINEAR_SOLVER\"] = \"ma86\"   # options: ma57 | ma77 | ma86")
        println("Pkg.develop(path=ENV[\"HSL_JLL_PATH\"])")
        println("Pkg.resolve()")
        println("LinearAlgebra.BLAS.set_num_threads(1)")
    elseif !isempty(pkg_roots)
        println("ENV[\"HSL_JLL_PATH\"] = raw\"", first(pkg_roots), "\"")
        println("ENV[\"IPOPT_LINEAR_SOLVER\"] = \"ma86\"   # options: ma57 | ma77 | ma86")
        println("Pkg.develop(path=ENV[\"HSL_JLL_PATH\"])")
        println("Pkg.resolve()")
        println("LinearAlgebra.BLAS.set_num_threads(1)")
    else
        println("No HSL_jll root found automatically; set HSL_JLL_PATH manually.")
        println("ENV[\"IPOPT_LINEAR_SOLVER\"] = \"ma86\"   # options: ma57 | ma77 | ma86")
        println("LinearAlgebra.BLAS.set_num_threads(1)")
    end

    banner("Expected runtime trace in solver")
    println("[ipopt] linear_solver requested=ma57 active=ma57 hsl_ready=true")
    println("[ipopt] linear_solver requested=ma77 active=ma77 hsl_ready=true")
    println("[ipopt] linear_solver requested=ma86 active=ma86 hsl_ready=true")
    println("If hsl_ready=false, fallback to mumps is expected.")

    println("\nDone.")
end

main()
