using LinearAlgebra
const IPOPT_PRINT_LEVEL = 0
const IPOPT_TOL = 1e-8
const IPOPT_ACCEPTABLE_TOL = 1e-6
const IPOPT_ACCEPTABLE_ITER = 5
const IPOPT_MAX_ITER = 50
const IPOPT_CONSTR_VIOL_TOL = 1e-6
const IPOPT_COMPL_INF_TOL = 1e-6
const IPOPT_LINEAR_SOLVER = "ma86"

ENV["IPOPT_HSLLIB"] = raw"C:\Users\ctorrealba\CoinHSL.v2023\bin\libhsl.dll"
ENV["PATH"] = dirname(ENV["IPOPT_HSLLIB"]) * ';' * get(ENV, "PATH", "")

include("pFBA_KKT_flux_Zenteno_vargam_simultaneous_sparsepatch.jl")

for s in ("ma57", "ma77", "ma86")
    println("\n=== request ", s, " ===")
    active, hsllib, ready = _configure_ipopt_linear_solver(s)
    println("active=", active)
    println("ready=", ready)
    println("hsllib=", hsllib)
end
