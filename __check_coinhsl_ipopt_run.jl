using JuMP, Ipopt
import MathOptInterface as MOI

ENV["IPOPT_HSLLIB"] = raw"C:\Users\ctorrealba\CoinHSL.v2023\bin\libhsl.dll"
ENV["PATH"] = dirname(ENV["IPOPT_HSLLIB"]) * ';' * get(ENV, "PATH", "")

include("pFBA_KKT_flux_Zenteno_vargam_simultaneous_sparsepatch.jl")

for s in ("ma57", "ma77", "ma86")
    println("\n=== IPOPT run with ", s, " ===")
    active, hsllib, _ = _configure_ipopt_linear_solver(s)
    m = Model(Ipopt.Optimizer)
    set_attribute(m, "print_level", 0)
    set_attribute(m, "linear_solver", active)
    hsllib !== nothing && set_attribute(m, "hsllib", hsllib)
    @variable(m, x >= 1.0)
    @objective(m, Min, (x - 2.0)^2)
    optimize!(m)
    term = termination_status(m)
    xv = value(x)
    println("termination=", term, " x=", xv, " active=", active)
end
