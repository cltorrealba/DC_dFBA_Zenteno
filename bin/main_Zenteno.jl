#=
    main_Zenteno.jl
    ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬ÃƒÂ¢Ã¢â‚¬ÂÃ¢â€šÂ¬
    SimulaciÃƒÆ’Ã‚Â³n de dFBA por colocaciÃƒÆ’Ã‚Â³n ortogonal + KKT embebido (MPCC/pFBA)
    Adaptado al modelo cinÃƒÆ’Ã‚Â©tico de Zenteno con yeast-GEM v9.0.2

    AnÃƒÆ’Ã‚Â¡logo a main.jl (Scerevisiae_iND750) pero con:
      - 6 estados: X, N, G, F, E, O2
      - Medio mÃƒÆ’Ã‚Â­nimo vÃƒÆ’Ã‚Â­nico (23 aminoÃƒÆ’Ã‚Â¡cidos, vitaminas, iones)
      - CinÃƒÆ’Ã‚Â©tica de Zenteno (temperatura, inhibiciÃƒÆ’Ã‚Â³n, muerte celular)
      - Bypass anaerÃƒÆ’Ã‚Â³bico para cofactores (Heme/CoQ)

    Uso:
      cd Yeast_83/estima
      julia main_Zenteno.jl
      # O con GEM v8:
      GEM_VERSION=v8 julia main_Zenteno.jl
=#

#ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢#
# Paquetes
using JuMP
using Ipopt
using LinearAlgebra
using DelimitedFiles
using Printf

#ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢ÃƒÂ¢Ã‹â€ Ã¢â‚¬â„¢#
# Incluir la funciÃƒÆ’Ã‚Â³n MPCC
include("pFBA_KKT_flux_Zenteno.jl")

#=============================================================================
  1. CARGAR GEM
=============================================================================#
const GEM_VERSION = Symbol(get(ENV, "GEM_VERSION", "v9"))
const ESTIMA_DIR  = @__DIR__

println(">>> Cargando GEM ($(GEM_VERSION))...")
if GEM_VERSION == :v9
    S     = readdlm(joinpath(ESTIMA_DIR, "S_v9.csv"), ',')
    lbraw = readdlm(joinpath(ESTIMA_DIR, "lb_v9.csv"), ',')
    ubraw = readdlm(joinpath(ESTIMA_DIR, "ub_v9.csv"), ',')
    const RXN_IDS = readlines(joinpath(ESTIMA_DIR, "rxn_ids_v9.txt"))
    const MET_IDS = readlines(joinpath(ESTIMA_DIR, "met_ids_v9.txt"))
else
    S     = readdlm(joinpath(ESTIMA_DIR, "S_hen.csv"), ',')
    lbraw = readdlm(joinpath(ESTIMA_DIR, "lb_hen.csv"), ',')
    ubraw = readdlm(joinpath(ESTIMA_DIR, "ub_hen.csv"), ',')
    const RXN_IDS = readlines(joinpath(ESTIMA_DIR, "rxn_ids.txt"))
    const MET_IDS = readlines(joinpath(ESTIMA_DIR, "met_ids.txt"))
end

vlb = lbraw isa AbstractVector ? Float64.(lbraw) : Float64.(lbraw[:,1])
vub = ubraw isa AbstractVector ? Float64.(ubraw) : Float64.(ubraw[:,1])
nm  = size(S, 1)
nv  = size(S, 2)

println("   GEM cargado: $(nm) metabolitos, $(nv) reacciones")

#=============================================================================
  2. DICCIONARIOS E ÃƒÆ’Ã‚ÂNDICES DE REACCIONES CLAVE
=============================================================================#
const RXN_INDEX = Dict{String,Int}(rxn => i for (i, rxn) in pairs(RXN_IDS))
const MET_INDEX = Dict{String,Int}(met => i for (i, met) in pairs(MET_IDS))

get_rxn(name::AbstractString) = get(RXN_INDEX, String(name)) do
    error("ReacciÃƒÆ’Ã‚Â³n '$(String(name))' no encontrada en RXN_INDEX")
end
get_met(name::AbstractString) = get(MET_INDEX, String(name)) do
    error("Metabolito '$(String(name))' no encontrado en MET_INDEX")
end

# --- Reacciones clave ---
const eth = get_rxn("r_1761")   # intercambio etanol
const obj = get_rxn("r_4041")   # biomasa
const glu = get_rxn("r_1714")   # uptake glucosa
const fru = get_rxn("r_1709")   # uptake fructosa
const o2  = get_rxn("r_1992")   # intercambio O2
const IDX_ATPM = get_rxn("r_4046")  # mantenimiento ATP

println("   ÃƒÆ’Ã‚Ândices: obj=$(obj) glu=$(glu) fru=$(fru) eth=$(eth) o2=$(o2)")

#=============================================================================
  3. FUENTES DE NITRÃƒÆ’Ã¢â‚¬Å“GENO
=============================================================================#
# 8 fuentes cinÃƒÆ’Ã‚Â©ticamente acopladas (modelo de Zenteno)
const KINETIC_N_SOURCE_IDS = [
    "r_1654",  # NH4+
    "r_1879",  # Arg
    "r_1891",  # Gln
    "r_1889",  # Glu
    "r_1906",  # Ser
    "r_1911",  # Thr
    "r_1873",  # Ala
    "r_1912",  # Trp
]

# 23 fuentes del medio vÃƒÆ’Ã‚Â­nico (para bounds del GEM v9)
const ALL_N_SOURCE_IDS = [
    "r_1654", "r_1904", "r_1891", "r_1889", "r_1880", "r_1881",
    "r_1873", "r_1879", "r_1810", "r_1906", "r_1911", "r_1914",
    "r_1899", "r_1897", "r_1903", "r_1913", "r_1912", "r_1893",
    "r_1900", "r_1902", "r_1883", "r_1987", "r_1800",
]

# ÃƒÆ’Ã‚Ândices de N cinÃƒÆ’Ã‚Â©ticamente acoplados
const NITROGEN_SOURCES = [get_rxn(rid) for rid in KINETIC_N_SOURCE_IDS]
const UPTAKE_IDXS = vcat([glu, fru], NITROGEN_SOURCES)
const n_up = length(UPTAKE_IDXS)

# Selectores que NO dependen de nv (se pueden crear ahora)
const IS_GLU = [x == glu ? 1.0 : 0.0 for x in UPTAKE_IDXS]
const IS_FRU = [x == fru ? 1.0 : 0.0 for x in UPTAKE_IDXS]
const IS_NIT = [1.0 - IS_GLU[i] - IS_FRU[i] for i in 1:n_up]

# ÃƒÆ’Ã‚Âtomos de N por molÃƒÆ’Ã‚Â©cula y perfil de consumo (por ÃƒÆ’Ã‚Â­ndice de reacciÃƒÆ’Ã‚Â³n)
const idx_NH4 = get_rxn("r_1654")
const idx_Arg = get_rxn("r_1879")
const idx_Gln = get_rxn("r_1891")
const idx_Glu_met = get_rxn("r_1889")
const idx_Ser = get_rxn("r_1906")
const idx_Thr = get_rxn("r_1911")
const idx_Ala = get_rxn("r_1873")
const idx_Trp = get_rxn("r_1912")

const N_atoms_dict = Dict(idx_NH4 => 1.0, idx_Arg => 4.0, idx_Gln => 2.0, idx_Glu_met => 1.0,
                          idx_Ser => 1.0, idx_Thr => 1.0, idx_Ala => 1.0, idx_Trp => 2.0)
const N_profile_dict = Dict(idx_NH4 => 0.40, idx_Arg => 0.20, idx_Gln => 0.10, idx_Glu_met => 0.05,
                            idx_Ser => 0.05, idx_Thr => 0.05, idx_Ala => 0.05, idx_Trp => 0.10)

# Constantes de masa molar (independientes de nv)
const MW_N   = 0.014007  # g/mmol
const MW_GLU = 0.180156  # g/mmol
const MW_FRU = 0.180156  # g/mmol
const MW_ETH = 0.046070  # g/mmol
const MW_O2  = 0.031998  # g/mmol

# NOTA: SELECT_UPTAKE, N_atoms_vec, N_profile_vec y N_frac se crean
# DESPUÃƒÆ’Ã¢â‚¬Â°S de inyectar bypass (ver secciÃƒÆ’Ã‚Â³n 4b) para usar el nv final.

#=============================================================================
  4. INYECTAR BYPASS ANAERÃƒÆ’Ã¢â‚¬Å“BICO (Heme/CoQ)
=============================================================================#
println(">>> Inyectando reacciones bypass (Heme/CoQ)...")
const METS_ANAEROBIC_BYPASS = [
    "s_3714[c]", "s_1198[c]", "s_1203[c]",
    "s_1207[c]", "s_1212[c]", "s_0529[c]"
]
n_bypass = length(METS_ANAEROBIC_BYPASS)
S_bypass = zeros(nm, n_bypass)
for (k, mid) in enumerate(METS_ANAEROBIC_BYPASS)
    if haskey(MET_INDEX, mid)
        S_bypass[get_met(mid), k] = 1.0
    else
        clean = replace(mid, "[c]" => "")
        haskey(MET_INDEX, clean) && (S_bypass[get_met(clean), k] = 1.0)
    end
end

nv_old = nv
S  = hcat(S, S_bypass)
nv = size(S, 2)
nm = size(S, 1)
vlb = vcat(vlb, zeros(n_bypass))
vub = vcat(vub, zeros(n_bypass))

const IDX_BYPASS_START = nv_old + 1
const IDX_BYPASS_END   = nv
println("   Bypass: reacciones $(IDX_BYPASS_START):$(IDX_BYPASS_END)")

# Eliminar cofactores de reacciÃƒÆ’Ã‚Â³n r_4598 (biomasa requiere estos en aerobiosis)
rxn_cofactor = get_rxn("r_4598")
for mid in METS_ANAEROBIC_BYPASS
    target = haskey(MET_INDEX, mid) ? mid : replace(mid, "[c]" => "")
    haskey(MET_INDEX, target) && (S[get_met(target), rxn_cofactor] = 0.0)
end

#--- 4b. Vectores que dependen de nv FINAL (post-bypass) ---
# Selector de uptake: SELECT_UPTAKE[mc,k] = 1 si reacciÃƒÆ’Ã‚Â³n mc es la fuente k
SELECT_UPTAKE = [Float64(mc == UPTAKE_IDXS[k]) for mc in 1:nv, k in 1:n_up]

# Vectores de N indexados por reacciÃƒÆ’Ã‚Â³n (nv elementos)
N_atoms_vec = ones(nv)
N_profile_vec = zeros(nv)
for k in keys(N_atoms_dict);   1 <= k <= nv && (N_atoms_vec[k] = N_atoms_dict[k]);   end
for k in keys(N_profile_dict); 1 <= k <= nv && (N_profile_vec[k] = N_profile_dict[k]); end

# FracciÃƒÆ’Ã‚Â³n de N pre-calculada: convierte gN/gDW/h ÃƒÂ¢Ã¢â‚¬Â Ã¢â‚¬â„¢ mmol_metabolito/gDW/h
N_frac = zeros(n_up)
for k in 1:n_up
    mc = UPTAKE_IDXS[k]
    if IS_NIT[k] > 0.5
        N_frac[k] = N_profile_vec[mc] / (N_atoms_vec[mc] * MW_N)
    end
end

#--- 4c. Productos con cotas dinÃƒÆ’Ã‚Â¡micas (etanol + biomasa) ---
const PRODUCT_IDXS = [eth, obj]
const n_prod = length(PRODUCT_IDXS)
SELECT_PRODUCT = [Float64(mc == PRODUCT_IDXS[k]) for mc in 1:nv, k in 1:n_prod]
const IS_ETH_prod = [1.0, 0.0]  # k=1 ÃƒÂ¢Ã¢â‚¬Â Ã¢â‚¬â„¢ etanol
const IS_OBJ_prod = [0.0, 1.0]  # k=2 ÃƒÂ¢Ã¢â‚¬Â Ã¢â‚¬â„¢ biomasa
println("   Productos dinÃƒÆ’Ã‚Â¡micos: etanol (idx $(eth)), biomasa (idx $(obj))")

#=============================================================================
  5. CONFIGURAR MODELO: BOUNDS Y MEDIO VÃƒÆ’Ã‚ÂNICO
=============================================================================#
println(">>> Configurando medio vÃƒÆ’Ã‚Â­nico y bounds...")

# --- O2: anaerobiosis estricta ---
# Sin intercambio de O2: v_o2 = 0 en toda la simulacion
vlb[o2] = 0.0
vub[o2] = 0.0
@assert vlb[o2] == 0.0 && vub[o2] == 0.0 "Modo anaerobico requiere v_o2 = 0"

# --- Suplementos anaerÃƒÆ’Ã‚Â³bicos abiertos ---
for rid in ["r_1757", "r_1915", "r_1994", "r_2106", "r_2134", "r_2137", "r_2189"]
    haskey(RXN_INDEX, rid) && (vlb[get_rxn(rid)] = -1000.0)
end

# --- Shuttles cerrados (fuerza metabolismo fermentativo) ---
haskey(RXN_INDEX, "r_0713") && (vlb[get_rxn("r_0713")] = 0.0)
haskey(RXN_INDEX, "r_0714") && (vlb[get_rxn("r_0714")] = 0.0)
haskey(RXN_INDEX, "r_0487") && (vub[get_rxn("r_0487")] = 0.0)

# --- Fuentes de carbono ---
vlb[glu] = -1000.0; vub[glu] = 0.0
vlb[fru] = -1000.0; vub[fru] = 0.0

# --- TODAS las fuentes de nitrÃƒÆ’Ã‚Â³geno (23 AA del medio vÃƒÆ’Ã‚Â­nico) ---
for rid in ALL_N_SOURCE_IDS
    haskey(RXN_INDEX, rid) || continue
    idx = get_rxn(rid)
    vlb[idx] = -1000.0; vub[idx] = 0.0
end

# --- Vitaminas ---
const VITAMIN_IDS = ["r_2067", "r_1671", "r_1967", "r_2028", "r_1548", "r_2038", "r_1947"]
for rid in VITAMIN_IDS
    haskey(RXN_INDEX, rid) || continue
    idx = get_rxn(rid)
    vlb[idx] = -1000.0; vub[idx] = 0.0
end

# --- Iones esenciales ---
const ION_IDS = ["r_2005", "r_2060", "r_2020", "r_1861", "r_4596", "r_4597",
                 "r_4594", "r_4600", "r_4595", "r_4593"]
for rid in ION_IDS
    haskey(RXN_INDEX, rid) || continue
    idx = get_rxn(rid)
    vlb[idx] = -1000.0; vub[idx] = 0.0
end

# --- Productos abiertos ---
const PRODUCT_IDS = ["r_1761", "r_1808", "r_1634", "r_2056", "r_1549", "r_1546",
                     "r_1552", "r_1765", "r_1867", "r_1866", "r_1862", "r_1865"]
for rid in PRODUCT_IDS
    haskey(RXN_INDEX, rid) || continue
    idx = get_rxn(rid)
    vlb[idx] = 0.0; vub[idx] = 1000.0
end

# --- Agua y H+ libres ---
haskey(RXN_INDEX, "r_2100") && (vlb[get_rxn("r_2100")] = -1000.0; vub[get_rxn("r_2100")] = 1000.0)
haskey(RXN_INDEX, "r_1832") && (vlb[get_rxn("r_1832")] = -1000.0; vub[get_rxn("r_1832")] = 1000.0)

# --- ATP mantenimiento flexible ---
if length(vlb) >= IDX_ATPM
    vlb[IDX_ATPM] = 0.7
    vub[IDX_ATPM] = 0.7
    println("   ATPM (idx $(IDX_ATPM)): lb=$(vlb[IDX_ATPM]), ub=$(vub[IDX_ATPM])")
end

#=============================================================================
  6. PARÃƒÆ’Ã‚ÂMETROS CINÃƒÆ’Ã¢â‚¬Â°TICOS DE ZENTENO
=============================================================================#
const ZENTENO_PARAM_SET = lowercase(get(ENV, "ZENTENO_PARAM_SET", "paper2010"))
if ZENTENO_PARAM_SET == "paper2010"
    # Valores originales de Zenteno et al. (2010), Tabla 2.
    # Unidades en este modelo: g/L, g/g y 1/h.
    const MU0_nom    = 0.18      # l0 (1/h)
    const YXN_nom    = 19.69     # YX/N (kgBio/kgN)
    const YXG_nom    = 1.60      # YX/G (kgBio/kgG)
    const YXF_nom    = 1.60      # YX/F (kgBio/kgF)
    const YEG_nom    = 0.49      # YE/G (kgE/kgG)
    const YEF_nom    = 0.49      # YE/F (kgE/kgF)
    const Kn0_nom    = 0.01      # KN0 (kgN/m3)
    const Kg0_nom    = 7.5       # KG0 (kgG/m3)
    const Kf0_nom    = 7.5       # KF0 (kgF/m3)
    const Kig0_nom   = 55.0      # KIG0 = S/4 con S_ref=220 g/L
    const Kie0_nom   = 40.0      # KIE0 (kgE/m3)
    const Kd0_nom    = 0.00044   # kD0 (1/h)
    const betaG0_nom = 0.225     # bG0 (kgE/kgBio/h)
    const betaF0_nom = 0.225     # bF0 (kgE/kgBio/h)
elseif ZENTENO_PARAM_SET == "legacy_calibrated"
    # Set previo calibrado localmente
    const MU0_nom    = 0.141665
    const YXN_nom    = 9.80576
    const YXG_nom    = 0.394345
    const YXF_nom    = 0.18622
    const YEG_nom    = 0.14133
    const YEF_nom    = 0.96932
    const Kn0_nom    = 0.226882
    const Kg0_nom    = 3.1514
    const Kf0_nom    = 2.97625
    const Kig0_nom   = 29.5276
    const Kie0_nom   = 2.99809
    const Kd0_nom    = 0.0000311736
    const betaG0_nom = 1.41182
    const betaF0_nom = 8.49482
else
    error("ZENTENO_PARAM_SET invalido: $(ZENTENO_PARAM_SET). Usa 'paper2010' o 'legacy_calibrated'.")
end
println("   Set de parametros Zenteno: $(ZENTENO_PARAM_SET)")

# ParÃƒÆ’Ã‚Â¡metros usados en la simulaciÃƒÆ’Ã‚Â³n (fijos, sin estimaciÃƒÆ’Ã‚Â³n)
const MU0     = MU0_nom
const YEG     = YEG_nom
const YEF     = YEF_nom
const YXN     = YXN_nom
const MRATE_0 = 0.01     # Tasa de mantenimiento base (g sustrato/gDW/h)

#=============================================================================
  7. FUNCIONES AUXILIARES (registradas en JuMP)
=============================================================================#
const R_GAS = 8.314  # Compartido con pFBA como 'R'
const R = R_GAS
const EPS = 1e-9

# --- Temperatura dinÃƒÆ’Ã‚Â¡mica ---
const T_BASE  = try parse(Float64, get(ENV, "T_CONST", "293.15")) catch; 293.15 end
const T_STEPS = [36.0, 96.0]
const T_DELTAS = [5.0, 3.0]
const T_STEEP = 0.5

function dynamic_temperature(t)
    val = T_BASE
    for idx in eachindex(T_STEPS)
        sigmoid = 1.0 / (1.0 + exp(-T_STEEP * (t - T_STEPS[idx])))
        val += T_DELTAS[idx] * sigmoid
    end
    return val
end

# --- Tasa de muerte celular ---
function death_rate_T(E, T_val)
    Td = -0.0001 * E^3 + 0.0049 * E^2 - 0.1279 * E + 315.89
    s = 0.5 * (1.0 + tanh(0.5 * (T_val - Td)))
    base = Kd0_nom * exp(0.0415 * E + (130000.0 * (T_val - 305.65)) / (305.65 * R * T_val))
    return base * s
end

# --- InyecciÃƒÆ’Ã‚Â³n de nutrientes (pulso gaussiano) ---
const SQRT_2PI = sqrt(2.0 * pi)
function smooth_injection(t, t_shot, dose, width=1.0)
    abs(t - t_shot) > 5 * width && return 0.0
    return (dose / (width * SQRT_2PI)) * exp(-0.5 * ((t - t_shot) / width)^2)
end

# ConfiguraciÃƒÆ’Ã‚Â³n de pulsos (neutros por defecto)
const T_INJ_1 = 0.0;  const DOSE_1 = 0.0;  const WIDTH_1 = 5.0
const T_INJ_2 = 0.0;  const DOSE_2 = 0.0;  const WIDTH_2 = 5.0

#=============================================================================
  8. CONSTANTES CINÃƒÆ’Ã¢â‚¬Â°TICAS Y FÃƒÆ’Ã‚ÂSICAS
=============================================================================#

# Parametros numericos del MPCC
const O2_SAT = 0.26  # mmol/L
const EPS_FLUX   = try parse(Float64, get(ENV, "EPS_FLUX", "1e-6")) catch; 0.0 end
const O2_ZERO_TOL = try parse(Float64, get(ENV, "O2_ZERO_TOL", "1e-7")) catch; 1e-7 end

const LINEAR_SOLVER = lowercase(get(ENV, "IPOPT_LINEAR_SOLVER", "mumps"))

#=============================================================================
  9. CONFIGURACIÃƒÆ’Ã¢â‚¬Å“N pFBA / MPCC
=============================================================================#
nc = 6   # estados: X, N, G, F, E, O2

# Vector gradiente del objetivo interno (max biomasa)
d  = zeros(nv); d[obj] = -1.0

# Peso de parsimonia (L2) en el pFBA
w  = 1e-20

# Pesos de complementariedad
phi1 = 1e0  # FO_L
phi2 = 1e0  # FO_upt
phi3 = 1e0  # FO_U
phi4 = 1e0  # FO_prod (etanol + biomasa)

# Escalado de concentraciones: c_physical = c_model ÃƒÆ’Ã¢â‚¬â€ cs
# Normaliza concentraciones a O(1) para mejor conditioning
cs = [1.0,      # X: biomasa (gDW/L) ~0.5-2
      0.2,      # N: nitrÃƒÆ’Ã‚Â³geno (gN/L) ~0.1-0.3
      100.0,    # G: glucosa (g/L) ~0-120
      100.0,    # F: fructosa (g/L) ~0-120
      10.0,     # E: etanol (g/L) ~0-60
      0.01]     # O2: oxigeno disuelto (g/L) ~0-0.01

# Escalado de flujos: v_physical = v_model ÃƒÆ’Ã¢â‚¬â€ vs[mc]
# Normaliza flujos grandes para reducir spread en Jacobiano
FLUX_SCALE_TARGET = 50.0  # MÃƒÆ’Ã‚Â¡s conservador que 10.0
vs = ones(nv)
let n_scaled = 0
  for mc in 1:nv
    bound_range = max(abs(vlb[mc]), abs(vub[mc]))
    if bound_range > FLUX_SCALE_TARGET
      vs[mc] = bound_range / FLUX_SCALE_TARGET
      n_scaled += 1
    end
  end
  println("[SCALE] Flujos escalados: $(n_scaled)/$(nv) (target=$(FLUX_SCALE_TARGET))")
end

#=============================================================================
  10. CONDICIONES INICIALES Y DISCRETIZACIÃƒÆ’Ã¢â‚¬Å“N
=============================================================================#
# Condiciones iniciales: X(gDW/L), N(gN/L), G(g/L), F(g/L), E(g/L), O2(g/L)
X0 = 0.5
N0 = 0.14
G0 = 110.0
F0 = 110.0
E0 = 0.0
O2_0 = 0.0

c0 = [X0, N0, G0, F0, E0, O2_0]

# ParÃƒÆ’Ã‚Â¡metros de discretizaciÃƒÆ’Ã‚Â³n (reducido para exploraciÃƒÆ’Ã‚Â³n rÃƒÆ’Ã‚Â¡pida)
nfe    = 6           # nÃƒÆ’Ã‚Âºmero de elementos finitos (usar 12+ para producciÃƒÆ’Ã‚Â³n)
ncp    = 3           # puntos de colocaciÃƒÆ’Ã‚Â³n (Radau IIA)
th     = 120.0       # horizonte temporal (horas)
h      = th / nfe    # longitud nominal del elemento finito
hm     = fill(h, nfe)
var_h  = 0.3         # variaciÃƒÆ’Ã‚Â³n permitida en hv (Ãƒâ€šÃ‚Â±50%)

println("\n>>> DiscretizaciÃƒÆ’Ã‚Â³n: nfe=$(nfe), ncp=$(ncp), th=$(th)h")
println(">>> Condiciones iniciales: X0=$(X0), N0=$(N0), G0=$(G0), F0=$(F0), E0=$(E0), O2_0=$(O2_0)")
println(">>> Modo anaerobico estricto: v_o2 fijado en 0 y O2 inicial en 0")
println(">>> GEM final: $(nm) metabolitos, $(nv) reacciones\n")

#=============================================================================
  11. RESOLVER MPCC
=============================================================================#
println(">>> Resolviendo MPCC (pFBA + KKT embebido)...")
sol, solv, solcd, soll, solalL, solalU, solag, solap, solh, diag = pFBA_KKT_flux_Zenteno_O2minimal(
    c0;
    eps_flux = EPS_FLUX,
)
if haskey(diag, :solver)
    println(">>> Solver status: term=$(diag.solver.term), primal=$(diag.solver.primal), dual=$(diag.solver.dual)")
    @printf(">>> Max violaciones internas: uptake=%.3e, product=%.3e, max|v_o2|=%.3e\n",
            diag.max_uptake_violation, diag.max_product_violation, diag.o2_abs_max)
    if diag.solver.term != JuMP.MOI.OPTIMAL && diag.solver.term != JuMP.MOI.LOCALLY_SOLVED
        @warn "La solucion MPCC no es OPTIMAL/LOCALLY_SOLVED; revisar saturaciones con cautela."
    end
end

#=============================================================================
  12. CONSTRUIR VECTORES DE TIEMPO Y RESULTADOS
=============================================================================#
const radau_pts = [0.15505, 0.64495, 1.00000]

# Tiempos de inicio de cada FE
ts = Vector{Float64}(undef, nfe + 1)
ts[1] = 0.0
for i in 2:nfe+1
    ts[i] = ts[i-1] + solh[i-1]
end

# Grilla completa de tiempo (nfe*ncp + 1 puntos: inicio + todos los puntos de colocaciÃƒÆ’Ã‚Â³n)
n_pts = nfe * ncp + 1
tsn = Vector{Float64}(undef, n_pts)
xk  = Matrix{Float64}(undef, nc, n_pts)

# Punto inicial
tsn[1]  = 0.0
xk[:,1] = [X0, N0, G0, F0, E0, O2_0]

# Puntos de colocaciÃƒÆ’Ã‚Â³n
for i in 1:nfe
    for j in 1:ncp
        kk = (i-1) * ncp + j + 1
        xk[:,kk] = sol[:,i,j]
        tsn[kk]  = ts[i] + radau_pts[j] * solh[i]
    end
end

#=============================================================================
  13. ESCRIBIR RESULTADOS
=============================================================================#
writedlm("xk_Zenteno.csv", xk)
writedlm("tsn_Zenteno.csv", tsn)
println("\n>>> Resultados escritos: xk_Zenteno.csv, tsn_Zenteno.csv")

#=============================================================================
  14. RESUMEN DE RESULTADOS
=============================================================================#
println("\n" * "="^60)
println("   RESUMEN DE RESULTADOS")
println("="^60)
for (s, name) in enumerate(("X [gDW/L]", "N [gN/L]", "G [g/L]", "F [g/L]", "E [g/L]", "O2 [g/L]"))
    vals = xk[s, :]
    @printf("   %-12s: inicio=%8.4f  fin=%8.4f  min=%8.4f  max=%8.4f\n",
            name, vals[1], vals[end], minimum(vals), maximum(vals))
end
println("   Horizonte  : $(tsn[1]) -> $(tsn[end]) h")
println("   hv range   : $(minimum(solh)) -> $(maximum(solh)) h")

# Flujos clave en el ultimo FE
println("\n   Flujos clave (ultimo FE):")
@printf("     v_obj (biomasa)  = %10.6f 1/h\n",        solv[obj, nfe])
@printf("     v_glu (glucosa)  = %10.4f mmol/gDW/h\n",  solv[glu, nfe])
@printf("     v_fru (fructosa) = %10.4f mmol/gDW/h\n",  solv[fru, nfe])
@printf("     v_eth (etanol)   = %10.4f mmol/gDW/h\n",  solv[eth, nfe])
v_o2_last = abs(solv[o2, nfe]) <= O2_ZERO_TOL ? 0.0 : solv[o2, nfe]
@printf("     v_o2  (oxigeno)  = %10.4f mmol/gDW/h\n",  v_o2_last)
println("="^60)

#=============================================================================
  15. GRAFICOS: CONCENTRACION VS FLUJO/RESTRICCION
=============================================================================#

using Plots

tfe = ts[2:end]
idx_glu_upt = findfirst(==(glu), UPTAKE_IDXS)
idx_fru_upt = findfirst(==(fru), UPTAKE_IDXS)
idx_obj_prod = findfirst(==(obj), PRODUCT_IDXS)
idx_eth_prod = findfirst(==(eth), PRODUCT_IDXS)
@assert idx_glu_upt !== nothing "No se encontro glu en UPTAKE_IDXS"
@assert idx_fru_upt !== nothing "No se encontro fru en UPTAKE_IDXS"
@assert idx_obj_prod !== nothing "No se encontro obj en PRODUCT_IDXS"
@assert idx_eth_prod !== nothing "No se encontro eth en PRODUCT_IDXS"
idx_glu_upt = something(idx_glu_upt)
idx_fru_upt = something(idx_fru_upt)
idx_obj_prod = something(idx_obj_prod)
idx_eth_prod = something(idx_eth_prod)

vX_fe = max.(0.0, solv[obj, :])                              # 1/h
limX_fe = max.(0.0, diag.L_prod[idx_obj_prod, :])            # 1/h

vG_fe = max.(0.0, -solv[glu, :])                             # mmol/gDW/h
limG_fe = max.(0.0, diag.L_upt[idx_glu_upt, :])              # mmol/gDW/h

vF_fe = max.(0.0, -solv[fru, :])                             # mmol/gDW/h
limF_fe = max.(0.0, diag.L_upt[idx_fru_upt, :])              # mmol/gDW/h

vE_fe = max.(0.0, solv[eth, :])                              # mmol/gDW/h
limE_fe = max.(0.0, diag.L_prod[idx_eth_prod, :])            # mmol/gDW/h
vO_fe = [abs(solv[o2, i]) <= O2_ZERO_TOL ? 0.0 : max(0.0, -solv[o2, i]) for i in 1:nfe]  # mmol/gDW/h

vN_fe = zeros(nfe)                                            # gN/gDW/h
limN_fe = zeros(nfe)                                          # gN/gDW/h
for i in 1:nfe
    for k in 1:n_up
        if IS_NIT[k] > 0.5
            idx = UPTAKE_IDXS[k]
            vN_fe[i] += max(0.0, -solv[idx, i]) * N_atoms_vec[idx] * MW_N
            limN_fe[i] += max(0.0, diag.L_upt[k, i]) * N_atoms_vec[idx] * MW_N
        end
    end
end

safe_ratio(v, l) = v ./ max.(l, 1e-12)
mean_ratio(v, l) = sum(safe_ratio(v, l)) / length(v)
mean_slack(v, l) = sum(l .- v) / length(v)
max_violation(v, l) = maximum(v .- l)
@printf("   Saturacion promedio v/lim: X=%.3f N=%.3f G=%.3f F=%.3f E=%.3f\n",
        mean_ratio(vX_fe, limX_fe),
        mean_ratio(vN_fe, limN_fe),
        mean_ratio(vG_fe, limG_fe),
        mean_ratio(vF_fe, limF_fe),
        mean_ratio(vE_fe, limE_fe))
@printf("   Slack promedio lim-v:      X=%.3e N=%.3e G=%.3e F=%.3e E=%.3e\n",
        mean_slack(vX_fe, limX_fe),
        mean_slack(vN_fe, limN_fe),
        mean_slack(vG_fe, limG_fe),
        mean_slack(vF_fe, limF_fe),
        mean_slack(vE_fe, limE_fe))
@printf("   Max violacion (v-lim):     X=%.3e N=%.3e G=%.3e F=%.3e E=%.3e\n",
        max_violation(vX_fe, limX_fe),
        max_violation(vN_fe, limN_fe),
        max_violation(vG_fe, limG_fe),
        max_violation(vF_fe, limF_fe),
        max_violation(vE_fe, limE_fe))
println("   O2: modo anaerobico estricto (v_o2 = 0, O2 inicial = 0)")

# Columna 1: concentraciones
pXc = plot(tsn, xk[1,:], label="X (gDW/L)", lw=2, color=:blue, xlabel="Tiempo (h)", ylabel="Concentracion")
pNc = plot(tsn, xk[2,:], label="N (gN/L)", lw=2, color=:green, xlabel="Tiempo (h)", ylabel="Concentracion")
pGc = plot(tsn, xk[3,:], label="G (g/L)", lw=2, color=:red, xlabel="Tiempo (h)", ylabel="Concentracion")
pFc = plot(tsn, xk[4,:], label="F (g/L)", lw=2, color=:orange, xlabel="Tiempo (h)", ylabel="Concentracion")
pEc = plot(tsn, xk[5,:], label="E (g/L)", lw=2, color=:purple, xlabel="Tiempo (h)", ylabel="Concentracion")
pOc = plot(tsn, xk[6,:], label="O2 (g/L)", lw=2, color=:cyan, xlabel="Tiempo (h)", ylabel="Concentracion")

# Columna 2: flujo v vs restriccion cinetica
pXv = plot(tfe, vX_fe, label="v_obj (1/h)", lw=2, color=:blue, xlabel="Tiempo (h)", ylabel="Flujo/limite")
plot!(pXv, tfe, limX_fe, label="lim_obj (1/h)", lw=2, color=:black, ls=:dash)

pNv = plot(tfe, vN_fe, label="vN efectivo (gN/gDW/h)", lw=2, color=:green, xlabel="Tiempo (h)", ylabel="Flujo/limite")
plot!(pNv, tfe, limN_fe, label="limN (gN/gDW/h)", lw=2, color=:black, ls=:dash)

pGv = plot(tfe, vG_fe, label="-v_glu (mmol/gDW/h)", lw=2, color=:red, xlabel="Tiempo (h)", ylabel="Flujo/limite")
plot!(pGv, tfe, limG_fe, label="lim_glu (mmol/gDW/h)", lw=2, color=:black, ls=:dash)

pFv = plot(tfe, vF_fe, label="-v_fru (mmol/gDW/h)", lw=2, color=:orange, xlabel="Tiempo (h)", ylabel="Flujo/limite")
plot!(pFv, tfe, limF_fe, label="lim_fru (mmol/gDW/h)", lw=2, color=:black, ls=:dash)

pEv = plot(tfe, vE_fe, label="v_eth (mmol/gDW/h)", lw=2, color=:purple, xlabel="Tiempo (h)", ylabel="Flujo/limite")
plot!(pEv, tfe, limE_fe, label="lim_eth (mmol/gDW/h)", lw=2, color=:black, ls=:dash)

pOv = plot(tfe, vO_fe, label="-v_o2 (mmol/gDW/h)", lw=2, color=:cyan, xlabel="Tiempo (h)", ylabel="Flujo")
plot!(pOv, tfe, zeros(nfe), label="lim_o2 = 0 (anaerobico)", lw=1.5, color=:black, ls=:dot)

plt = plot(pXc, pXv, pNc, pNv, pGc, pGv, pFc, pFv, pEc, pEv, pOc, pOv,
           layout=(6,2), size=(1400,1300), title="dFBA Zenteno (yeast-GEM $(GEM_VERSION))")
savefig(plt, "zenteno_dfba_result.png")
println(">>> Grafico guardado: zenteno_dfba_result.png")
