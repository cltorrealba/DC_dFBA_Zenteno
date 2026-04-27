# Activar leucina como fuente N cinética (alineado con el manejo de otros aminoácidos)
print("\n" + "=" * 100)
print("Ajuste: leucina disponible como nutriente cinético")
print("=" * 100)

LEU_ID = "r_1899"

# 1) Incluir leucina en la lista cinética
if "KINETIC_N_SOURCE_IDS" not in globals():
    raise RuntimeError("KINETIC_N_SOURCE_IDS no está definido en sesión.")

if LEU_ID not in KINETIC_N_SOURCE_IDS:
    KINETIC_N_SOURCE_IDS.append(LEU_ID)

# 2) Asegurar composición de N para leucina
N_atoms_map = dict(globals().get("N_atoms_map", {}))
N_atoms_map[LEU_ID] = 1.0  # leucina tiene 1 átomo de N
globals()["N_atoms_map"] = N_atoms_map

N_profile_map = dict(globals().get("N_profile_map", {}))

# Dar a leucina una disponibilidad tipo aminoácido (0.05),
# compensando desde NH4 para mantener el total ~1.0 en el perfil.
target_leu_profile = 0.05
if N_profile_map:
    if LEU_ID not in N_profile_map:
        nh4_id = "r_1654"
        if nh4_id in N_profile_map:
            N_profile_map[nh4_id] = max(0.0, float(N_profile_map[nh4_id]) - target_leu_profile)
        N_profile_map[LEU_ID] = target_leu_profile

    # Completar faltantes y normalizar solo sobre fuentes cinéticas
    for rid in KINETIC_N_SOURCE_IDS:
        N_profile_map.setdefault(rid, 0.0)

    prof_sum = sum(float(N_profile_map.get(rid, 0.0)) for rid in KINETIC_N_SOURCE_IDS)
    if prof_sum > 0:
        for rid in KINETIC_N_SOURCE_IDS:
            N_profile_map[rid] = float(N_profile_map.get(rid, 0.0)) / prof_sum
else:
    # Fallback: distribución uniforme si no existía perfil previo
    uniform = 1.0 / max(1, len(KINETIC_N_SOURCE_IDS))
    N_profile_map = {rid: uniform for rid in KINETIC_N_SOURCE_IDS}

globals()["N_profile_map"] = N_profile_map

# 3) Recalcular fracciones de N por fuente y límites cinéticos
if "MW_N" in globals():
    N_frac_map = {}
    for rid in KINETIC_N_SOURCE_IDS:
        n_atoms = float(N_atoms_map.get(rid, 1.0))
        n_prof = float(N_profile_map.get(rid, 0.0))
        N_frac_map[rid] = n_prof / max(n_atoms * float(MW_N), 1e-12)
    globals()["N_frac_map"] = N_frac_map

    if "vn" in globals():
        lim_n_map = {rid: float(vn) * float(N_frac_map[rid]) for rid in KINETIC_N_SOURCE_IDS}
        globals()["lim_n_map"] = lim_n_map

# 4) Aplicar bounds de leucina como uptake
if "rxn_by_id" in globals() and LEU_ID in rxn_by_id:
    r = rxn_by_id[LEU_ID]
    old_lb, old_ub = float(r.lower_bound), float(r.upper_bound)
    r.lower_bound = max(old_lb, -1000.0)
    r.upper_bound = min(old_ub, 0.0)
    print(f"[OK] {LEU_ID} | {r.name} | ({old_lb:.3f}, {old_ub:.3f}) -> ({r.lower_bound:.3f}, {r.upper_bound:.3f})")
else:
    print(f"[WARN] {LEU_ID} no está disponible en rxn_by_id; revisa mapeo de reacciones.")

print(f"Fuentes N cinéticas activas: {len(KINETIC_N_SOURCE_IDS)}")
print("Incluye leucina:", LEU_ID in KINETIC_N_SOURCE_IDS)
if "N_profile_map" in globals():
    print("Perfil leucina N_profile_map[r_1899] =", N_profile_map.get(LEU_ID, None))
if "lim_n_map" in globals():
    print("Límite cinético leucina lim_n_map[r_1899] =", globals()["lim_n_map"].get(LEU_ID, None))