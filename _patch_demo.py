import json

f = open('test.ipynb', 'r', encoding='utf-8')
nb = json.load(f)
f.close()

demo_idx = None
for i, c in enumerate(nb['cells']):
    src = ''.join(c['source'])
    if 'DEMO' in src and 'Validacion biologica del forward multifase' in src:
        demo_idx = i
        print(f"Found DEMO cell at index {i}, lines={len(c['source'])}")
        break

if demo_idx is None:
    print("ERROR: DEMO cell not found")
    exit(1)

new_code = '''# ======================================================================================
# DEMO  Validacion biologica del forward multifase (Step 1)
# ======================================================================================
# Ejecutar las celdas anteriores primero (modelo, cinetica, solver, PARAM_SPECS, forward)

print("="*70)
print("DEMO: Forward multifase con parametros nominales")
print("="*70)

# Run forward con theta nominal (forward_simulate aplica apply_theta internamente)
res = forward_simulate(theta_nominal_unit, verbose=True)

if res.get("status") != "ok":
    print(f"\\n[DEMO] ERROR: forward_simulate fallo: {res.get('error', '?')}")
else:
    sdf = res["states_df"]
    fdf = res["flux_df"]
    t_h  = sdf["t_h"].to_numpy()
    X    = sdf["X_gDW_L"].to_numpy()
    G    = sdf["G_g_L"].to_numpy()
    F    = sdf["F_g_L"].to_numpy()
    E    = sdf["E_g_L"].to_numpy()
    Nf   = sdf["N_free_gN_L"].to_numpy()
    Prot = sdf["Prot_g_L"].to_numpy()
    Carb = sdf["Carb_g_L"].to_numpy()

    # Aromas from states_df columns (format: {rid}_mg_L)
    aroma_cols = [c for c in sdf.columns if c.endswith("_mg_L") and not c.startswith("AA_")]
    ia_col = f"{ISOAMYL_ACETATE_RID}_mg_L"

    # Phases from flux_df
    phases = fdf["phase"].tolist() if "phase" in fdf.columns else []

    fig, axes = plt.subplots(3, 2, figsize=(14, 12), sharex=True)

    # (0,0) Biomasa
    axes[0,0].plot(t_h, X, 'k-', lw=2)
    axes[0,0].set_ylabel("Biomasa X (g/L)")
    axes[0,0].set_title("Biomasa")

    # (0,1) Azucares + Etanol
    axes[0,1].plot(t_h, G, 'b-', lw=1.5, label="Glucosa")
    axes[0,1].plot(t_h, F, 'c-', lw=1.5, label="Fructosa")
    axes[0,1].plot(t_h, E, 'r-', lw=1.5, label="Etanol")
    axes[0,1].set_ylabel("Conc (g/L)")
    axes[0,1].legend(fontsize=8)
    axes[0,1].set_title("Sustratos / Etanol")

    # (1,0) Nitrogeno
    axes[1,0].plot(t_h, Nf, 'g-', lw=2)
    axes[1,0].axhline(0.02, color='gray', ls='--', lw=0.8, label="N_LIMITED_THRESHOLD")
    axes[1,0].set_ylabel("N libre (g/L)")
    axes[1,0].legend(fontsize=8)
    axes[1,0].set_title("Nitrogeno asimilable")

    # (1,1) Composicion celular
    axes[1,1].plot(t_h, Prot, 'm-', lw=1.5, label="Prot content")
    axes[1,1].plot(t_h, Carb, 'y-', lw=1.5, label="Carb content")
    axes[1,1].set_ylabel("Fraccion masa (g/L)")
    axes[1,1].legend(fontsize=8)
    axes[1,1].set_title("Composicion celular")

    # (2,0) Isoamyl acetate
    if ia_col in sdf.columns:
        ia_mg = sdf[ia_col].to_numpy()
        axes[2,0].plot(t_h, ia_mg, 'darkorange', lw=2)
        axes[2,0].set_ylabel("Isoamyl acetate (mg/L)")
        axes[2,0].set_title(f"Isoamyl acetate ({ISOAMYL_ACETATE_RID})")
    else:
        axes[2,0].text(0.5, 0.5, f"Col {ia_col} not found\\nAvailable: {aroma_cols[:5]}",
                        ha='center', va='center', transform=axes[2,0].transAxes, fontsize=8)

    # (2,1) All aromas
    for ac in aroma_cols:
        axes[2,1].plot(t_h, sdf[ac].to_numpy(), lw=1, label=ac.replace("_mg_L",""))
    axes[2,1].set_ylabel("Aroma (mg/L)")
    axes[2,1].legend(fontsize=6, ncol=2)
    axes[2,1].set_title("Todos los aromas")

    for ax in axes[-1, :]:
        ax.set_xlabel("Tiempo (h)")

    fig.suptitle("Forward multifase - Validacion biologica Step 1", fontsize=14, y=1.01)
    fig.tight_layout()
    plt.show()

    # Phase counts from flux_df
    if phases:
        from collections import Counter
        pc = Counter(phases)
        print(f"\\n[DEMO] Phase counts: {dict(pc)}")

    # Key biological checks
    print("\\n" + "="*70)
    print("BIOLOGICAL CHECKS:")
    print(f"  Biomasa final:        {X[-1]:.3f} g/L")
    print(f"  Glucosa final:        {G[-1]:.4f} g/L")
    print(f"  Fructosa final:       {F[-1]:.4f} g/L")
    print(f"  Etanol final:         {E[-1]:.2f} g/L")
    print(f"  N libre final:        {Nf[-1]:.4f} g/L")
    if ia_col in sdf.columns:
        ia_mg = sdf[ia_col].to_numpy()
        ia_max = float(ia_mg.max())
        ia_final = float(ia_mg[-1])
        n_depl_idx = int(np.argmax(Nf < 0.02)) if np.any(Nf < 0.02) else len(Nf)-1
        ia_at_n_depl = float(ia_mg[:n_depl_idx+1].max())
        print(f"  Isoamyl acetate max:  {ia_max:.4f} mg/L")
        print(f"  Isoamyl acetate final:{ia_final:.4f} mg/L")
        print(f"  Isoamyl acetate at N depletion: {ia_at_n_depl:.4f} mg/L")
        ia_check = 'YES' if ia_at_n_depl > 1e-4 else 'NO'
        print(f"  >> Produced during growth: {ia_check}")
    print("="*70)
'''

new_lines = []
for line in new_code.split('\n'):
    new_lines.append(line + '\n')
if new_lines and new_lines[-1] == '\n':
    new_lines[-1] = ''

nb['cells'][demo_idx]['source'] = new_lines

with open('test.ipynb', 'w', encoding='utf-8') as f:
    json.dump(nb, f, ensure_ascii=False, indent=1)

print(f"Patched DEMO cell at index {demo_idx} ({len(new_lines)} lines)")
print("Saved test.ipynb")
