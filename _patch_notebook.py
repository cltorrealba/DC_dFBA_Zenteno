"""
Patch test.ipynb to add multiphase architecture (Step 1).

Modifications:
  1. Cell 63 (PARAM_SPECS): Add Grupo 4 (overflow) and Grupo 5 (phase) params
  2. Insert new cell after cell 45: multiphase solver from cell_46b_multiphase.py
  3. Cell 65 → 66 (forward_simulate): Add multiphase dispatch as first preference
  4. Insert demo cell after forward_simulate

Creates a backup test.ipynb.bak before modifying.
"""
import json
import shutil
import os
import uuid

NOTEBOOK = "test.ipynb"
BACKUP   = "test.ipynb.bak"
MULTIPHASE_SRC = os.path.join("_cell_extracts", "cell_46b_multiphase.py")

# ── 0. Backup ──────────────────────────────────────────────────────────────────
shutil.copy2(NOTEBOOK, BACKUP)
print(f"[patch] Backup: {BACKUP}")

with open(NOTEBOOK, "r", encoding="utf-8") as f:
    nb = json.load(f)

cells = nb["cells"]
print(f"[patch] Original cells: {len(cells)}")

# ── 1. Patch PARAM_SPECS (cell 63) ─────────────────────────────────────────────
c63 = cells[63]
src63 = c63["source"]

# Find the line with the closing bracket of PARAM_SPECS list: standalone "]"
# It comes right after K_AA_UPTAKE_GROWTH line
insert_idx = None
for i, line in enumerate(src63):
    if "K_AA_UPTAKE_GROWTH" in line:
        # The "]" is on the next line
        insert_idx = i + 1
        break

if insert_idx is None:
    raise RuntimeError("Could not find K_AA_UPTAKE_GROWTH in cell 63")

# Verify the line at insert_idx is the closing bracket
assert src63[insert_idx].strip() == "]", f"Expected ']' at line {insert_idx}, got: {src63[insert_idx]!r}"

# New lines to insert BEFORE the "]"
new_params = [
    "\n",
    "    # === Grupo 4: Overflow (multifase, celda 46b) ===\n",
    '    ParamSpec("ALPHA_GROWTH",              _g("ALPHA_GROWTH",              0.95),  0.80,  0.999, False, "overflow", False),\n',
    '    ParamSpec("W_ATPM_GROWTH",             _g("W_ATPM_GROWTH",             1.0),   0.1,   5.0,   False, "overflow", False),\n',
    '    ParamSpec("W_ISOAMYLOL_GROWTH",        _g("W_ISOAMYLOL_GROWTH",        0.5),   0.01,  3.0,   True,  "overflow", False),\n',
    '    ParamSpec("W_ISOAMYL_ACETATE_GROWTH",  _g("W_ISOAMYL_ACETATE_GROWTH",  0.3),   0.01,  3.0,   True,  "overflow", False),\n',
    '    ParamSpec("W_ETHYL_ACETATE_GROWTH",    _g("W_ETHYL_ACETATE_GROWTH",    0.2),   0.01,  3.0,   True,  "overflow", False),\n',
    "\n",
    "    # === Grupo 5: Phase transitions (multifase, celda 46b) ===\n",
    '    ParamSpec("LAG_END_H",           _g("LAG_END_H",           4.0),   1.0,   10.0,  False, "phase", False),\n',
    '    ParamSpec("N_LIMITED_THRESHOLD",  _g("N_LIMITED_THRESHOLD",  0.02),  0.005, 0.1,   True,  "phase", False),\n',
    '    ParamSpec("DECAY_START_H",       _g("DECAY_START_H",       60.0),  40.0,  70.0,  False, "phase", False),\n',
]

src63[insert_idx:insert_idx] = new_params
c63["source"] = src63
print(f"[patch] PARAM_SPECS: inserted {len(new_params)} lines (Grupo 4+5)")


# ── 2. Insert multiphase solver cell after cell 45 ────────────────────────────
with open(MULTIPHASE_SRC, "r", encoding="utf-8") as f:
    mp_code = f.read()

# Convert to line-based source (ipynb format)
mp_lines = []
for line in mp_code.split("\n"):
    mp_lines.append(line + "\n")
# Remove trailing empty-newline artifact
if mp_lines and mp_lines[-1] == "\n":
    mp_lines[-1] = ""

new_cell_mp = {
    "cell_type": "code",
    "execution_count": None,
    "id": uuid.uuid4().hex[:8],
    "metadata": {},
    "outputs": [],
    "source": mp_lines,
}

# Insert after cell 45 (balanced solver) → becomes new cell 46
cells.insert(46, new_cell_mp)
print(f"[patch] Inserted multiphase solver cell at index 46 ({len(mp_lines)} lines)")

# After insert: old cell 63 → now 64, old cell 65 → now 66


# ── 3. Patch forward_simulate (now cell 66 after insert) ──────────────────────
# Cell indices shifted by 1 after insertion
c_fs = cells[66]
src_fs = c_fs["source"]

# Find the dispatch pattern: look for "_run_dfba_colloc_balanced" in globals() check
dispatch_idx = None
for i, line in enumerate(src_fs):
    if '_run_dfba_colloc_balanced' in line and 'globals()' in line:
        dispatch_idx = i
        break

if dispatch_idx is None:
    raise RuntimeError("Could not find _run_dfba_colloc_balanced dispatch in forward_simulate")

# Insert multiphase dispatch BEFORE the balanced dispatch
# The balanced block is:
#     if "_run_dfba_colloc_balanced" in globals():
#         res = _run_dfba_colloc_balanced(...)
# We change to:
#     if "_run_dfba_colloc_multiphase" in globals():
#         res = _run_dfba_colloc_multiphase(...)
#     elif "_run_dfba_colloc_balanced" in globals():
#         res = _run_dfba_colloc_balanced(...)
# So we change "if" to "elif" on the balanced line, and add multiphase block before it.

# Get the indentation of the existing "if" line
indent = ""
for ch in src_fs[dispatch_idx]:
    if ch in " \t":
        indent += ch
    else:
        break

# Change existing "if" to "elif"
src_fs[dispatch_idx] = src_fs[dispatch_idx].replace("if ", "elif ", 1)

# Find the res= line for balanced to know what args to copy
# Look for the res = _run_dfba_colloc_balanced(...) block
res_line_idx = None
for j in range(dispatch_idx + 1, min(dispatch_idx + 10, len(src_fs))):
    if "_run_dfba_colloc_balanced(" in src_fs[j]:
        res_line_idx = j
        break

# Extract the call block (may span multiple lines)
# We'll build the multiphase block with the same args pattern
multiphase_block = [
    f'{indent}if "_run_dfba_colloc_multiphase" in globals():\n',
    f'{indent}    res = _run_dfba_colloc_multiphase(\n',
    f'{indent}        scenario_name="opt_forward",\n',
    f'{indent}        aerobic_mode=False, o2_init_g_l=0.0,\n',
    f'{indent}        nfe=nfe, ncp=ncp, t_end=t_end, n_iter=n_iter,\n',
    f'{indent}    )\n',
]

src_fs[dispatch_idx:dispatch_idx] = multiphase_block
c_fs["source"] = src_fs
print(f"[patch] forward_simulate: added multiphase dispatch ({len(multiphase_block)} lines)")

# Also patch diagnose_forward: add multiphase check
# Find diagnose_forward function
diag_idx = None
for i, line in enumerate(src_fs):
    if "def diagnose_forward" in line:
        diag_idx = i
        break

if diag_idx is not None:
    # Find the line that checks _run_dfba_colloc_balanced in diagnose_forward
    for j in range(diag_idx, min(diag_idx + 40, len(src_fs))):
        if "_run_dfba_colloc_balanced" in src_fs[j] and "disponible" in src_fs[j]:
            # Add multiphase check right before this line
            d_indent = ""
            for ch in src_fs[j]:
                if ch in " \t":
                    d_indent += ch
                else:
                    break
            new_diag_line = f'{d_indent}("_run_dfba_colloc_multiphase disponible", "_run_dfba_colloc_multiphase" in globals()),\n'
            src_fs.insert(j, new_diag_line)
            print(f"[patch] diagnose_forward: added multiphase check")
            break
    c_fs["source"] = src_fs


# ── 4. Insert demo cell ──────────────────────────────────────────────────────
demo_code = r'''# ======================================================================================
# DEMO — Validacion biologica del forward multifase (Step 1)
# ======================================================================================
# Ejecutar las celdas anteriores primero (modelo, cinetica, solver, PARAM_SPECS, forward)

print("="*70)
print("DEMO: Forward multifase con parametros nominales")
print("="*70)

# Aplicar parametros nominales
apply_theta(theta_nominal_unit, verbose=True)

# Run forward
res = forward_simulate(verbose=True)

if res is None:
    print("\n[DEMO] ERROR: forward_simulate retorno None")
else:
    t_h  = res["t_h"]
    X    = res["states"][:, IDX_X]
    G    = res["states"][:, IDX_G]
    F    = res["states"][:, IDX_F]
    E    = res["states"][:, IDX_E]
    Nf   = res["states"][:, IDX_NFREE]
    Prot = res["states"][:, IDX_PROT]
    Carb = res["states"][:, IDX_CARB]

    # Aromas
    n_aa    = len(AA_RIDS)
    n_aroma = len(AROMA_RIDS)
    aroma_states = res["states"][:, IDX_AROMA0 : IDX_AROMA0 + n_aroma]

    # Find isoamyl acetate index
    ia_idx = None
    for k, rid in enumerate(AROMA_RIDS):
        if rid == ISOAMYL_ACETATE_RID:
            ia_idx = k
            break

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
    axes[1,1].set_ylabel("Fraccion masa")
    axes[1,1].legend(fontsize=8)
    axes[1,1].set_title("Composicion celular")

    # (2,0) Isoamyl acetate
    if ia_idx is not None:
        axes[2,0].plot(t_h, aroma_states[:, ia_idx]*1e3, 'darkorange', lw=2)
        axes[2,0].set_ylabel("Isoamyl acetate (mg/L)")
        axes[2,0].set_title(f"Isoamyl acetate ({ISOAMYL_ACETATE_RID})")
    else:
        axes[2,0].text(0.5, 0.5, "Isoamyl acetate not found", ha='center', va='center',
                        transform=axes[2,0].transAxes)

    # (2,1) All aromas
    for k in range(n_aroma):
        axes[2,1].plot(t_h, aroma_states[:, k]*1e3, lw=1, label=AROMA_RIDS[k])
    axes[2,1].set_ylabel("Aroma (mg/L)")
    axes[2,1].legend(fontsize=6, ncol=2)
    axes[2,1].set_title("Todos los aromas")

    for ax in axes[-1, :]:
        ax.set_xlabel("Tiempo (h)")

    fig.suptitle("Forward multifase - Validacion biologica Step 1", fontsize=14, y=1.01)
    fig.tight_layout()
    plt.show()

    # Phase annotation
    if "phase_log" in res:
        print("\n[DEMO] Phase log (first 10 entries):")
        for entry in res["phase_log"][:10]:
            print(f"  t={entry.get('t',0):.1f}h  phase={entry.get('phase','?')}")

    # Key biological checks
    print("\n" + "="*70)
    print("BIOLOGICAL CHECKS:")
    print(f"  Biomasa final:        {X[-1]:.3f} g/L")
    print(f"  Glucosa final:        {G[-1]:.4f} g/L")
    print(f"  Fructosa final:       {F[-1]:.4f} g/L")
    print(f"  Etanol final:         {E[-1]:.2f} g/L")
    print(f"  N libre final:        {Nf[-1]:.4f} g/L")
    if ia_idx is not None:
        ia_max = aroma_states[:, ia_idx].max() * 1e3
        ia_final = aroma_states[:, ia_idx][-1] * 1e3
        # Check: is isoamyl acetate produced DURING growth (before N depletion)?
        n_depl_idx = np.argmax(Nf < 0.02) if np.any(Nf < 0.02) else len(Nf)-1
        ia_at_n_depl = aroma_states[:n_depl_idx+1, ia_idx].max() * 1e3
        print(f"  Isoamyl acetate max:  {ia_max:.4f} mg/L")
        print(f"  Isoamyl acetate final:{ia_final:.4f} mg/L")
        print(f"  Isoamyl acetate at N depletion: {ia_at_n_depl:.4f} mg/L")
        print(f"  >> Produced during growth: {'YES' if ia_at_n_depl > 1e-4 else 'NO'}")
    print("="*70)
'''

demo_lines = []
for line in demo_code.split("\n"):
    demo_lines.append(line + "\n")
if demo_lines and demo_lines[-1] == "\n":
    demo_lines[-1] = ""

new_cell_demo = {
    "cell_type": "code",
    "execution_count": None,
    "id": uuid.uuid4().hex[:8],
    "metadata": {},
    "outputs": [],
    "source": demo_lines,
}

# Insert after forward_simulate cell (now at index 66, so insert at 67)
# But we also need to insert after diagnose_forward+ester+synthetic data which are all in cell 66
cells.insert(67, new_cell_demo)
print(f"[patch] Inserted demo cell at index 67 ({len(demo_lines)} lines)")

# ── 5. Write ──────────────────────────────────────────────────────────────────
nb["cells"] = cells
with open(NOTEBOOK, "w", encoding="utf-8") as f:
    json.dump(nb, f, ensure_ascii=False, indent=1)

print(f"\n[patch] Done. Total cells: {len(cells)}")
print("[patch] Changes:")
print("  - Cell 63 (cell 64 in VSCode): PARAM_SPECS + Grupo 4 (overflow) + Grupo 5 (phase)")
print("  - Cell 46 (new): Multiphase solver from cell_46b_multiphase.py")
print("  - Cell 66 (cell 67 in VSCode): forward_simulate with multiphase dispatch first")
print("  - Cell 67 (new): Demo/validation cell")
print(f"[patch] Backup saved as {BACKUP}")
