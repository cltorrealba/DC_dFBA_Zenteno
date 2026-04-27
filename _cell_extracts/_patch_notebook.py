"""Patch notebook: insert multiphase cell, update PARAM_SPECS and forward_simulate."""
import json, uuid, sys

NB_PATH = r'c:\Users\ctorrealba\OneDrive - Viña Concha y Toro S.A\Documentos\Doctorado\Artículos\Artículo - Estimación_dFBA\DC_dFBA_Zenteno\01_preparacion_modelo_cobrapy.ipynb'
NEW_CELL_PATH = r'c:\Users\ctorrealba\OneDrive - Viña Concha y Toro S.A\Documentos\Doctorado\Artículos\Artículo - Estimación_dFBA\DC_dFBA_Zenteno\_cell_extracts\cell_45b_multiphase.py'

with open(NB_PATH, "r", encoding="utf-8") as f:
    nb = json.load(f)
with open(NEW_CELL_PATH, "r", encoding="utf-8") as f:
    new_cell_source = f.read()

# ── 1) INSERT NEW CELL after cell 45 (id=0ccbc00c) ──
insert_after_id = "0ccbc00c"
insert_idx = None
for i, c in enumerate(nb["cells"]):
    if c.get("id") == insert_after_id:
        insert_idx = i + 1
        break
assert insert_idx is not None, f"Cell {insert_after_id} not found"

new_cell = {
    "cell_type": "code",
    "execution_count": None,
    "id": uuid.uuid4().hex[:8],
    "metadata": {},
    "outputs": [],
    "source": new_cell_source.splitlines(keepends=True),
}
nb["cells"].insert(insert_idx, new_cell)
print(f"[OK] Inserted multiphase cell at index {insert_idx} (id={new_cell['id']})")

# ── 2) UPDATE PARAM_SPECS (cell id=d6fbabd2) ──
for c in nb["cells"]:
    if c.get("id") == "d6fbabd2":
        src = "".join(c["source"])
        marker = "    # === Grupo 3: Turnover / N-reciclaje ==="
        assert marker in src, "Turnover marker not found in PARAM_SPECS"
        new_params = (
            "\n"
            '    # === Grupo 2b: Growth overflow ===\n'
            '    ParamSpec("ALPHA_GROWTH",              _g("ALPHA_GROWTH", 0.95),              0.85,  0.999, False, "growth_overflow", True),\n'
            '    ParamSpec("W_ATPM_GROWTH",             _g("W_ATPM_GROWTH", 1.0),             0.01,  10.0,  True,  "growth_overflow", True),\n'
            '    ParamSpec("W_ISOAMYLOL_GROWTH",        _g("W_ISOAMYLOL_GROWTH", 0.5),        0.01,  10.0,  True,  "growth_overflow", True),\n'
            '    ParamSpec("W_ISOAMYL_ACETATE_GROWTH",  _g("W_ISOAMYL_ACETATE_GROWTH", 0.3),  0.01,  10.0,  True,  "growth_overflow", True),\n'
            "\n"
        )
        src = src.replace(marker, new_params + marker)
        c["source"] = src.splitlines(keepends=True)
        print("[OK] Added 4 growth_overflow params to PARAM_SPECS")
        break

# ── 3) UPDATE forward_simulate (cell id=4527a2cf) ──
for c in nb["cells"]:
    if c.get("id") == "4527a2cf":
        src = "".join(c["source"])

        # The block we want to replace (exact whitespace match from the file)
        old_solver_comment = "# Preferimos el solver balanceado (celda 46); fallback al basico si no existe."
        assert old_solver_comment in src, "Old solver comment not found in forward_simulate"

        # Strategy: find the line with the comment, replace from there through the
        # raise RuntimeError line
        lines = src.split("\n")
        start_idx = None
        end_idx = None
        for li, line in enumerate(lines):
            if old_solver_comment in line:
                start_idx = li
            if start_idx is not None and "raise RuntimeError" in line and "solver dFBA" in line:
                end_idx = li
                break

        assert start_idx is not None and end_idx is not None, "Could not find solver block boundaries"

        # Build replacement lines
        indent = "                "
        replacement = [
            indent + "# Solver multifase con overflow; fallback a balanced/basico.",
            indent + 'if "_run_dfba_colloc_multiphase" in globals():',
            indent + "    res = _run_dfba_colloc_multiphase(",
            indent + '        scenario_name="opt_forward",',
            indent + "        aerobic_mode=False,",
            indent + "        o2_init_g_l=0.0,",
            indent + "        nfe=nfe, ncp=ncp, t_end=t_end, n_iter=n_iter,",
            indent + "    )",
            indent + 'elif "_run_dfba_colloc_balanced" in globals():',
            indent + "    res = _run_dfba_colloc_balanced(",
            indent + '        scenario_name="opt_forward",',
            indent + "        aerobic_mode=False,",
            indent + "        o2_init_g_l=0.0,",
            indent + "        nfe=nfe, ncp=ncp, t_end=t_end, n_iter=n_iter,",
            indent + "        growth_lock_frac=1.0, w_atpm_growth=0.0, w_aa_growth=0.0,",
            indent + "        use_pfba_mix=True,",
            indent + "    )",
            indent + 'elif "_run_dfba_colloc" in globals():',
            indent + "    res = _run_dfba_colloc(",
            indent + '        scenario_name="opt_forward",',
            indent + "        aerobic_mode=False, o2_init_g_l=0.0,",
            indent + "        nfe=nfe, ncp=ncp, t_end=t_end, n_iter=n_iter,",
            indent + "    )",
            indent + "else:",
            indent + '    raise RuntimeError("No hay solver dFBA disponible.")',
        ]

        new_lines = lines[:start_idx] + replacement + lines[end_idx + 1:]
        src = "\n".join(new_lines)
        c["source"] = src.splitlines(keepends=True)
        print(f"[OK] Updated forward_simulate solver selection (replaced lines {start_idx}-{end_idx})")

        # Also add multiphase check to diagnose_forward
        old_check = '    checks.append(("_run_dfba_colloc_balanced disponible", "_run_dfba_colloc_balanced" in globals()))'
        new_check = (
            '    checks.append(("_run_dfba_colloc_multiphase disponible", "_run_dfba_colloc_multiphase" in globals()))\n'
            + old_check
        )
        src2 = "".join(c["source"])
        if old_check in src2 and "_run_dfba_colloc_multiphase" not in src2.split("diagnose_forward")[1]:
            src2 = src2.replace(old_check, new_check, 1)
            c["source"] = src2.splitlines(keepends=True)
            print("[OK] Added multiphase check to diagnose_forward")
        break

# ── SAVE ──
with open(NB_PATH, "w", encoding="utf-8") as f:
    json.dump(nb, f, indent=1, ensure_ascii=False)

print(f"\n[DONE] Notebook saved with {len(nb['cells'])} cells")
