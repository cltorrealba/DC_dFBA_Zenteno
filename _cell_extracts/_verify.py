import json
NB = r'c:\Users\ctorrealba\OneDrive - Viña Concha y Toro S.A\Documentos\Doctorado\Artículos\Artículo - Estimación_dFBA\DC_dFBA_Zenteno\01_preparacion_modelo_cobrapy.ipynb'
with open(NB, 'r', encoding='utf-8') as f:
    nb = json.load(f)

print(f"Total cells: {len(nb['cells'])}")

# Check new cell at index 46
for i in [45, 46, 47]:
    c = nb['cells'][i]
    src = ''.join(c['source'])
    fl = src.strip().split('\n')[0][:100]
    nlines = src.count('\n') + 1
    print(f"Cell {i} (id={c.get('id','?')}): lines={nlines}  {fl}")

print()

# Check PARAM_SPECS
for c in nb['cells']:
    if c.get('id') == 'd6fbabd2':
        src = ''.join(c['source'])
        for kw in ['ALPHA_GROWTH', 'W_ATPM_GROWTH', 'W_ISOAMYLOL_GROWTH', 'W_ISOAMYL_ACETATE_GROWTH', 'growth_overflow']:
            print(f"  PARAM_SPECS has '{kw}': {kw in src}")
        break

print()

# Check forward_simulate
for c in nb['cells']:
    if c.get('id') == '4527a2cf':
        src = ''.join(c['source'])
        print(f"  forward_simulate has '_run_dfba_colloc_multiphase': {'_run_dfba_colloc_multiphase' in src}")
        print(f"  forward_simulate has diagnose check: {'_run_dfba_colloc_multiphase disponible' in src}")
        # Show the solver selection block
        lines = src.split('\n')
        for li, line in enumerate(lines):
            if 'Solver multifase' in line or '_run_dfba_colloc_multiphase' in line:
                print(f"  line {li}: {line.rstrip()[:100]}")
        break
