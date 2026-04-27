import json
f = open('test.ipynb', 'r', encoding='utf-8')
nb = json.load(f)
f.close()

print(f"Total cells: {len(nb['cells'])}")

c_ps = nb['cells'][64]
src_ps = ''.join(c_ps['source'])
print(f"\n[CHECK 1] PARAM_SPECS cell (index 64):")
print(f"  Has ALPHA_GROWTH: {'ALPHA_GROWTH' in src_ps}")
print(f"  Has W_ISOAMYL_ACETATE_GROWTH: {'W_ISOAMYL_ACETATE_GROWTH' in src_ps}")
print(f"  Has LAG_END_H: {'LAG_END_H' in src_ps}")
print(f"  Has DECAY_START_H: {'DECAY_START_H' in src_ps}")
print(f"  Has Grupo 4: {'Grupo 4' in src_ps}")
print(f"  Has Grupo 5: {'Grupo 5' in src_ps}")

c_mp = nb['cells'][46]
src_mp = ''.join(c_mp['source'])
print(f"\n[CHECK 2] Multiphase solver cell (index 46):")
print(f"  Lines: {len(c_mp['source'])}")
print(f"  Has _run_dfba_colloc_multiphase: {'def _run_dfba_colloc_multiphase' in src_mp}")
print(f"  Has _fba_at_state_colloc_multiphase: {'def _fba_at_state_colloc_multiphase' in src_mp}")
print(f"  Has _infer_phase: {'def _infer_phase' in src_mp}")
print(f"  First line: {c_mp['source'][0][:60]!r}")

c_fs = nb['cells'][67]
src_fs = ''.join(c_fs['source'])
print(f"\n[CHECK 3] forward_simulate cell (index 67):")
print(f"  Lines: {len(c_fs['source'])}")
print(f"  Has forward_simulate: {'def forward_simulate' in src_fs}")
print(f"  Has _run_dfba_colloc_multiphase dispatch: {'_run_dfba_colloc_multiphase' in src_fs}")
print(f"  Has elif _run_dfba_colloc_balanced: {'elif' in src_fs and '_run_dfba_colloc_balanced' in src_fs}")
print(f"  First line: {c_fs['source'][0][:60]!r}")

c_demo = nb['cells'][68]
src_demo = ''.join(c_demo['source'])
print(f"\n[CHECK 4] Demo cell (index 68):")
print(f"  Lines: {len(c_demo['source'])}")
print(f"  Has DEMO: {'DEMO' in src_demo}")
print(f"  Has BIOLOGICAL CHECKS: {'BIOLOGICAL CHECKS' in src_demo}")
print(f"  Has isoamyl_acetate: {'isoamyl_acetate' in src_demo.lower() or 'Isoamyl acetate' in src_demo}")
print(f"  First line: {c_demo['source'][0][:60]!r}")
