import json
f = open('test.ipynb', 'r', encoding='utf-8')
nb = json.load(f)
f.close()

for i, c in enumerate(nb['cells']):
    src = ''.join(c['source'])
    if 'def forward_simulate' in src:
        print(f"forward_simulate at index {i}, lines={len(c['source'])}")
        print(f"  Has multiphase: {'_run_dfba_colloc_multiphase' in src}")
        print(f"  Has balanced: {'_run_dfba_colloc_balanced' in src}")
    if 'DEMO' in src and 'BIOLOGICAL' in src:
        print(f"Demo cell at index {i}, lines={len(c['source'])}")
    if 'isoamyl_acetate' in src.lower():
        print(f"isoamyl ref at index {i}")

for idx in [67, 68, 69, 70]:
    if idx < len(nb['cells']):
        c = nb['cells'][idx]
        first = c['source'][0][:60] if c['source'] else 'EMPTY'
        print(f"\nIndex {idx}: type={c['cell_type']}, lines={len(c['source'])}")
        print(f"  First: {repr(first)}")
