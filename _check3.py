import json
f = open('test.ipynb', 'r', encoding='utf-8')
nb = json.load(f)
f.close()

c_fs = nb['cells'][66]
src = c_fs['source']

# Find the dispatch section
for i, line in enumerate(src):
    if '_run_dfba_colloc_multiphase' in line or '_run_dfba_colloc_balanced' in line:
        # Print surrounding context
        start = max(0, i-1)
        end = min(len(src), i+8)
        for j in range(start, end):
            print(f"  L{j:3d}: {src[j]}", end='')
        print("  ---")
