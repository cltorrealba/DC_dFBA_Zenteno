import json

f = open('test.ipynb', 'r', encoding='utf-8')
nb = json.load(f)
f.close()

c_fs = nb['cells'][66]
src = c_fs['source']

# Fix L422: bare tuple -> checks.append(...)
fixed = False
for i, line in enumerate(src):
    if '_run_dfba_colloc_multiphase disponible' in line and 'checks.append' not in line:
        indent = '    '
        src[i] = indent + 'checks.append(("_run_dfba_colloc_multiphase disponible", "_run_dfba_colloc_multiphase" in globals()))\n'
        fixed = True
        print(f"Fixed line {i}: {src[i]}", end='')
        break

if not fixed:
    print("Line already fixed or not found")

c_fs['source'] = src
nb['cells'][66] = c_fs

with open('test.ipynb', 'w', encoding='utf-8') as f:
    json.dump(nb, f, ensure_ascii=False, indent=1)

print("Saved test.ipynb")
