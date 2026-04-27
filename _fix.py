import json

f = open('test.ipynb', 'r', encoding='utf-8')
nb = json.load(f)
f.close()

c_fs = nb['cells'][66]
src = c_fs['source']

# Fix line 422: bare tuple -> checks.append(...)
old = src[422]
print(f"OLD L422: {old!r}")

# Replace the bare tuple with checks.append(...)
# Remove trailing comma and newline, wrap in checks.append()
indent = old[:len(old) - len(old.lstrip())]
src[422] = indent + 'checks.append(("_run_dfba_colloc_multiphase disponible", "_run_dfba_colloc_multiphase" in globals()))\n'
print(f"NEW L422: {src[422]!r}")

nb['cells'][66]['source'] = src

f = open('test.ipynb', 'w', encoding='utf-8')
json.dump(nb, f, ensure_ascii=False)
f.close()

print("Fix applied and saved.")
