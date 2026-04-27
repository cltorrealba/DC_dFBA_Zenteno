import json

f = open('test.ipynb', 'r', encoding='utf-8')
nb = json.load(f)
f.close()

c_fs = nb['cells'][66]
src = c_fs['source']

# Print lines around the diagnose_forward checks (around L420-430)
for i in range(415, min(len(src), 440)):
    if 'disponible' in src[i] or 'checks' in src[i]:
        print(f"  L{i:3d}: {src[i]}", end='')

print("\n--- Now fixing ---")

# Find and fix the bare tuple line
for i, line in enumerate(src):
    if '_run_dfba_colloc_multiphase disponible' in line and 'checks.append' not in line:
        # Check the pattern of surrounding lines for correct format
        # Print context
        for j in range(max(0,i-2), min(len(src), i+3)):
            print(f"  L{j}: {src[j]}", end='')
        
        # Get the indentation
        stripped = line.lstrip()
        indent = line[:len(line) - len(stripped)]
        
        # Check if this is inside a checks = [...] list literal or checks.append()
        # Look backwards for the pattern
        for j in range(i-1, max(0, i-10), -1):
            if 'checks' in src[j]:
                print(f"\n  Pattern context at L{j}: {src[j]}", end='')
                break
        
        break

print("\n\nDone analysis")
