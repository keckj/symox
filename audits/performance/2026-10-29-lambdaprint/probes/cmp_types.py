import json, sys
a = json.load(open(sys.argv[1])); b = json.load(open(sys.argv[2]))
diff_type = []; diff_val = []
for k in a:
    va, vb = a[k], b.get(k)
    ta = va.split(' ')[0]; tb = (vb or '').split(' ')[0]
    if ta != tb: diff_type.append((k, va, vb))
    elif va != vb: diff_val.append((k, va, vb))
print('total', len(a), 'type diffs', len(diff_type), 'value-repr diffs', len(diff_val))
for d in diff_type: print('T', d)
for d in diff_val: print('V', d)
