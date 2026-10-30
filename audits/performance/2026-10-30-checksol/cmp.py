import json, sys, collections
cases={c['id']:c for c in json.load(open(sys.argv[1]))}
a={r['id']:r for r in json.load(open(sys.argv[2]))}
b={r['id']:r for r in json.load(open(sys.argv[3]))}
def outc(r):
    if 'setup' in r: return ('setup',r['setup'])
    if 'err' in r: return ('err',r['err'],r.get('msg'))
    return ('r',r['r'])
cnt=collections.Counter(); bad=[]
for i in a:
    ra,rb=a[i],b.get(i)
    if rb is None: cnt['missing']+=1; continue
    oa,ob=outc(ra),outc(rb)
    if 'setup' in ra or 'setup' in rb:
        cnt['setup']+=1; 
        if oa!=ob: bad.append((i,'setup',oa,ob))
        continue
    if ra.get('err')=='Timeout' or rb.get('err')=='Timeout': cnt['timeout']+=1; bad.append((i,'timeout',oa,ob)); continue
    same = oa==ob and ra.get('warn')==rb.get('warn') and ra.get('after')==rb.get('after')
    if same: cnt['agree']+=1
    else:
        kind='result' if oa!=ob else ('warn' if ra.get('warn')!=rb.get('warn') else 'after')
        cnt['dis_'+kind]+=1; bad.append((i,kind,oa,ob,ra.get('warn'),rb.get('warn'),ra.get('after'),rb.get('after')))
print(dict(cnt))
seen=collections.Counter()
for x in bad:
    c=cases[x[0]]; key=str((c.get('f'),c.get('sym'),c.get('sol')))
    seen[key]+=1
    if seen[key]<=1: print(x[1], c.get('f'), c.get('sym'), c.get('sol'), c.get('flags'), x[2:4], (x[4:] if len(x)>4 and x[1]!='result' else ''))
print('distinct bases', len(seen))
