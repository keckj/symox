import json,re,sys
b=json.load(open(sys.argv[1])); n=json.load(open(sys.argv[2]))
def norm(s):
    m={}
    return re.sub(r"dummy_index=(\d+)",lambda g:"d%d"%m.setdefault(g.group(1),len(m)),s)
same=diff=0
for k in b:
    if norm(b[k][0])==norm(n[k][0]):
        same+=1
        r=b[k][1]/max(n[k][1],1e-6)
        if r>1.3 or r<0.77: print(f"{k:45s} base {b[k][1]:.4f} new {n[k][1]:.4f} x{r:.2f}")
    else:
        diff+=1; print("DIFF",k,"\n  base",b[k][0][:400],"\n  new ",n[k][0][:400])
print("same",same,"diff",diff)
