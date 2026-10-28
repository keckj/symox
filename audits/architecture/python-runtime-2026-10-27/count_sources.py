#!/usr/bin/env python3
"""Reproducible lexical production SLOC; never a mathematical-parity metric."""
import ast,collections,io,json,re,sys,tokenize
from pathlib import Path

root=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).parent/'input'
output=Path(__file__).parent

def python_count(text):
    tree=ast.parse(text)
    docs=set()
    for node in ast.walk(tree):
        if isinstance(node,(ast.Module,ast.ClassDef,ast.FunctionDef,ast.AsyncFunctionDef)) and node.body:
            first=node.body[0]
            if isinstance(first,ast.Expr) and isinstance(first.value,ast.Constant) and isinstance(first.value.value,str):
                docs.update(range(first.lineno,first.end_lineno+1))
    code=set()
    ignored={tokenize.ENCODING,tokenize.ENDMARKER,tokenize.NL,tokenize.NEWLINE,tokenize.INDENT,tokenize.DEDENT,tokenize.COMMENT}
    for token in tokenize.generate_tokens(io.StringIO(text).readline):
        if token.type not in ignored:
            rows=(token.start[0],) if token.type==tokenize.STRING else range(token.start[0],token.end[0]+1)
            code.update(i for i in rows if i not in docs)
    loops=sum(isinstance(n,(ast.For,ast.AsyncFor,ast.While,ast.comprehension)) for n in ast.walk(tree))
    functions=[]
    for n in ast.walk(tree):
        if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef)):
            calls=[]
            for call in ast.walk(n):
                if isinstance(call,ast.Call):
                    name=ast.unparse(call.func)
                    if '_core.' in name or name.startswith('_native_'):
                        calls.append({'line':call.lineno,'call':name})
            functions.append({'name':n.name,'line':n.lineno,'end':n.end_lineno,'loops':sum(isinstance(v,(ast.For,ast.While,ast.comprehension)) for v in ast.walk(n)),'native_calls':calls})
    return len(code),len(docs),loops,functions

def rust_without_comments_and_literals(text):
    # Preserve line numbers; literals remain one significant marker on their
    # first line, with contents masked so braces cannot affect test parsing.
    chars=list(text);i=0;n=len(text)
    def mask(a,b,marker=False):
        for j in range(a,b):
            if chars[j]!='\n':chars[j]=' '
        if marker:chars[a]='S'
    while i<n:
        if text.startswith('//',i):
            end=text.find('\n',i);end=n if end<0 else end;mask(i,end);i=end;continue
        if text.startswith('/*',i):
            end=i+2;depth=1
            while end<n and depth:
                if text.startswith('/*',end):depth+=1;end+=2
                elif text.startswith('*/',end):depth-=1;end+=2
                else:end+=1
            mask(i,end);i=end;continue
        raw=re.match(r'(?:br|r)(#{0,255})"',text[i:])
        if raw:
            close='"'+raw.group(1);end=text.find(close,i+len(raw.group(0)));end=n if end<0 else end+len(close)
            mask(i,end,True);i=end;continue
        if text[i]=='"' or text.startswith('b"',i):
            end=i+(2 if text.startswith('b"',i) else 1)
            while end<n:
                if text[end]=='\\':end+=2
                elif text[end]=='"':end+=1;break
                else:end+=1
            mask(i,min(end,n),True);i=end;continue
        char=re.match(r"b?'(?:\\(?:u\{[0-9a-fA-F_]+\}|x[0-9a-fA-F]{2}|.)|[^'\\\n])'",text[i:])
        if char:
            end=i+len(char.group());mask(i,end,True);i=end;continue
        i+=1
    return ''.join(chars)

def rust_count(text):
    clean=rust_without_comments_and_literals(text);chars=list(clean);removed=set()
    pattern=re.compile(r'#\s*\[\s*(?:cfg\s*\(\s*test\s*\)|test)\s*\]')
    for match in pattern.finditer(clean):
        if match.start() in removed:continue
        start=match.start();opening=clean.find('{',match.end());semi=clean.find(';',match.end())
        if semi>=0 and (opening<0 or semi<opening):end=semi+1
        elif opening>=0:
            level=1;end=opening+1
            while end<len(clean) and level:
                level+=(clean[end]=='{')-(clean[end]=='}');end+=1
        else:continue
        for i in range(start,end):
            removed.add(i)
            if chars[i]!='\n':chars[i]=' '
    before=sum(bool(line.strip()) for line in clean.splitlines())
    after=sum(bool(line.strip()) for line in ''.join(chars).splitlines())
    return after,before-after

records=[]
for prefix,suffix in [('python/symox','.py'),('crates/symox-core/src','.rs'),('crates/symox-py/src','.rs')]:
    for path in sorted((root/prefix).rglob('*'+suffix)):
        relative=path.relative_to(root);parts=relative.parts
        if any(p in ('tests','test','testing','__pycache__','target','docs','benches','benchmarks','examples') for p in parts):continue
        if suffix=='.rs' and (path.name=='tests.rs' or path.stem.endswith(('_tests','_test'))):continue
        text=path.read_text();record={'path':str(relative),'physical_lines':len(text.splitlines()),'language':'Python' if suffix=='.py' else 'Rust'}
        if suffix=='.py':
            sloc,docs,loops,functions=python_count(text);record.update(sloc=sloc,docstring_lines=docs,loop_nodes=loops,functions=functions)
        else:
            sloc,test=rust_count(text);record.update(sloc=sloc,inline_test_sloc_removed=test)
        records.append(record)
aggregate={}
for lang in ('Python','Rust'):
    selected=[r for r in records if r['language']==lang];aggregate[lang]={'files':len(selected),'physical_lines':sum(r['physical_lines'] for r in selected),'sloc':sum(r['sloc'] for r in selected)}
def region(path):
    for prefix in ('python/symox/solvers/ode/','python/symox/solvers/','python/symox/integrals/','python/symox/polys/matrices/','python/symox/polys/_dense/','python/symox/polys/','python/symox/stats/','python/symox/matrices/'):
        if path.startswith(prefix):return prefix
    return 'other Python'
regions=collections.defaultdict(lambda:{'files':0,'sloc':0,'physical_lines':0})
for record in records:
    if record['language']=='Python':
        bucket=regions[region(record['path'])];bucket['files']+=1;bucket['sloc']+=record['sloc'];bucket['physical_lines']+=record['physical_lines']
(output/'counts.json').write_text(json.dumps({'root':str(root),'aggregate':aggregate,'python_regions':dict(regions),'files':records},indent=2)+'\n')
(output/'counts.tsv').write_text('language\tsloc\tphysical_lines\tpath\n'+''.join(f"{r['language']}\t{r['sloc']}\t{r['physical_lines']}\t{r['path']}\n" for r in records))
print(json.dumps({'aggregate':aggregate,'python_regions':dict(regions)},indent=2))
