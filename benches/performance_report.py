#!/usr/bin/env python3
"""Summarize comparable timings or recheck retained answers without retiming."""
import argparse
import copy
import hashlib
from pathlib import Path
import random
import statistics
import time

from bench_algorithms import read,write,verify_pair
from performance_cases import Case


def summaries(report):
    rows=[]
    for name,row in sorted(report["results"].items()):
        statuses=sorted({row["validations"][trial["validation_key"]]["status"] for trial in row["trials"]})
        result=dict(case=name,statuses=statuses,trials=len(row["trials"]),comparable=False)
        if statuses==["verified"] and report.get("environment_unchanged") is True:
            trials=[trial["backends"] for trial in row["trials"]]
            result["comparable"]=True
            for mode in ("first","warm"):
                ratios=[t["sympy"][mode+"_seconds"]/t["symox"][mode+"_seconds"] for t in trials]
                rng=random.Random(int.from_bytes(hashlib.sha256((name+mode).encode()).digest()[:8],"big"))
                boot=sorted(statistics.median(rng.choices(ratios,k=len(ratios))) for _ in range(1000))
                result[mode]=dict(speedup=statistics.median(ratios),
                    bootstrap_95_percent=[boot[24],boot[974]] if len(ratios)>=3 else None,
                    native_seconds=statistics.median(t["symox"][mode+"_seconds"] for t in trials),
                    reference_seconds=statistics.median(t["sympy"][mode+"_seconds"] for t in trials))
        else:
            result["backend_states"]=[{b:{key:s[key] for key in ("status","phase","error","timeout_seconds") if key in s}
                                        for b,s in trial["backends"].items()} for trial in row["trials"]]
        rows.append(result)
    return rows


def markdown(report,rows,threshold):
    info=report["identity"]["backends"]
    comparable=[r for r in rows if r["comparable"]]
    weak=sorted((r for r in comparable if min(r[m]["speedup"] for m in ("first","warm"))<=threshold),
                key=lambda r:min(r[m]["speedup"] for m in ("first","warm")))
    lines=["# Public algorithm performance screening","",
        f"Reference: SymPy {info['sympy']['version']}. Symox distribution: {info['symox']['distribution_version']}.",
        f"Package fingerprint: `{info['symox']['fingerprint']['sha256']}`.","",
        f"{len(rows)} workloads; {len(comparable)} verified and comparable; {len(rows)-len(comparable)} excluded.",
        f"{len(weak)} comparable workloads have a first-call or repeated-call speedup ≤ {threshold:g}.","",
        "Speedup means reference time / native time; values below 1 mean Symox is slower.",
        "Fresh processes alternate backend order on one CPU. Inputs and result verification are outside timing.",
        "Repeated calls retain both libraries' caches. These measurements are not a universal parity or complexity certificate.",
        "Three-process screenings are provisional; bootstrap intervals describe observed variation, not host isolation.","",
        "## Weak or marginal cases","",
        "| Workload | First speedup | Repeated speedup | Native first, ms | Native repeated, ms | Processes |",
        "|---|---:|---:|---:|---:|---:|"]
    for r in weak:
        lines.append(f"| `{r['case']}` | {r['first']['speedup']:.3f} | {r['warm']['speedup']:.3f} | {r['first']['native_seconds']*1000:.6g} | {r['warm']['native_seconds']*1000:.6g} | {r['trials']} |")
    lines.extend(["","## Excluded results","","These statuses require diagnosis; a validation timeout or alternate formatting is not evidence of an incorrect native algorithm.","",
                  "| Workload | Status |","|---|---|"])
    for r in rows:
        if not r["comparable"]:lines.append(f"| `{r['case']}` | {', '.join(r['statuses'])} |")
    lines.extend(["","## All comparable cases","","| Workload | First speedup | Repeated speedup |","|---|---:|---:|"])
    for r in comparable:lines.append(f"| `{r['case']}` | {r['first']['speedup']:.3f} | {r['warm']['speedup']:.3f} |")
    return "\n".join(lines)+"\n"


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input",type=Path)
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--summary-json",type=Path)
    parser.add_argument("--threshold",type=float,default=1.5)
    parser.add_argument("--revalidate",action="store_true")
    parser.add_argument("--case",action="append")
    parser.add_argument("--validation-timeout",type=float,default=60)
    args=parser.parse_args();report=read(args.input)
    if args.revalidate:
        if args.out.resolve()==args.input.resolve():parser.error("Keep the original evidence; choose a new output")
        args.work_dir=args.out.parent/(args.out.name+".work");args.work_dir.mkdir(parents=True,exist_ok=True)
        args.symox_package=None
        source_sha=hashlib.sha256((Path(__file__).parent/"performance_validation.py").read_bytes()).hexdigest()
        for name,row in report["results"].items():
            if args.case and not any(name==key or name.startswith(key+".") for key in args.case):continue
            visited=set()
            for trial in row["trials"]:
                key=trial["validation_key"]
                if key in visited:continue
                visited.add(key)
                states=copy.deepcopy(trial["backends"])
                for state in states.values():
                    for field in ("result","first_result"):
                        if field+"_sha256" in state:state[field]=report["outputs"][state[field+"_sha256"]]
                result=verify_pair(args,Case(**row["case"]),states["symox"],states["sympy"])
                report.setdefault("validation_updates",[]).append(dict(case=name,key=key,
                    previous=row["validations"][key],replacement=result,source_sha256=source_sha,unix_time=time.time()))
                row["validations"][key]=result
                print(name,result["status"],flush=True)
        write(args.out,report)
    else:
        rows=summaries(report);args.out.parent.mkdir(parents=True,exist_ok=True)
        args.out.write_text(markdown(report,rows,args.threshold))
        if args.summary_json:write(args.summary_json,dict(cases=rows,threshold=args.threshold))


if __name__=="__main__":main()
