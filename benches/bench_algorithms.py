#!/usr/bin/env python3
"""Reproducible public-API latency comparison with independent result checks.

Run --help for screening, confirmation and immutable-package options. Native
implementation code is never edited by this runner. Each backend/trial has its
own process; only one measurement process runs at a time.
"""
from __future__ import annotations

import argparse
import gc
import gzip
import hashlib
import importlib
import json
import math
import os
from pathlib import Path
import platform
import resource
import statistics
import subprocess
import sys
import time

from performance_cases import Case, build, catalog

HERE=Path(__file__).resolve().parent


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(",",":"),allow_nan=False).encode()).hexdigest()


def write(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    pending=path.with_suffix(path.suffix+".pending")
    data=(json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+"\n").encode()
    pending.write_bytes(gzip.compress(data,mtime=0) if path.suffix==".gz" else data)
    pending.replace(path)


def read(path):
    path=Path(path);data=path.read_bytes()
    return json.loads(gzip.decompress(data) if path.suffix==".gz" else data)


def package_fingerprint(path):
    result=hashlib.sha256();count=0
    for file in sorted(Path(path).rglob("*")):
        if file.suffix not in (".py",".so"):continue
        result.update(str(file.relative_to(path)).encode()+b"\0")
        result.update(hashlib.sha256(file.read_bytes()).digest());count+=1
    return dict(sha256=result.hexdigest(),files=count)


def source_fingerprint():
    return {name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
            for name in ("bench_algorithms.py","performance_cases.py")}


def configuration(args):
    return {key:getattr(args,key) for key in
            ("samples","sample_seconds","warm_budget","timeout","memory_mib","cpu")}


def pin(cpu):
    if cpu is not None and hasattr(os,"sched_setaffinity"):
        os.sched_setaffinity(0,{cpu})


def probe(backend):
    module=importlib.import_module(backend)
    import importlib.metadata as metadata
    details=dict(library=backend,version=module.__version__,package=str(Path(module.__file__).resolve().parent),
                 python=sys.version,executable=sys.executable,
                 fingerprint=package_fingerprint(Path(module.__file__).resolve().parent))
    details["distribution_version"]=metadata.version(backend)
    details["dependencies"]={name:metadata.version(name) for name in ("mpmath","numpy")}
    if backend=="sympy":
        from sympy.external.gmpy import GROUND_TYPES
        from mpmath.libmp import BACKEND
        details.update(ground_types=GROUND_TYPES,mpmath_backend=BACKEND)
        if module.__version__!="1.14.0":raise ValueError("Reference must be SymPy 1.14.0")
    return details


class BlockReference:
    def find_spec(self,fullname,path=None,target=None):
        if fullname=="sympy" or fullname.startswith("sympy."):
            raise ModuleNotFoundError("Reference SymPy import forbidden during native measurement: "+fullname)


def evaluated_numeric(value):
    """Check native output before any exporter can evaluate held functions."""
    if isinstance(value,bool):return False
    if isinstance(value,(int,float,complex)):return True
    if type(value).__module__.startswith("mpmath") and (hasattr(value,"_mpf_") or hasattr(value,"_mpc_")):return True
    name=getattr(value,"class_name",type(value).__name__)
    if name in ("Float","Integer","Rational","Zero","One","NegativeOne","Half","ImaginaryUnit"):
        return True
    return name in ("Add","Mul") and all(evaluated_numeric(arg) for arg in value.args)


def worker(args):
    pin(args.cpu)
    if args.memory_mib:
        cap=args.memory_mib*1024*1024
        resource.setrlimit(resource.RLIMIT_AS,(cap,cap))
    spec=next(c for c in catalog() if c.id==args.worker_case)
    state=dict(status="running",backend=args.backend,case=spec.record(),
               samples=[],phase="import",pid=os.getpid(),gc_enabled=gc.isenabled())
    output=Path(args.worker_out)
    def checkpoint(phase):
        state.update(phase=phase,phase_started=time.perf_counter())
        write(output,state)
    guard=BlockReference() if args.backend=="symox" else None
    try:
        if guard:
            if any(name=="sympy" or name.startswith("sympy.") for name in sys.modules):
                raise RuntimeError("Reference preloaded before native measurement")
            sys.meta_path.insert(0,guard)
        checkpoint("import")
        module=importlib.import_module(args.backend)
        if args.backend=="sympy" and module.__version__!="1.14.0":
            raise ValueError("Reference must be SymPy 1.14.0")
        state.update(version=module.__version__,module_path=str(Path(module.__file__).resolve()))
        checkpoint("setup");start=time.perf_counter();function=build(spec,module)
        state["setup_seconds"]=time.perf_counter()-start
        state["setup_peak_rss_kib"]=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        checkpoint("first_call");start=time.perf_counter_ns();first=function()
        state["first_seconds"]=(time.perf_counter_ns()-start)*1e-9
        # Retain the first result for both backends. This deliberately measures
        # a repeated identical query with an already populated process cache.
        checkpoint("warm_calibration")
        loops=1;latest=first
        if state["first_seconds"]<args.warm_budget:
            while True:
                start=time.perf_counter_ns()
                for _ in range(loops):latest=function()
                elapsed=(time.perf_counter_ns()-start)*1e-9
                if elapsed>=args.sample_seconds/2 or loops>=32768:break
                loops*=2
        state["batch_loops"]=loops
        checkpoint("warm_measurement")
        warm_start=time.perf_counter()
        for _ in range(args.samples):
            start=time.perf_counter_ns()
            for _ in range(loops):latest=function()
            elapsed=(time.perf_counter_ns()-start)*1e-9
            state["samples"].append(dict(seconds_per_call=elapsed/loops,batch_seconds=elapsed,loops=loops))
            checkpoint("warm_measurement")
            if time.perf_counter()-warm_start>=args.warm_budget:break
        state["warm_seconds"]=statistics.median(sample["seconds_per_call"] for sample in state["samples"])
        state["measurement_peak_rss_kib"]=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        state["result_python_type"]=type(latest).__module__+"."+type(latest).__qualname__
        if spec.validation in ("numeric","numeric_roots"):
            values=[first,latest] if spec.validation=="numeric" else [*first,*latest]
            state["evaluated_numeric_output"]=all(evaluated_numeric(value) for value in values)
        if guard:sys.meta_path.remove(guard);guard=None
        checkpoint("serialization")
        from performance_validation import encode
        # Both the first and warmed answers are checked. A state-dependent
        # change in output cannot be hidden by only exporting the last answer.
        state["first_result"]=encode(first);state["result"]=encode(latest)
        state["first_result_sha256"]=digest(state["first_result"])
        state["result_sha256"]=digest(state["result"])
        state.update(status="ok",phase="complete")
    except MemoryError:
        state.update(status="memory_limit",error="Address-space budget exhausted")
    except BaseException as error:
        state.update(status="error",error=f"{type(error).__name__}: {str(error)[:800]}")
    finally:
        if guard:sys.meta_path.remove(guard)
        state["peak_rss_kib"]=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        write(output,state)


def environment(args):
    allowed=sorted(os.sched_getaffinity(0)) if hasattr(os,"sched_getaffinity") else []
    cpuinfo=Path("/proc/cpuinfo").read_text() if Path("/proc/cpuinfo").exists() else ""
    model=next((line.split(":",1)[1].strip() for line in cpuinfo.splitlines() if line.startswith("model name")),platform.processor())
    governor=Path(f"/sys/devices/system/cpu/cpu{args.cpu}/cpufreq/scaling_governor")
    return dict(platform=platform.platform(),cpu_model=model,allowed_cpus=allowed,pinned_cpu=args.cpu,
                cpu_governor=governor.read_text().strip() if governor.exists() else None,
                load_average=list(os.getloadavg()),clock=time.get_clock_info("perf_counter")._asdict()
                if hasattr(time.get_clock_info("perf_counter"),"_asdict") else
                dict(resolution=time.get_clock_info("perf_counter").resolution,monotonic=True),
                gc_policy="Keep Python GC enabled for both backends",
                cache_policy="Fresh process per trial; first call after input setup, then repeated identical calls; first result retained",
                thread_environment={name:"1" for name in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS")})


def child_env(args):
    env=dict(os.environ,PYTHONHASHSEED="0",OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1",
             MKL_NUM_THREADS="1",NUMEXPR_NUM_THREADS="1",MPLBACKEND="Agg")
    env.pop("SYMOX_DOWNSTREAM_BACKEND",None)
    if args.symox_package:env["PYTHONPATH"]=str(Path(args.symox_package).resolve().parent)
    return env


def run_trial(args,case,backend,index):
    path=args.work_dir/"workers"/f"{case.id}.{backend}.{index}.json"
    command=[sys.executable,str(Path(__file__).resolve()),"--worker-case",case.id,"--backend",backend,
             "--worker-out",str(path),"--samples",str(args.samples),"--sample-seconds",str(args.sample_seconds),
             "--warm-budget",str(args.warm_budget),"--memory-mib",str(args.memory_mib)]
    if args.cpu is not None:command.extend(["--cpu",str(args.cpu)])
    begin=time.perf_counter()
    try:
        process=subprocess.run(command,env=child_env(args),capture_output=True,text=True,timeout=args.timeout)
        state=read(path) if path.exists() else dict(status="worker_error",error=process.stderr[-1200:])
        if process.returncode:state.update(status="worker_error",exit_code=process.returncode,error=process.stderr[-1200:])
    except subprocess.TimeoutExpired:
        state=read(path) if path.exists() else {}
        phase=state.get("phase","startup")
        state.update(status="timeout",timeout_seconds=args.timeout)
        if phase=="first_call" and state.get("phase_started"):
            state["first_lower_bound_seconds"]=time.perf_counter()-state["phase_started"]
    state["wall_seconds"]=time.perf_counter()-begin
    state["backend"]=backend;state["trial_index"]=index
    return state


def verify_pair(args,case,a,b):
    if a.get("status")!="ok" or b.get("status")!="ok":
        return dict(status="not_comparable",detail="Both calls must finish and export their answers")
    if case.validation in ("numeric","numeric_roots"):
        if not a.get("evaluated_numeric_output") or not b.get("evaluated_numeric_output"):
            return dict(status="unevaluated_numeric",detail="A numerical API returned held work; exporting/evaluating it is not part of its measured cost")
    payload=args.work_dir/"validation-input.json"
    write(payload,dict(case=case.record(),native=a["result"],reference=b["result"],
                       native_first=a["first_result"],reference_first=b["first_result"]))
    command=[sys.executable,str(Path(__file__).resolve()),"--validate",str(payload)]
    try:
        result=subprocess.run(command,env=child_env(args),capture_output=True,text=True,timeout=args.validation_timeout)
    except subprocess.TimeoutExpired:
        return dict(status="validation_timeout",timeout_seconds=args.validation_timeout)
    if result.returncode:return dict(status="validation_error",error=result.stderr[-1200:])
    try:return json.loads(result.stdout)
    except ValueError:return dict(status="validation_error",error=result.stdout[-1200:])


def validation_worker(path):
    from performance_validation import verify
    data=read(path);case=Case(**data["case"])
    checks=[verify(case,data["native"],data["reference"])]
    if (data["native_first"]!=data["native"] or data["reference_first"]!=data["reference"]):
        checks.append(verify(case,data["native_first"],data["reference_first"]))
    bad=next((check for check in checks if check["status"]!="verified"),None)
    return bad or dict(status="verified",checks=checks)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path)
    parser.add_argument("--work-dir",type=Path)
    parser.add_argument("--symox-package",type=Path,help="Immutable package directory containing __init__.py")
    parser.add_argument("--provenance",type=Path,help="Optional validated-wheel baseline record")
    parser.add_argument("--case",action="append",help="Case ID or family prefix; repeat to select several")
    parser.add_argument("--list",action="store_true")
    parser.add_argument("--rounds",type=int,default=3)
    parser.add_argument("--samples",type=int,default=3)
    parser.add_argument("--sample-seconds",type=float,default=0.01)
    parser.add_argument("--warm-budget",type=float,default=0.15)
    parser.add_argument("--timeout",type=float,default=20)
    parser.add_argument("--validation-timeout",type=float,default=20)
    parser.add_argument("--memory-mib",type=int,default=4096)
    parser.add_argument("--cpu",type=int)
    parser.add_argument("--resume",action="store_true")
    parser.add_argument("--refine-below",type=float,help="On resume, add rounds only when first-call or warm speedup is below this threshold")
    parser.add_argument("--worker-case",help=argparse.SUPPRESS)
    parser.add_argument("--worker-out",type=Path,help=argparse.SUPPRESS)
    parser.add_argument("--backend",choices=("symox","sympy"),help=argparse.SUPPRESS)
    parser.add_argument("--probe",choices=("symox","sympy"),help=argparse.SUPPRESS)
    parser.add_argument("--validate",type=Path,help=argparse.SUPPRESS)
    args=parser.parse_args()
    if args.probe:print(json.dumps(probe(args.probe)));return
    if args.validate:print(json.dumps(validation_worker(args.validate)));return
    if args.worker_case:worker(args);return
    if min(args.rounds,args.samples)<1 or min(args.sample_seconds,args.warm_budget,args.timeout,args.validation_timeout)<=0:
        parser.error("Counts and time budgets must be positive")
    selected=[c for c in catalog() if not args.case or any(c.id==key or c.id.startswith(key+".") for key in args.case)]
    if not selected:parser.error("No matching cases")
    if args.list:
        for case in selected:print(json.dumps(case.record()))
        return
    if args.out is None or args.work_dir is None:parser.error("--out and --work-dir are required")
    args.out=args.out.resolve();args.work_dir=args.work_dir.resolve();args.work_dir.mkdir(parents=True,exist_ok=True)
    if args.cpu is None and hasattr(os,"sched_getaffinity"):args.cpu=min(os.sched_getaffinity(0))
    backend_info={}
    for name in ("symox","sympy"):
        process=subprocess.run([sys.executable,str(Path(__file__).resolve()),"--probe",name],
                               env=child_env(args),capture_output=True,text=True,check=True)
        backend_info[name]=json.loads(process.stdout)
    baseline=None
    if args.provenance:
        baseline_record=read(args.provenance);candidate=baseline_record["candidate"]
        for name,expected in candidate["wheel_package_files_sha256"].items():
            file=Path(backend_info["symox"]["package"]).parent/name
            if hashlib.sha256(file.read_bytes()).hexdigest()!=expected:raise ValueError("Baseline package differs from its validated wheel")
        baseline={key:candidate[key] for key in ("wheel_sha256","production_fix_commit","build_source_snapshot_commit","actual_symox_version")}
    identity=dict(timing_sources=source_fingerprint(),backends=backend_info,configuration=configuration(args),cases=[c.record() for c in selected],baseline=baseline)
    if args.resume and args.out.exists():
        report=read(args.out)
        if report["identity"]!=identity:raise ValueError("Resume refused: package, cases, environment or timing protocol changed")
    else:
        report=dict(schema=1,identity=identity,environment=environment(args),results={},outputs={},
                    validation_source_sha256=hashlib.sha256((HERE/"performance_validation.py").read_bytes()).hexdigest(),
                    limitations=["Selected deterministic workloads, not universal algorithmic parity or an asymptotic-complexity proof.",
                                 "First call is after input construction; not a cold operating-system page cache.",
                                 "Repeated-call results measure populated process caches; no asymmetric cache clearing.",
                                 "Peak RSS is whole-process high-water memory, including imports and result export.",
                                 "Concurrent host activity can affect timings; order alternates and repeated processes expose variation."])
    report.setdefault("invocations",[]).append(dict(argv=sys.argv,started_unix=time.time(),rounds=args.rounds))
    for case in selected:
        row=report["results"].setdefault(case.id,dict(case=case.record(),trials=[],validations={}))
        if row["trials"]:
            if any(t["backends"][name].get("status")!="ok" for t in row["trials"] for name in ("symox","sympy")):continue
            if args.refine_below is not None:
                ratios=[statistics.median(t["backends"]["sympy"][mode]/t["backends"]["symox"][mode]
                                          for t in row["trials"]) for mode in ("first_seconds","warm_seconds")]
                if min(ratios)>args.refine_below:continue
        for index in range(len(row["trials"]),args.rounds):
            order=("symox","sympy") if index%2==0 else ("sympy","symox")
            trial={name:run_trial(args,case,name,index) for name in order}
            a,b=trial["symox"],trial["sympy"]
            # Exported values are interned once in the report, not duplicated
            # for every trial or replaced with hashes that cannot be checked.
            pair_key=digest([a.get("first_result_sha256"),a.get("result_sha256"),b.get("first_result_sha256"),b.get("result_sha256")])
            if pair_key not in row["validations"]:
                row["validations"][pair_key]=verify_pair(args,case,a,b)
            for state in trial.values():
                for field in ("result","first_result"):
                    value=state.pop(field,None)
                    if value is not None:report["outputs"][digest(value)]=value
            row["trials"].append(dict(backends=trial,validation_key=pair_key,backend_order=list(order)))
            write(args.out,report)
            print(case.id,f"round {index+1}/{args.rounds}",row["validations"][pair_key]["status"],flush=True)
            # Repeated timeouts/errors are not measurements of fast algorithms.
            # One failed attempt is enough for a screening gap; rerun explicitly
            # with a larger budget if the case is important.
            if a.get("status")!="ok" or b.get("status")!="ok":break
    final={}
    for name in ("symox","sympy"):
        process=subprocess.run([sys.executable,str(Path(__file__).resolve()),"--probe",name],
                               env=child_env(args),capture_output=True,text=True,check=True)
        final[name]=json.loads(process.stdout)
    report["environment_unchanged"]=final==backend_info
    report["environment"]["load_average_end"]=list(os.getloadavg())
    write(args.out,report)
    if not report["environment_unchanged"]:raise RuntimeError("Backend environment changed during measurement")


if __name__=="__main__":main()
