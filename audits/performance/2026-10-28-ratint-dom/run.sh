#!/bin/bash
# usage: run.sh <python dir before> <python dir after> > raw.jsonl
H=$(dirname "$0")/harness.py
for i in $(seq 0 10); do
  for side in before after; do
    if [ $side = before ]; then P=$1; else P=$2; fi
    (cd /tmp && PYTHONPATH=$P timeout 600 /home/keckj/Documents/symox/.venv/bin/python $H $side $i)
  done
done
