#!/bin/bash
# run_case.sh NAME THREADS ARGS...   (BIN=... selects another build, e.g. src/pf_mig_96)
#   Starts runs/NAME with src/pf_mig_128 and the given ARGS, or, if runs/NAME/checkpoint.bin
#   exists, resumes it bit-exactly (CSV rows written after the checkpoint are dropped first,
#   null bytes from a crash are stripped). Every command is appended to runs/NAME/run.sh.
NAME=$1; TH=$2; shift 2; ARGS="$*"
ROOT=$(cd "$(dirname "$0")/.." && pwd); D=$ROOT/runs/$NAME; mkdir -p $D; cd $D
BIN=${BIN:-$ROOT/src/pf_mig_128}
if [ -f checkpoint.bin ]; then
  tc=$(python3 -c "import struct;f=open('checkpoint.bin','rb');f.read(20);print(struct.unpack('d',f.read(8))[0])")
  python3 - "$tc" <<'PY'
import sys, os
tc = float(sys.argv[1])
for fn, col in (("bubble_metrics.csv", 1), ("sphere.csv", 1), ("gb_bubble.csv", 0), ("snapshots.csv", 0)):
    if not os.path.exists(fn): continue
    raw = open(fn, 'rb').read().replace(b'\x00', b'').decode()
    L = [l for l in raw.splitlines() if l.strip()]
    keep = [L[0]]
    for l in L[1:]:
        try:
            if float(l.split(',')[col]) <= tc + 1e-9: keep.append(l)
        except (ValueError, IndexError):
            pass
    open(fn, 'w').write('\n'.join(keep) + '\n')
PY
  # remove maps written after the checkpoint (they will be rewritten)
  for f in maps/t*.bin; do [ -e "$f" ] || continue; t=$(basename $f .bin); t=${t#t}; python3 -c "import sys;sys.exit(0 if float('$t')>$tc+1e-9 else 1)" && rm -f "$f"; done
  echo "OMP_NUM_THREADS=$TH $BIN restartIn=checkpoint.bin $ARGS   # resumed at t=$tc" >> run.sh
  OMP_NUM_THREADS=$TH exec $BIN restartIn=checkpoint.bin $ARGS >> run.log 2>&1
else
  echo "OMP_NUM_THREADS=$TH $BIN $ARGS" > run.sh
  OMP_NUM_THREADS=$TH exec $BIN $ARGS > run.log 2>&1
fi
