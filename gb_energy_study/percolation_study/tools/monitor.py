"""Background monitor: (1) converts each c snapshot (float32 [z][y][x]) into a column-thickness
map (# cells with c>=0.5, uint8 [y][x]) saved in thick/, then deletes the 4 MB raw file;
(2) stops a run once number_of_projected_bubbles == 1 has held for HOLD time units
(the run is resumable from checkpoint.bin, so nothing is lost)."""
import os, glob, time, signal, numpy as np, csv
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NX, NY, NZ = 128, 128, 64; HOLD = 3.0
TCAP = float(os.environ.get('TCAP', '1e9'))   # hard time cap (stop + resumable)
while True:
    active = 0
    for d in sorted(glob.glob(f"{ROOT}/runs/N*_bGB*")):
        os.makedirs(f"{d}/thick", exist_ok=True)
        snaps = {}
        if os.path.exists(f"{d}/snapshots.csv"):
            for r in csv.DictReader(open(f"{d}/snapshots.csv")): snaps[r["c_file"]] = float(r["time_dimensionless"])
        for f in glob.glob(f"{d}/morphology_step_*.raw"):
            name = os.path.basename(f)
            if name not in snaps or os.path.getsize(f) != 4*NX*NY*NZ or time.time() - os.path.getmtime(f) < 3: continue
            c = np.fromfile(f, np.float32).reshape(NZ, NY, NX)
            np.save(f"{d}/thick/t{snaps[name]:07.2f}.npy", (c >= 0.5).sum(0).astype(np.uint8))
            os.remove(f)
        pidf = f"{d}/pid"
        if not os.path.exists(pidf): continue
        pid = int(open(pidf).read())
        alive = os.path.exists(f"/proc/{pid}")
        if not alive: continue
        active += 1
        rows = list(csv.reader(open(f"{d}/bubble_metrics.csv")))[1:]
        t1 = None
        for r in rows:
            n, t = int(r[4]), float(r[1])
            if n == 1 and t1 is None: t1 = t
            if n != 1: t1 = None
        if t1 is not None and float(rows[-1][1]) >= t1 + HOLD:
            os.kill(pid, signal.SIGTERM)
            open(f"{d}/STOPPED_single_bubble.txt", "w").write(f"count==1 from t={t1} held to t={rows[-1][1]}; stopped by monitor\n")
        elif rows and float(rows[-1][1]) >= TCAP:
            os.kill(pid, signal.SIGTERM)
            open(f"{d}/STOPPED_time_cap.txt", "w").write(f"reached time cap t={rows[-1][1]} with {rows[-1][4]} bubbles; stopped by monitor (resumable)\n")
    open(f"{ROOT}/monitor_heartbeat", "w").write(time.ctime() + f" active={active}\n")
    time.sleep(30)
