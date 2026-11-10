# Percolation begin and percolation end vs contact angle (N = 33 and N = 25)

**Question.** For a fixed bubble arrangement, is the projected coverage at which percolation *begins* (first merger) and at which it *ends* (only one bubble left) the same for every contact angle, with the angle changing only the time? Is this true for 25 bubbles as well as for 33?

**Short answer (all values measured from the simulations below).**

* **Percolation begin.** N = 33: confirmed, coverage 0.390–0.394 for all four angles (1.0 % spread). N = 25: **not** the same. Coverage 0.320 → 0.344, rising steadily with θ (7.2 % spread), much more than the output resolution.
* **Percolation end.** **Not** angle-independent. Only θ = 104° and 120° ever reach a single bubble, at coverage 0.704 / 0.716 (N = 33) and 0.640 / 0.650 (N = 25). The two angles agree within about 1.5 %, but the single-bubble time differs by a factor of 1.9.
  * θ = 128° and 148° **never reach a single bubble** by t = 150. That is 5–6× longer than the 104° case needed. 4–8 islands remain, with smallest gaps of 3.5–11 cells at t = 150.
  * Their coverage peaks at 0.50–0.63, below the single-bubble coverage of the low-angle cases, and then *falls*.
  * For high angles the end state is a set of separate, rounded lenses, not one connected bubble.

---

## 1. What was run

Everything is identical to the original study (`../src/phasefield_gb_bubbles3d.cpp`, unmodified source, recompiled with the same flags), except that each run continues until one bubble remains.

| Fixed | Value |
|---|---|
| Model | conserved Cahn–Hilliard (gas) + Allen–Cahn (grains), no source, no clamp |
| Grid | 128 × 128 × 64, periodic x, y, no-flux z |
| A0, bS, j, jGB | 4, 1, 2, 2 |
| Ms, Mgb, Mb, LGB | 0.05, 0.5, 5, 1 |
| cInit, r0, seed, minSep | 0.12, 4.5, 11, 19 |
| Output | metrics every 0.1, thickness maps every 1.0 |
| Varied | bGB = 1.2, 1.5, 2.0, 3.0 → measured θ = 148°, 128°, 120°, 104° (±6°, from the original 2-D lens study) |

**Initial conditions.** `ic/N33/ic.bin` is byte-for-byte the same centre list as the original `common_ic.bin`. `ic/N25/ic.bin` was made by the unmodified generator with `nBubbles=25`, seed 11. Its centres are exactly the first 25 of the 33 centres (checked with `diff`). The same subset approach was used for the earlier N = 20 study.

**Stopping rule.** A background monitor stopped a run once the projected bubble count had been 1 for 3.0 time units, or at a hard cap of t = 150. All runs can be resumed from `checkpoint.bin`.

**Definitions** (unchanged from the original study):
* A column is covered when ≥ 3 cells have c ≥ 0.5.
* Components are 8-connected and periodic. Coverage = covered columns / 128².
* **Percolation begin** = first 0.1-spaced output with fewer components than at t = 0.
* **Percolation end** = first output with exactly one component.

## 2. Verification (done before any result was used)

* The N = 33 re-runs reproduce the archived original runs exactly: step, bubble count, total and largest area, coverage, gas and free energy are identical in every row with a matching time stamp. That is 253–380 rows per case, up to t = 40, with 0 mismatches.
* Gas conservation: max |mass drift| ≤ 1.3e-13 in all 8 runs. Free energy never increases (0 increases in 274–1502 outputs). No non-finite values. Bubble count never above the initial value. No duplicated or out-of-order rows.
* **Incident.** On 9 Nov the computer shut down while four runs were writing; each CSV ended in about 240 null bytes.
  * Repair: the null bytes were stripped and the runs resumed from the last checkpoint. Original files are kept in `runs/*/crash_backup/`.
  * Check: the restarted runs reproduced, bit for bit, the rows written after the checkpoint that had survived the crash. The restart is therefore exact.
* "Single bubble" is a final state, not a transient: in all 4 cases that reach it, the count stays at 1 for the remaining 3 time units.

## 3. Results

### 3.1 Percolation begin (first merger)

| θ | N = 33 coverage | N = 33 time | N = 25 coverage | N = 25 time |
|---|---|---|---|---|
| 148° | 0.3941 | 27.6 | 0.3435 | 26.7 |
| 128° | 0.3940 | 16.9 | 0.3344 | 16.6 |
| 120° | 0.3921 | 11.3 | 0.3259 | 11.2 |
| 104° | 0.3903 | 8.0 | 0.3198 | 7.9 |
| range / mean | **0.0038 / 0.393 = 1.0 %** | 3.45× | **0.0237 / 0.331 = 7.2 %** | 3.38× |

Resolution: the coverage at the previous output is at most 0.0027 lower than the event value. The N = 25 spread is 9× larger than that and strictly monotone in θ.

**Conclusion.** The near-constant onset coverage seen for 33 bubbles does **not** hold for 25 bubbles. The 25-bubble result (7.2 %, rising with θ) is like the earlier 20-bubble result (6.1 %, rising with θ). N = 33 is the exception, not the rule. The time to onset changes about 3.4× in both populations.

### 3.2 Percolation end (single bubble left)

| θ | N = 33: single bubble? | coverage | time | N = 25: single bubble? | coverage | time |
|---|---|---|---|---|---|---|
| 104° | yes | **0.7042** | 24.1 | yes | **0.6403** | 29.6 |
| 120° | yes | **0.7158** | 45.6 | yes | **0.6498** | 55.7 |
| 128° | **no**, 4 islands at t = 150 | (max ever 0.635, t = 57) | – | **no**, 4 islands at t = 150 | (max ever 0.581, t = 73) | – |
| 148° | **no**, 8 islands at t = 150 | (max ever 0.533, t = 63) | – | **no**, 6 islands at t = 150 | (max ever 0.500, t = 81) | – |

* **Between 104° and 120°** the end coverage differs by 0.012 (1.6 %, N = 33) and 0.010 (1.5 %, N = 25), while the time nearly doubles (×1.89 and ×1.88). For these two angles alone, the statement "same coverage, different time" holds approximately. The end coverage depends on N (0.70–0.72 vs 0.64–0.65), as the onset coverage does.
* **For 128° and 148°** the statement fails. The highest coverage these runs ever reach is below the single-bubble coverage of the low-angle runs with the same N. After that peak, coverage *decreases*:
  * N = 33: 0.635 → 0.542 at 128° and 0.533 → 0.455 at 148°.
  * N = 25: 0.581 → 0.529 at 128° and 0.500 → 0.430 at 148°.
  * Gas uptake meanwhile continues slowly (porosity still rising, 0.097–0.109 at t = 150).
  * The lenses become more compact: the same gas occupies less GB area. They cannot become one bubble "at the same coverage, just later", because the coverage moves away from that value.

### 3.3 Why the high-angle runs stop merging

From the thickness maps (every 1.0 time unit; `results/gaps.csv`, `results/events_*.json`):

* **Gap between islands.**
  * At 104° and 120° the smallest gap between the last two islands shrinks to 1–2 cells while coverage is still growing, and they merge.
  * At 128° the smallest gap grows after t ≈ 50: from 4.7 to 5.7 cells (N = 33) and from 4 to 11 cells (N = 25).
  * At 148° it stays mostly at 3–5 cells. Briefly (t = 148–149, N = 25) it closes to 1.8 cells; that pair then merges at t = 150.
* **How islands are lost after t ≈ 60.**
  * 128°: only by occasional mergers (N = 25: 2 mergers after t = 60; N = 33: none, 4 islands from t = 50 to t = 150).
  * 148°: mostly by **dissolution** of small islands (Ostwald ripening). N = 33: 5 islands disappear without merging between t = 89 and 119, at 18 → 6 columns. N = 25: 3 disappear between t = 103 and 145.
  * No dissolution happens at 104° or 120°; there every loss of an island is a merger.
* **Interpretation (hypothesis, not separately tested).** At high θ, a merged pair relaxes toward a rounder, more compact lens, so it pulls back from its neighbours. Once the supersaturation is used up, nothing drives the footprints together again. The remaining route to a single bubble is ripening, which is much slower and works by islands disappearing, not by islands joining. At low θ the lenses stay flat and spread out, so the network closes while the gas is still being absorbed.

### 3.4 Other milestones (for completeness)

| θ | N = 33: A1/Atot > 0.5 (cov, t) | N = 25: A1/Atot > 0.5 (cov, t) | N = 33: spanning cluster (cov, t) | N = 25: spanning cluster (cov, t) |
|---|---|---|---|---|
| 148° | not reached | not reached | not reached | not reached |
| 128° | 0.528, 28.2 | 0.558, 53.4 | not reached | 0.562, 55 |
| 120° | 0.496, 16.7 | 0.513, 25.6 | 0.632, 28 | 0.543, 29 |
| 104° | 0.479, 10.9 | 0.503, 17.1 | 0.608, 17 | 0.518, 18 |

"Spanning" means one component wraps the periodic box (checked on 1.0-spaced maps). These intermediate milestones also depend on θ: higher θ needs more coverage.

### 3.5 Sensitivity to the detection rule (1.0-spaced maps, so times are rounded up to the next integer)

Coverage at the event is measured with the original (≥ 3 cells) definition, at the map where each rule first fires:

| N | θ | begin: ≥1 cell | begin: ≥3 | begin: ≥5 | end: ≥1 cell | end: ≥3 |
|---|---|---|---|---|---|---|
| 33 | 148° | 0.390 | 0.397 | 0.403 | – | – |
| 33 | 128° | 0.384 | 0.396 | 0.419 | – | – |
| 33 | 120° | 0.370 | 0.404 | 0.444 | 0.707 | 0.717 |
| 33 | 104° | 0.365 | 0.390 | 0.481 | 0.680 | 0.715 |
| 25 | 148° | 0.340 | 0.346 | 0.360 | – | – |
| 25 | 128° | 0.327 | 0.340 | 0.363 | – | – |
| 25 | 120° | 0.307 | 0.341 | 0.377 | 0.645 | 0.651 |
| 25 | 104° | 0.299 | 0.322 | 0.398 | 0.619 | 0.643 |

* The ≥ 5-cell rule never gives a single component within the saved window. Runs were stopped 3 time units after the ≥ 3 rule fired, so this is "not resolved", not a result.
* With thin (≥ 1 cell) contact as the criterion, onset coverage rises with θ in both populations. This agrees with the earlier fine-resolution study for N = 33 (13 % at first touch). The single-bubble coverage of 104° vs 120° differs by 4 % (vs 1.6 % with the 3-cell rule). The flat N = 33 onset under the 3-cell rule is the same compensation effect found before.

## 4. Conclusions

1. **Percolation begin, N = 33**: approximately independent of θ (1.0 %), as found before.
2. **Percolation begin, N = 25**: depends on θ, 7.2 %, monotone, like N = 20 (6.1 %). The flat 33-bubble result does not carry over to other bubble counts in this layout.
3. **Percolation end**: approximately the same coverage for θ = 104° and 120° (1.5–1.6 %, with the time ×1.9). It is **not reached at all** for θ = 128° and 148° within t = 150, because the islands stop approaching each other and their footprints shrink. The end state itself depends on the contact angle, so "the end coverage is the same, only the time differs" is **not** supported.
4. In every case the angle changes the *time* strongly (×3.4 at onset; ×1.9 between the two angles that reach the end).

## 5. Limitations

* One bubble arrangement (seed 11) per N; N = 25 is a subset of N = 33, so the two are not independent samples.
* θ carries ±6°. There are four deterministic runs per N, so no statistical claim can be made.
* t = 150 is a finite horizon. The 128°/148° systems could still reach one bubble much later through ripening. If they do, it would be through islands dissolving, at a coverage below the low-angle end values, unless the current downward trend in coverage reverses.
* The rule-sensitivity table uses 1.0-spaced maps. In the fastest case (104°) coverage changes by up to ~0.01 per time unit near the events, so those numbers carry about that uncertainty. The main tables use 0.1-spaced outputs.

## 6. Files

| Path | Content |
|---|---|
| `runs/N{33,25}_bGB{1.2,1.5,2.0,3.0}/` | `bubble_metrics.csv` (every 0.1), `thick/t*.npy` (column-thickness maps every 1.0), `run.sh` (exact commands incl. restarts), `run.log`, `params*.txt`, `checkpoint.bin` (resumable), `STOPPED_*.txt` (why the run stopped) |
| `results/events_summary.csv` | all event times/coverages per run |
| `results/rules.csv`, `results/gaps.csv`, `results/events_*.json` | rule sensitivity, smallest gap vs time, merger/dissolution log |
| `results/fig1_count_and_coverage.png` | bubble count vs coverage; coverage vs time with begin/end markers |
| `results/fig2_events_vs_angle.png` | coverage and time of each event vs θ |
| `results/fig3_final_maps.png` | thickness maps at the single-bubble time (or t = 150) |
| `results/fig4_late_stage.png` | smallest gap, coverage, porosity vs time |
| `tools/run_case.sh`, `tools/monitor.py` | driver and stop-at-single-bubble monitor |
| `tools/analyze.py`, `tools/rules.py`, `tools/gaps.py`, `tools/plot.py` | analysis, run in this order |

To reproduce:
```
bin/pf3d initOut=ic/N25/ic.bin initOnly=1 r0=4.5 cInit=0.12 nBubbles=25 outDir=ic/N25
tools/run_case.sh 25 1.2 3        # etc.; TCAP=150 python3 tools/monitor.py
cd tools && python3 analyze.py && python3 rules.py && python3 gaps.py && python3 plot.py
```
