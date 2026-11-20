# Flat GB migrating (Δf = 0.1) while 33 bubbles on it grow

Question: a flat GB at z = 20 is driven upward by a bulk energy difference Δf = 0.1 while 33 gas bubbles sitting on it grow and coalesce. Does the bubble population pin the GB, or does the GB drag the bubbles? How does coalescence (A1/Atot vs coverage) change?

## Setup

| | |
|---|---|
| Solver | `src/pf_gb_bubbles_mig.cpp` = Flat_Gb's `phasefield_gb_migration3d.cpp` plus items 9–11 in its header |
| Grid | 128 × 128 × 64, periodic x, y; no-flux z walls; top wall at z = 63 |
| Model parameters | A0 = 4, bS = 1, j = jGB = 2, Ms = 0.05, Mgb = 0.5, Mb = 5, LGB = 1 (unchanged) |
| Contact angle | bGB = 2.0 (θ = 120°), the same reference case as Flat_Gb |
| Bubbles | 33 bubbles, same xy centres as `gb_energy_study/percolation_study/ic/N33` (checked with `diff`); r0 = 4.5, seed 11 |
| Initial height | GB and bubble centres at z = 20 |
| Gas | **cInit = 0.18** (see below) |
| Initial condition | One shared file, `ic/N33_z20_c0.18.bin`, so both bubble runs start from the same fields |

Runs (`runs/`):

| Run | Content |
|---|---|
| `bubbles_only` | 33 bubbles, Δf = 0 (GB cannot migrate), t = 0–120 |
| `gb_only_df0.1` | no bubbles, Δf = 0.1, t = 0–400. 32 × 32 × 64: the GB is exactly planar, so the xy size is irrelevant |
| `coupled_df0.1` | 33 bubbles, Δf = 0.1, t = 0–350 |

### Two changes from earlier setups, and why

1. **Bias form (biasMode = 1).** Flat_Gb used h = n1²/(n1²+n2²). That function is undefined at n1 = n2 = 0, which is exactly the bubble interior. There dh/dn ~ 1/|n| and d²h/dn² ~ 1/n², so the explicit update becomes singular and dt → 0.
   - New form: F += Δf · g(n2) with g = n²(3 − 2n).
   - It vanishes with zero slope at n = 0 and at n = 1. Grain I, grain II and the bubble therefore all stay exact equilibrium states, and f_II − f_I = Δf exactly.
   - The bubble interior gets the energy of grain I, the growing grain. A bubble left behind in grain I would therefore feel the original, unbiased model.
   - With Δf = 0 the code is byte-identical to the original solver: the checkpoint after 632 steps from the N33 initial condition matches.
2. **More gas.** With cInit = 0.12 the earlier N = 33 runs reach at most 0.72 coverage (θ = 120°), so 0.8–0.9 is impossible. Pilots on a 64² quadrant with 8 bubbles:
   - cInit = 0.16: coverage 0.80 only at t ≈ 44;
   - cInit = 0.20: coverage reaches 1.0 by t = 49;
   - chosen: cInit = 0.18.

## Check: coverage 0.8–0.9 happens long before the free GB reaches the top

| Independent run | Event | time |
|---|---|---|
| bubbles only | coverage 0.8 | **16.3** |
| bubbles only | coverage 0.9 | **31.8** |
| bubbles only | maximum coverage 0.974 | t = 118 |
| GB only | z = 40 | 130 |
| GB only | z = 60 | **260** |
| GB only | GB leaves the box | 272 |

The free GB velocity is 0.154 cells per time unit, close to Flat_Gb's Moelans-bias value of 0.158. Coverage 0.9 is reached about 8× earlier than the GB reaches the top. ✔

## Results (`results/`)

* `fig1_gb_position.png`: GB position vs time
* `fig2_largest_fraction_vs_coverage.png`: A1/Atot vs coverage
* `fig3_xz_slices.png`: xz slices and top views (supporting)
* `summary.json`: all numbers

**1. Coverage < ~0.94 (t ≲ 40): no pinning; the GB drags the bubbles.**
- The GB moves at 0.135 cells per time unit, against 0.154 for the free GB over the same window, i.e. 12 % slower. It reaches z = 25 at t = 33.8, against about t = 33 for the free GB.
- In no column and at no time does the GB break away from a bubble: `cols_gb_above_bubble` = 0 throughout. The bubbles are carried along instead:
  - their tops advance with the GB (mean top 22.9 → 29.1 by t = 40);
  - their bottoms lag (17.2 → 18.8);
  - the bubble centroid rises 20 → 23.9.
- Gas moves through the matrix (bulk mobility Mb = 5), so a bubble can relocate faster than the GB moves.

**2. Coverage 0.94 → 1 (t ≈ 40–100): the GB slows down.**
- Only small GB necks remain between holes in the gas sheet. They fall behind the free GB: 28.6 against 33.8 at t = 90.

**3. t ≥ 100: coverage = 1; the GB no longer exists.**
- The bubbles have merged into one continuous gas sheet (about 10 cells thick) that separates grain I from grain II.
- Without the bias the same bubbles stop at 0.974. Moving the GB removes the last GB necks.
- From then on the only way to consume grain II is to transfer gas across the sheet. The sheet's upper face (the lower edge of grain II) creeps at **0.027 cells per time unit, about 6× slower than the free GB**. It is at z = 40 at t = 350, while the free GB left the box at t = 272.
- The creep is visibly cell-by-cell. This is the same lattice-pinning effect Flat_Gb found below Δf ≈ 0.03, because the effective driving force is now small.
- **In this sense the system ends up pinned. The mechanism is not Zener-type pinning by isolated bubbles, which never happened. It is the GB being replaced by a continuous gas layer.**

**4. A1/Atot vs coverage is essentially unchanged by GB motion.**
- Both runs follow the same staircase: first mergers at coverage ≈ 0.47; single bubble at coverage ≈ 0.75 (t = 13.0 bubbles only, t = 11.4 coupled).
- The moving GB only speeds up coverage growth: 0.9 at t = 25.3 instead of 31.8.

## Verification

- Gas conservation: |mass drift| ≤ 5e-14 in all runs. Free energy never increased (0 increases in 1200, 3500 and 4000 outputs). No non-finite values.
- The coupled run was stopped at t = 169 and resumed with 12 threads from the checkpoint. The restart was not separately re-verified here. In the percolation study the same checkpoint mechanism reproduced the rows after a restart bit for bit.
- Minimum c over the run is −0.028 (coupled) and −0.009 (bubbles only). This is ordinary Cahn–Hilliard undershoot in the depleted matrix; c is not clamped.

## Limitations

- One parameter set (θ = 120°, Δf = 0.1, cInit = 0.18) and one bubble arrangement.
- The complete sheet (coverage = 1) is a consequence of the gas amount; with cInit = 0.12 coverage would stay ≤ 0.72.
- In the bias choice the bubble is given grain I's energy. Converting grain II into gas therefore gains Δf, the same as GB migration. This is physical for a stored-energy driving force, but it helps the bubbles follow the GB. A different, still continuous choice would only shift the gas solubility. The bias cannot be continuous and keep a "grain II-like" bubble at the same time.
- GB height in the bubbles-only run is measured only over the shrinking set of bubble-free columns. The small rise from 20 to 21.5 reflects the remaining GB necks, not migration.

## Reproduce

```
g++ -O3 -march=native -fopenmp -std=c++17 src/pf_gb_bubbles_mig.cpp -o src/pf_mig_128
g++ -O3 -march=native -fopenmp -std=c++17 -DGRID_NX=32 -DGRID_NY=32 -DGRID_NZ=64 src/pf_gb_bubbles_mig.cpp -o src/pf_mig_32
# exact commands (incl. the restart) are in runs/*/run.sh ; tools/run_case.sh resumes from checkpoint.bin
python3 tools/analyze.py
```

---

# Part 2: round boundary — a spherical grain shrinking inside another grain

Question: replace the flat GB by a closed spherical one. An inner grain (grain II, radius R0 = 32) sits in the centre of an outer grain (grain I). It shrinks by curvature alone (Δf = 0) or by curvature plus the bulk driving force Δf, until it disappears. First without bubbles, then with 25 bubbles placed on the spherical GB.

## Setup

| | |
|---|---|
| Grid | 96 × 96 × 96 (build `src/pf_mig_96`), periodic x, y, no-flux z; box centre (47.5, 47.5, 47.5) |
| Inner grain | n2 = 1 inside R0 = 32, tanh profile; inner grain is the high-energy grain (bias F += Δf·g(n2), biasMode = 1, as in Part 1) |
| Δf | 0, 0.025, 0.05, 0.1 |
| Other parameters | identical to Part 1 (A0 = 4, bS = 1, bGB = 2.0 → θ = 120°, j = jGB = 2, Ms = 0.05, Mgb = 0.5, Mb = 5, LGB = 1) |
| Bubbles | 25 bubbles with centres exactly on the sphere r = 32, minimum chord 19.07 (same areal density as the 33 bubbles on the 128² flat GB within 3.5 %), r0 = 4.5, cInit = 0.18, seed 11; one shared IC `ic/sphere_R32_N25_c0.18.bin` for all Δf |
| Stopping | 10 time units after the inner grain has vanished (`stopAfterVanish=10`) |

Measures written to `runs/sphere_*/sphere.csv` every 0.1 time units:

* `r_gb_mean`: GB radius = distance from the centre to the n1 = n2 crossing, averaged over 6000 uniformly distributed rays (only rays where a real GB exists, n1 + n2 > 0.5 at the crossing).
* `R_vol` = (3 V_inner / 4π)^{1/3} with V_inner = number of cells with n2 > n1 and n1 + n2 > 0.5: the volume-equivalent radius of the inner grain. Without bubbles the two measures agree within 0.14 cells at all times. With bubbles `R_vol` is the radius of the grain that is left inside the gas, which is the quantity that matters once the GB is covered.
* `coverage`: fraction of rays that cross ≥ 3 cells of gas (c ≥ 0.5) = fraction of the spherical GB covered by gas; patches by 6-neighbour connectivity on the sphere.

## Results: no bubbles (`results/sphere_fig0_nobubbles_radius_vs_time.png`, `sphere_nobubbles_radius_vs_time.csv`)

| Δf | grain vanishes at t | sharp-interface ODE | R = 16 at |
|---|---|---|---|
| 0 | 136 | 128 | 102.5 |
| 0.025 | 110 | 103 | 79.7 |
| 0.05 | 93 | 87 | 65.4 |
| 0.1 | 71 | 67 | 48.4 |

* The curves follow dR/dt = −L κ (2/R + Δf/γ_GB) (γ_GB = 1.0724): at Δf = 0, R² ≈ R0² − 8t with the characteristic acceleration at the end; Δf adds a constant speed, so the grain disappears in half the time at Δf = 0.1.
* The simulation is 6–7 % slower than the sharp-interface prediction in all four cases; this is the same lattice-friction offset as for the flat GB (0.154 measured vs 0.187 predicted at Δf = 0.1).
* The grain stays spherical and centred (centroid offset 0; see the slices in `sphere_fig3_slices.png`).

## Results: 25 bubbles on the spherical GB (`results/sphere_fig1_radius_vs_time.png`, right panel; `sphere_fig2_coverage.png`)

| Δf | single patch at t | coverage 0.8 | coverage 0.9 | shell closed (coverage 1) | grain vanishes at t | without bubbles |
|---|---|---|---|---|---|---|
| 0 | 11.7 | 18.6 | 32.3 | 155 | **228** | 136 |
| 0.025 | 10.9 | 17.7 | 30.6 | 125 | **195** | 110 |
| 0.05 | 10.6 | 16.9 | 29.6 | 109 | **170** | 93 |
| 0.1 | 9.7 | 15.9 | 29.2 | 99 | **135** | 71 |

1. **Bubble growth and coalescence are nearly independent of Δf.** The coverage(t) curves and the A1/Atot-vs-coverage staircases of the four runs lie on top of each other (fig 2): first mergers at coverage ≈ 0.47, a single connected patch at coverage ≈ 0.6–0.68, exactly as on the flat GB. GB motion only shifts the timing by about 10 %.
2. **Early on, the inner grain shrinks faster with bubbles than without.** Up to t ≈ 20 the bubbles grow into both grains, so the volume of the inner grain drops (R_vol = 27.8 at t = 20 for Δf = 0, against 29.6 without bubbles). This is gas consuming the grain, not GB migration: the GB necks between the bubbles (dotted lines) move *slower* than the free GB (R_gb = 30.1 at t = 20).
3. **Then the gas shell slows everything down.** From t ≈ 50 onward the inner grain with bubbles is larger than without, and after the shell closes the inner grain is a solid sphere inside a gas shell 8–10 cells thick. It still shrinks, at Δf = 0 too, because the shell's inner surface is more strongly curved than its outer surface, so moving matrix material from the inner face to the outer face lowers the total surface energy; the gas is redistributed through the shell. Measured shrink rates at the same radius:

   | Δf | R 20 → 15, with shell | without bubbles | ratio | R 15 → 10, with shell | without | ratio |
   |---|---|---|---|---|---|---|
   | 0 | 0.095 | 0.214 | 2.3 | 0.134 | 0.301 | 2.2 |
   | 0.025 | 0.107 | 0.257 | 2.4 | 0.148 | 0.344 | 2.3 |
   | 0.05 | 0.120 | 0.300 | 2.5 | 0.162 | 0.390 | 2.4 |
   | 0.1 | 0.146 | 0.384 | 2.6 | 0.187 | 0.478 | 2.6 |

   The inner grain shrinks 2.2–2.6× slower once it is enclosed in gas, and the lifetime of the grain is 1.7–1.9× longer (228 vs 136 at Δf = 0; 135 vs 71 at Δf = 0.1). The Δf dependence survives: the shell-enclosed grain still shrinks faster with larger Δf.
4. **As on the flat GB, there is no Zener-type pinning by isolated bubbles.** The GB never detaches from the bubbles; the slowdown comes from the GB being replaced by a continuous gas layer. The difference from the flat case is that the shell, being closed, cannot stop the inner grain: capillarity of the shell's inner surface keeps consuming it, with or without Δf.

Verification: |mass drift| ≤ 4.4e-14, 0 free-energy increases (3000–4000 outputs per run), no non-finite values, in all 8 runs. The three bubble runs Δf = 0, 0.025, 0.05 were stopped at t = 147 / 123 / 142 and resumed from their checkpoints with 4 threads (`runs/*/run.sh`); the CSV rows are continuous (0.1 spacing, no gaps or duplicates). The restart itself was not re-verified bit for bit here.

Limitations: one bubble arrangement, one gas amount (cInit = 0.18, chosen so that the flat GB could reach 0.8–0.9 coverage; here it gives a closed shell in every case), one contact angle, one R0. The inner grain is only 3 interface widths across when R ≈ 10, so the last ~20 time units of each curve are not sharp-interface physics.

Reproduce: `g++ -O3 -march=native -fopenmp -std=c++17 -DGRID_NX=96 -DGRID_NY=96 -DGRID_NZ=96 src/pf_gb_bubbles_mig.cpp -o src/pf_mig_96`; commands in `runs/sphere_*/run.sh` (no bubbles: `sphereR=32 nBubbles=0 cInit=0 noise=0 ... stopAfterVanish=5`; bubbles: `initIn=ic/sphere_R32_N25_c0.18.bin sphereR=32 ... stopAfterVanish=10`); `tools/sphere_orchestrate.sh`; `python3 tools/analyze_sphere.py`.

# Part 3: is the coalescence path independent of GB motion? (robustness check, 2026-11-16)

In Part 1 the A1/Atot-vs-coverage curve was the same with and without GB motion, but the GB moved only 0.73 cells while the bubbles coalesced (t = 4.3–8.3), so the test was weak. Here it is repeated for another contact angle, another bubble number, and with a GB forced to move 3–10× faster. Each case is a pair of 128 × 128 × 64 runs from the same initial condition: GB fixed (Δf = 0) and GB free to migrate. A GB-only run (32 × 32 × 64, no bubbles) gives the free GB velocity.

| Case | θ (bGB) | N | Δf | LGB | free GB v | runs |
|---|---|---|---|---|---|---|
| reference (Part 1) | 120° (2.0) | 33 | 0.1 | 1 | 0.154 | `bubbles_only`, `coupled_df0.1`, `gb_only_df0.1` |
| contact angle changed | 128° (1.5) | 33 | 0.1 | 1 | 0.210 | `theta128_*` |
| bubble number changed | 120° (2.0) | 22 | 0.1 | 1 | 0.154 | `N22_*` (IC `ic/N22_z20_c0.18.bin` = the first 22 of the 33 sites) |
| low contact angle | 104° (3.0) | 33 | 0.3 | 1 | 0.343 | `theta104_*`; at Δf = 0.1 the free GB is lattice-pinned (stops at z = 20.58, `theta104_gb_only_df0.1`) |
| strong test A | 120° | 33 | 0.5 | 1 | 0.784 | `coupled_df0.5`, `gb_only_df0.5` (the bias also changes the gas thermodynamics) |
| strong test B | 120° | 33 | 0.1 | 10 | 1.537 | `L10_*` (same free energy, GB 10× more mobile) |

Kinetics-only spread: `bubbles_only` (LGB = 1) and `L10_bubbles_only` (LGB = 10) are both GB-fixed with identical thermodynamics; their paths differ by at most 0.012 in coverage (median 0.005).

Results (`results/fig4_percolation_robustness.png`, `fig5_all_paths.png`, `fig6_gb_motion_new_cases.png`, `fig7_shift_vs_gb_motion.png`, `robustness_summary.json`). Path shift = coverage(GB moving) − coverage(GB fixed) at the same A1/Atot, median over A1/Atot = 0.10…0.95.

| Case | coverage at A1/Atot = 0.5 (fixed / moving) | GB moved during coalescence (free GB) | path shift, median (25–75 %) | GB speed after coalescence / free | max detached columns |
|---|---|---|---|---|---|
| reference | 0.548 / 0.551 | 0.73 (0.70) | +0.002 (−0.003, +0.007) | 0.93 | 0 |
| θ = 128° | 0.553 / 0.542 | 1.13 (1.16) | −0.007 (−0.010, 0.000) | 0.75 | 0.2 % |
| N = 22 | 0.573 / 0.578 | 1.47 (1.21) | +0.005 (+0.003, +0.007) | 0.94 | 0 |
| θ = 104°, Δf = 0.3 | 0.559 / 0.575 | 0.47 (0.38) | +0.020 (+0.009, +0.025) | 1.03 | 0 |
| Δf = 0.5 | 0.548 / 0.567 | 2.04 (1.91) | +0.028 (+0.020, +0.033) | 0.81 | 18.7 % |
| LGB = 10 | 0.554 / 0.538 | 2.49 (7.99) | −0.020 (−0.020, −0.018) | 0.15 | 1.6 % |

1. **Both requested cases confirm the result.** θ = 128° and N = 22 shift the path by −0.007 and +0.005, inside the kinetics-only spread, although their GB moved 1.5–2× farther during coalescence than in the reference. The bubble number moves the whole curve (first merger at coverage 0.35 with 22 bubbles instead of 0.47); GB motion does not.
2. **Forcing the GB 3–10× faster gives a small real shift of about 0.02 in coverage**, twice the kinetics-only spread. Its sign follows the mechanism, not the displacement: a more mobile GB (LGB = 10) shifts the path to lower coverage, a larger bias (Δf = 0.3, 0.5) to higher coverage. θ = 104° at Δf = 0.3 shifts by +0.020 with only 0.47 cells of GB motion. At Δf = 0.5 five extra gas pockets appeared at t ≈ 3, and porosity at t = 4 was 0.068 against 0.046 with the GB fixed, so the bias changes the gas thermodynamics.
3. **The percolation point (A1/Atot = 0.5) stays at coverage 0.54–0.58 in every run.**
4. **The bubbles drag a fast GB hard.** With LGB = 10 the GB covers only 31 % of the free distance during coalescence and runs at 15 % of its free speed afterwards (0.225 vs 1.537), barely faster than with LGB = 1 (0.144). It stays attached (≤ 1.6 % of columns detached), so the coupled speed is set by how fast the bubbles can follow.
5. **Breakaway appears only at Δf = 0.5**, after coalescence: from t ≈ 10 the GB leaves the bubble layer behind, reaching 18.7 % of covered columns by t = 30 (detached = covered column with the real GB above the bubble top).

Verification: |mass drift| ≤ 1.5e-13, 0 free-energy increases, no non-finite values, continuous CSVs (0.1 spacing) in all 13 runs, including the restarted `theta128_*`, `N22_*`, `theta104_*`. Runs were shortened (t = 20–60); every one reached a single bubble before its end.

Limitations: one random seed and one site pattern (N = 22 is a subset of the 33), one gas content; the "single bubble" coverage is sensitive to the last small straggler (N = 22: 0.646 vs 0.781) and is not a path measure.

Reproduce: IC `./src/pf_mig_128 nBubbles=22 r0=4.5 cInit=0.18 noise=0.001 seed=11 minSep=19 zgb0=20 bGB=2.0 initOut=ic/N22_z20_c0.18.bin initOnly=1 outDir=ic`; runs `tools/run_case.sh NAME THREADS ARGS` with the arguments in `runs/NAME/run.sh` (e.g. strong test B: `initIn=ic/N33_z20_c0.18.bin bGB=2.0 dfBias=0.1 biasMode=1 LGB=10 ... tEnd=25`); `tools/robustness_orchestrate.sh`; `python3 tools/analyze_robustness.py`.
