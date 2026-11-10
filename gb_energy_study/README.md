# Grain-boundary energy study: bubble growth and coalescence on a planar GB

Modified from `phasefield_gb_bubbles3d.cpp` (conserved version). The physical model, the conservative flux-form Cahn-Hilliard update, the absence of a gas source and of post-step clamping of c, the adaptive time step (safety 0.5), grid, boundary conditions, and the existing columns of `bubble_metrics.csv` are unchanged.

## 1. Analysis of the existing free energy

All three bulk states (grain 1: c,n1,n2 = 0,1,0; grain 2: 0,0,1; bubble: 1,0,0) have f = -1/4, so flat interfaces coexist at mu = 0. Each interfacial energy was obtained as the excess of the full functional, sigma = ∫[A0 (f + 1/4) + gradient terms] dz, over an equilibrium 1-D profile found by gradient flow with fixed end states. This was done on the dx = 1 lattice the code uses, with the GB centred between cells as in the 3-D runs, and also on a fine grid (continuum limit). See `analysis/interface_energies.py` and `analysis/table.py`.

* **Bubble/grain energy gamma_PG.** Only c and the grain's own n vary across this interface. The other order parameter is zero, so the bGB n1²n2² term is identically zero there. gamma_PG depends on A0, bS, j and jGB but **not on bGB**. Values: gamma_PG = 1.068 on the lattice (1.063 to 1.072 depending on interface position) and 1.096 in the continuum.
* **GB energy gamma_GB.** It depends on A0, jGB and bGB. With c = 0 the potential is (n1²+n2²)²/4 - (n1²+n2²)/2 + (bGB-1)/2 n1²n2².
  * bGB = 1: the potential depends only on n1²+n2², so the GB has **no barrier**. Its energy falls as 1/width (0.24, 0.12, 0.06, 0.03 for boxes of 10, 20, 40 and 80 cells) and the GB spreads.
  * bGB < 1: the mixed n1 = n2 state is below either grain, so the **grains are unstable**.
  * The original baseline (bGB = 1) therefore has no well-defined GB energy and cannot be a reference case. bGB = 0.7 (the old "theta = 80") is not a two-grain system at all.
* **Old setTheta().** jGB = 2 bGB tied the gradient coefficient of the grain order parameters to bGB. Because jGB also acts across the bubble/grain interface, every old parameter set changed gamma_PG as well as gamma_GB. The labels 40/50/60/80 do not match Young's-law angles of the model (approx. 124°, 139°, ~180°, undefined).
* **Young's balance.** The bubble meets two identical grains, so the two bubble/grain energies are equal and gamma_GB = 2 gamma_PG cos(theta/2), with theta the dihedral angle measured through the bubble. This is the geometry of the model: a lens symmetric about z = 31.5.
* **Check of the calculation.** At bGB = 2 with jGB = j and bS = 1, the GB potential has exactly the same form as the bubble/grain potential. So gamma_GB = gamma_PG and theta = 120° exactly, which the calculation reproduces (1.0959 vs 1.0959 in the continuum).

## 2. Independent control of GB energy

A0 = 4, bS = 1, j = 2 and jGB = 2 (the baseline values) are held fixed and only bGB is varied. This keeps gamma_PG, the bulk thermodynamics and the bubble interface width exactly unchanged, and changes only gamma_GB, plus the GB width, which is unavoidable.

| case | bGB | gamma_GB (lattice) | GB 10-90% width (cells) | ∫16 n1²n2² dz (GB mobility channel) | Young theta |
|---|---|---|---|---|---|
| low | 1.2 | 0.587 | 5.0 | 10.2 | 148.1° |
| mod. low | 1.5 | 0.848 | 3.5 | 4.9 | 133.2° |
| mod. high | 2.0 | 1.063 | 2.7 | 2.4 | 120.3° |
| high | 3.0 | 1.257 | 2.3 | 1.0 | 107.9° |

bGB is the only changed parameter. It sets the height of the n1/n2 double-well barrier between the grains, i.e. the GB energy. Side effect: a higher barrier gives a narrower GB, which narrows the region where the Mgb mobility term acts by about 10x across the range. This is a kinetic confounder, addressed by the Ms = Mgb = 0 test.

## 3. Initial conditions and the critical radius

### 3.1 cInit = 0.30 is not usable: GB gas film

With cInit = 0.30, gas segregates to the GB and forms a continuous film even when **no bubbles are present**. In a 64³ bubble-free run the GB-core concentration rose from 0.30 to 0.64 by t = 1.5, laterally uniform. A planted r0 = 1.5 bubble (far below any critical size) therefore "grew". Film onset was mapped with exact 1-D columns (NX = NY = 1, same 64-cell reservoir per unit GB area); `preliminary/film`.

| onset time of film (first c ≥ 0.5) | 0.12 | 0.13 | 0.14 | 0.15 | 0.16 | 0.18 | 0.20 | 0.24 |
|---|---|---|---|---|---|---|---|---|
| bGB 1.2 | – | – | – | – | – | – | – | 27 |
| bGB 1.5 | – | – | – | – | – | – | 20 | 4 |
| bGB 2.0 | – | – | – | – | 34 | 9 | 4 | 2 |
| bGB 3.0 | – | 46 | 17 | 10 | 6 | 3 | 2 | 1 |

(– : none by t = 300.) Higher GB energy wets more easily, consistent with a film barrier proportional to 2 gamma_PG - gamma_GB. **cInit = 0.12** is the largest supersaturation that is film-free for all four cases. This is the smallest change that makes the planted bubbles the only gas sinks. Total gas remains conserved; no source is added.

### 3.2 Critical radius at cInit = 0.12

Classical estimate (parallel tangent, relaxed grain order parameter n² = 1 - 2 bS c²): dg_v = 0.445 per unit volume and R* = 2 gamma_PG / dg_v = 4.8 cells. For a GB lens the critical radius of *curvature* is the same, while the lens volume for a given curvature is smaller. Limitations: the interface 10-90% width (2.7 cells) is comparable to R*, so sharp-interface theory is unreliable here; and in a conserved system the supersaturation falls as bubbles grow.

Numerical test: a single bubble on the GB in a 48×48×64 box, bGB = 1.2 (most spherical lens, worst case). Values are cells with c ≥ 0.5 at t = 0.1 → 10 (`preliminary/crit2`):

| r0 | 2.5 | 3.0 | 3.4 | 3.8 | 4.5 | 5.0 |
|---|---|---|---|---|---|---|
| bGB 1.2 | 110 → 0 | 158 → 0 | 208 → 190 (shrinks) | 312 → 624 | 506 → 1252 | 642 → 1692 |
| bGB 3.0 | | | | 312 → 1648 | 506 → 2138 | 650 → 2578 |

The numerical critical radius is between 3.4 and 3.8 cells. r0 = 3.8 would leave less than 12% margin, so **r0 = 4.5** is used for every case. The largest footprint radius (bGB = 3.0) is then about 6 cells, against a minimum centre spacing of 19.0.

### 3.3 Common initial condition

`ic/common_ic.bin` was generated once (seed 11, r0 = 4.5, cInit = 0.12) and is loaded by every run. The bubble centres are byte-identical to the original seed-11 layout (`ic/common_ic.bin.centres.csv` vs `ic/legacy_centres.csv`). Initial gas = 139952.2624, 33 bubbles, minimum spacing 19.02 cells. The initial coverage and gas inventory are identical in all cases (first row of each `bubble_metrics.csv`).

## 4. Code changes

See the header of `src/phasefield_gb_bubbles3d.cpp`:
* bGB, bS, j, jGB and A0 are independent inputs; `theta=` is kept only as a legacy option and prints a warning; a warning is printed for bGB ≤ 1.
* `initOut`/`initIn`/`initOnly`/`centresIn` handle the shared initial condition.
* `outDir`, `params.txt`, `snapshots.csv` (snapshot time ↔ step ↔ file) and `grains_step_N.raw`.
* Appended metrics columns: `free_energy, gb_z_mean, gb_z_std, gb_delta, n_min, n_max, nonfinite`.
* Bit-exact checkpoint/restart (`checkpoint.bin`, `restartIn=`), verified byte-identical against an uninterrupted run.
* Compile-time grid override, used only for the test runs.

## 5. Reproducing

```
bash run_all.sh                          # build, IC, all stages (resumable)
python3 measure_contact_angles.py        # angles -> results/case_properties.json
python3 plot_gb_energy_study.py          # figures A-G + results/combined_metrics.csv
python3 plot_morphology.py times=0,9,21,39
```
Each case directory contains `run.sh` (exact commands, including restarts), `params.txt`, `run.log`, `diagnostics.txt`, `bubble_metrics.csv` and `snapshots.csv`. The `.raw` and `.bin` field files are not in the zip because of their size; `run_all.sh` regenerates them.

## 6. Validation results

### 6.1 Contact angle (measured, not prescribed)
2-D (x-z) lens on the GB, same kernel and lattice (`-DGRID_NY=1`), gas-conserving, equilibrated to t = 400 from two opposite starts: a circle (theta_init = 180°) and a flat lens (theta_init = 80°). The angle is measured by A, 4·atan(h/a), and B, a circle fit to each cap. The circle starts converge by t ≈ 250. The flat starts level off 8-11° lower (1.5 and 3.0) or are still rising slowly (1.2 and 2.0). This is contact-angle hysteresis from triple-line pinning on a lattice whose interfaces are only 1-2 cells wide. The final free energies of the two starts differ by about 1e-5 relative, which is within the effect of their 1-unit gas difference, so neither start is clearly the global minimum. The reported angle is the mean of the four final values ± half their range.

| bGB | gamma_GB | Young theta | circle A / B | flat A / B | **theta measured** |
|---|---|---|---|---|---|
| 1.2 | 0.587 | 148.1° | 153.9 / 151.3 | 145.3 / 141.3 | **148.0 ± 6.3°** |
| 1.5 | 0.848 | 133.2° | 134.9 / 132.4 | 125.3 / 121.2 | **128.4 ± 6.9°** |
| 2.0 | 1.063 | 120.3° | 125.9 / 123.8 | 117.0 / 113.7 | **120.1 ± 6.1°** |
| 3.0 | 1.257 | 107.9° | 110.0 / 108.0 | 100.3 / 97.2 | **103.9 ± 6.4°** |

Young's prediction lies inside the measured bracket in every case. Increasing gamma_GB decreases theta. Figures: `results/F_angle_vs_gb_energy.png`, `results/contact_angle_lenses.png`.

### 6.2 Conservation and numerics (main runs, t = 0-20)
| bGB | max abs mass drift | c range | F increases | dt mean | GB z | max bubbles |
|---|---|---|---|---|---|---|
| 1.2 | 1.0e-13 | [0.073, 1.040] | 0 / 200 | 4.00e-4 | 31.50 | 33 |
| 1.5 | 7.2e-14 | [0.054, 1.036] | 0 / 200 | 4.31e-4 | 31.50 | 33 |
| 2.0 | 8.9e-14 | [0.035, 1.033] | 0 / 200 | 4.62e-4 | 31.50 | 33 |
| 3.0 | 8.1e-14 | [0.024, 1.030] | 0 / 200 | 4.89e-4 | 31.50 | 33 |

dt is about half the stability limit in every case. No new bubbles nucleated (the count never exceeds 33). Away from bubbles the GB-plane concentration is 0.10-0.18, close to the matrix (0.07-0.12), so there is no GB film. Thin bridges seen in top view are necks between flat lens edges. A checkpoint restart (used for case 1.5 after an interruption) reproduced the interrupted trajectory exactly.

## 7. Results, t = 0-20

| bGB | theta meas. | coverage t=5 / 10 / 20 | bubbles t=20 | A1/Atot t=20 | porosity t=20 | first merger: t, coverage | A1/Atot>0.5: t, coverage |
|---|---|---|---|---|---|---|---|
| 1.2 | 148° | 0.206 / 0.258 / 0.342 | 33 | 0.036 | 0.0449 | none | none |
| 1.5 | 128° | 0.239 / 0.313 / 0.431 | 26 | 0.131 | 0.0529 | 16.9, 0.394 | none |
| 2.0 | 120° | 0.273 / 0.370 / 0.549 | 6 | 0.757 | 0.0629 | 11.3, 0.392 | 16.7, 0.496 |
| 3.0 | 104° | 0.309 / 0.452 / 0.654 | 2 | 0.938 | 0.0721 | 8.0, 0.390 | 10.9, 0.479 |

**Observed.** Higher gamma_GB gives a smaller dihedral angle, flatter lenses (`results/morphology_sideview.png`), faster coverage growth, faster gas uptake (porosity) and earlier coalescence in time. Plotted against coverage, the A1/Atot curves nearly coincide: the first merger occurs at coverage 0.390-0.394 in every case that coalesces, and the rise past 0.5 is at 0.48-0.50.

**Interpretation (hypotheses).**
1. With identical bubble centres, coalescence is set mainly by when footprints touch, i.e. by coverage. GB energy acts mainly through how fast coverage grows.
2. For a given gas volume a lower angle gives a larger footprint, more interface area and a larger radius of curvature. A larger radius means a lower Gibbs-Thomson chemical potential in the bubble, hence a larger driving force for uptake. This is consistent with the faster porosity growth at high gamma_GB. It is an equilibrium-wetting effect acting through kinetics.
3. The fast-diffusion channel along the GB (∫16 n1²n2² dz) becomes about 10x narrower from bGB 1.2 to 3.0, yet coverage grows faster. Enhanced GB transport is therefore not the cause of the trend; the Ms = Mgb = 0 runs test this directly.
4. bGB = 1.2 does not reach coverage 0.39 by t = 20, so whether its coalescence also starts there requires the t = 40 extension.

## 8. Mobility sensitivity (Ms = Mgb = 0, Mb = 5; bGB = 1.2 and 3.0; t = 0-20)

| case | t | coverage base / M0 | bubbles | A1/Atot | porosity |
|---|---|---|---|---|---|
| bGB 1.2 | 10 | 0.258 / 0.256 | 33 / 33 | 0.033 / 0.033 | 0.0301 / 0.0298 |
| bGB 1.2 | 20 | 0.342 / 0.340 | 33 / 33 | 0.036 / 0.036 | 0.0449 / 0.0445 |
| bGB 3.0 | 10 | 0.452 / 0.450 | 15 / 15 | 0.211 / 0.211 | 0.0460 / 0.0459 |
| bGB 3.0 | 20 | 0.654 / 0.653 | 2 / 2 | 0.938 / 0.938 | 0.0721 / 0.0719 |

With bGB = 3.0, the first merger is at t = 8.0 and coverage 0.389 with Ms = Mgb = 0, against t = 8.0 and coverage 0.390 at baseline. Removing the enhanced surface and GB mobility changes every metric by less than 1%. Transport is dominated by bulk diffusion (Mb = 5; the Mgb term adds at most 0.5 locally). The differences between the four GB-energy cases therefore come from the equilibrium wetting geometry, not from the change in GB width or the GB diffusion channel. Mass drift in the M0 runs is at most 1e-13. Figure: `results/G_mobility_sensitivity.png`.

## 9. Final results, t = 0-40 (all four cases extended by bit-exact restart)

| bGB | theta meas. | coverage t=10 / 20 / 30 / 40 | bubbles t=40 | A1/Atot t=40 | porosity t=40 | first merger (t, cov) | A1/Atot>0.5 (t, cov) | single network (t, cov) |
|---|---|---|---|---|---|---|---|---|
| 1.2 | 148° | 0.258 / 0.342 / 0.410 / 0.473 | 21 | 0.163 | 0.0672 | 27.6, 0.394 | – | – |
| 1.5 | 128° | 0.313 / 0.431 / 0.543 / 0.600 | 7 | 0.738 | 0.0822 | 16.9, 0.394 | 28.2, 0.527 | – |
| 2.0 | 120° | 0.370 / 0.549 / 0.647 / 0.700 | 2 | 0.939 | 0.0931 | 11.3, 0.392 | 16.7, 0.496 | – |
| 3.0 | 104° | 0.452 / 0.654 / 0.758 / 0.806 | 1 | 1.000 | 0.1026 | 8.0, 0.390 | 10.9, 0.479 | 24.1, 0.704 |

Diagnostics over 0-40: max abs mass drift ≤ 1.33e-13 in all cases, 0 free-energy increases in 400 outputs, bubble count never above 33, GB position 31.50.

**Observed (final).**
* Increasing gamma_GB lowers the measured dihedral angle (148° → 104°), flattens the lenses, raises the rate of coverage and porosity growth, and brings coalescence earlier in time (first merger at t = 27.6, 16.9, 11.3, 8.0).
* The coverage at which coalescence *begins* is the same in all four cases (0.390-0.394).
* After onset the A1/Atot-vs-coverage curves separate systematically. At equal coverage the network is more connected for lower theta: A1/Atot passes 0.5 at coverage 0.479 (bGB 3.0), 0.496 (2.0) and 0.527 (1.5). bGB 1.2 reaches only 0.16 at coverage 0.47.

**Interpretation (hypotheses, not separately tested).**
* Onset is geometric: with identical centres, the first footprints touch at about the same coverage regardless of shape.
* The post-onset ordering is consistent with lens shape. A merged pair of flat lenses (low theta) spreads its gas over a wider, thinner footprint that reaches further toward the next neighbour, so each merger adds more projected connectivity. Rounder lenses (high theta) contract toward a compact shape after merging. Kinetics may also contribute: at low theta the system passes through each coverage faster, leaving less time for the shape to relax.
* The mobility test (section 8) rules out the GB/surface mobility terms as the cause. The differences are attributable to the equilibrium wetting balance gamma_GB = 2 gamma_PG cos(theta/2), acting through bubble shape.
