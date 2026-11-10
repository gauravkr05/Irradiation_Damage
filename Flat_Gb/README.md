# Grain-boundary migration driven by unequal bulk free energies (bubble-free)

This is a separate experiment. Nothing in the original project (`gb_energy_study/`) was modified. There are no bubbles, no gas, and no bubble or coverage analysis: c = 0 everywhere at all times, and the max |c| recorded in every output is exactly 0.

## 1. The existing model

Source: `src/phasefield_gb_bubbles3d.cpp` (original).

**Free energy**

F = ∫ [A0 f + j/2 |∇c|² + jGB/2 (|∇n1|² + |∇n2|²)] dV

f = c⁴/4 − c²/2 + n1⁴/4 − n1²/2 + n2⁴/4 − n2²/2 + bGB/2 n1²n2² + bS c²(n1² + n2²)

**Evolution equations**
- Grains (explicit Allen–Cahn): ∂n_i/∂t = −LGB δF/δn_i, with n_i clamped to [0, 1] after each step.
- Gas (conserved Cahn–Hilliard): ∂c/∂t = ∇·(M∇μ), with μ = A0[c³ − c + 2bS c(n1² + n2²)] − j∇²c.

**Boundary conditions and time step**
- Periodic in x and y; no-flux (mirror) walls in z.
- dt = min(dtMax, 0.5·dt_CH, 0.5·dt_AC, time to next output).

**Why the original system has no driving force.** At c = 0 both pure grains have f = −1/4, i.e. A0·f = −1 per cell. The bGB term is symmetric in n1 and n2. A flat boundary therefore feels no net force; only curvature could move it.

**Orientation in the original code.** Grain I (n1) lies below the GB and Grain II (n2) above.

## 2. Modification: a controlled bulk-energy difference

File: `src/phasefield_gb_migration3d.cpp`; changes are documented in its header, item 7.

F → F + ∫ Δf · h(n1, n2) dV,   h = n1² / (n1² + n2²)   (parameter `dfBias` = Δf)

* **Pure states.** h = 1 in Grain I and 0 in Grain II, so **f_I − f_II = Δf** exactly. The bias is added once per cell, in the same units as A0·f.
* **The grain values are unchanged.** ∂h/∂n1 = 2n1n2²/s² and ∂h/∂n2 = −2n1²n2/s² (s = n1² + n2²) both vanish at the pure states, so n = 1 and n = 0 remain exact equilibria. This is the standard multi-order-parameter interpolation (Moelans-type).
* **Stability.** The bias changes the curvature of the potential in grain I by −2Δf, against A0(bGB − 1) = 4, so the grain states stay stable for Δf < 2. The bias curvature is also added to the dt_AC bound.
* **Scope.** h is independent of c, so the bias enters only the n1/n2 equations, the free-energy diagnostic and the time-step bound. With Δf = 0 every code path is the original one.
* **Bubble-free output (`gbOnly=1`).** All bubble and coverage measurements are skipped and `gb_position.csv` is written instead.

## 3. Setup

| Item | Value |
|---|---|
| Grid | 32×32×64, via the existing compile-time override (`-DGRID_NX=32 -DGRID_NY=32 -DGRID_NZ=64`); dx = 1, same discretisation |
| Parameters | bGB = 2.0 (reference case, θ = 120°); A0 = 4, bS = 1, j = 2, jGB = 2, LGB = 1, Ms = 0.05, Mgb = 0.5, Mb = 5, safety = 0.5, dtMax = 1e-3, outDt = 0.1 |
| Initial condition | `make_ic.py`: c = 0; **Grain II (favoured) below, Grain I above**, so the GB can migrate upward into the 44-cell side; tanh profile with w = 1.5 (same as the original generator); n1 = n2 at cell index z = 20 |
| Main bias | Δf = 0.05. That is 5 % of the bulk energy density, and Δf·δ/γ_GB ≈ 0.13, where δ ≈ 2.7 is the GB width. |
| Additional runs | zero-bias control; Δf = 0.025, 0.0375, 0.075, 0.1; Δf = 0.05 with an initial 1-cell sinusoidal GB ripple (planarity test) |
| GB position | In every x-y column, the first z where n1 − n2 changes sign, linearly interpolated. This gives mean, std, min and max, and counts of columns with no crossing or several crossings. All runs have exactly one crossing in every column. |

## 4. Results (`results/`)

GB energy from the control's excess energy: γ_GB = (F − N·f_bulk)/A = **1.0724**.
Continuum (sharp-interface) Allen–Cahn prediction for a flat GB: v = LGB·jGB·Δf/γ_GB = 1.865·Δf.

| Δf | z_GB: start → end (t_end) | Mean velocity v_fit, t ≥ 10 | R² | Measured / predicted | dF/dt vs −Δf·A·v |
|---:|---|---:|---:|---:|---|
| 0 (control) | 20.000 → 20.000 (250) | 0 | – | – | 0 / 0 |
| 0.025 | 20.000 → 20.627 (250), **pinned after t ≈ 60** | 0.00004 | – | ≈ 0 | – |
| 0.0375 | 20.000 → 30.03 (250) | 0.04001 ± 0.00003 | 0.9983 | 0.572 | −1.5357 / −1.5364 |
| **0.05** | **20.000 → 36.88 (250)** | **0.06766 ± 0.00002** | 0.9997 | **0.726** | −3.4655 / −3.4643 |
| 0.075 | 20.000 → 41.73 (190) | 0.11432 ± 0.00003 | 0.9999 | 0.817 | −8.7798 / −8.7799 |
| 0.1 | 20.000 → 42.19 (140) | 0.15799 ± 0.00004 | 0.9999 | 0.847 | −16.1774 / −16.1781 |
| 0.05, rippled | 20.000 → 36.98 (250) | 0.06761 ± 0.00002 | 0.9997 | 0.725 | −3.4628 / −3.4616 |

The velocity fits start at t = 10, after the initial relaxation. The interface shape relaxes by t ≈ 1.1 (control), and the first ~10 time units contain a ~0.6-cell adjustment from the non-equilibrium tanh start.

**Velocity is not constant within a cell.** It oscillates with a period of exactly one cell crossed: about 14.8 time units at Δf = 0.05, where the 1-unit-window velocity ranges 0.032–0.133. The GB width integral oscillates with the same period. Averaged over many cells, the velocity is constant (R² ≥ 0.998).

**Ripple.** Its amplitude decays at a rate of 0.0763, against the Allen–Cahn curvature-flow prediction LGB·jGB·k² = 0.0771. It then migrates like the planar GB.

Figures: `gb_position_vs_time.png` (main), `displacement_velocity_relaxation.png`, `velocity_vs_driving_force.png`, `gb_profiles.png`.
Data: `gb_position_all_runs.csv`, `velocity_summary.csv`, `velocity_history.csv`, `profile_checks.csv`, `summary.json`.

## 5. Validation

* **Original behaviour preserved.** With Δf = 0 the code reproduces the original code bit-for-bit (identical checkpoint and metrics, on a bubble test IC).
* **Bias magnitude.** Uniform domains give f_I = −0.950000 and f_II = −1.000000 per cell at Δf = 0.05, and both pure states stay exactly at n = 1/0.
* **Direction.** The GB moves upward into Grain I for every Δf above the threshold, so the favoured Grain II grows, as predicted by the functional.
* **Energetic consistency.** F decreases at exactly Δf × GB area × v (within 0.05 % for every moving case), so the motion is the bulk-energy-driven migration and not a transient.
* **Stationary control.** The control is stationary to all printed digits.
* **Wall distance.** The GB stays ≥ 20.8 cells below the top wall, and n2 at the top wall stays 0.
* **Planarity.** The boundary is planar throughout: std ≤ 1e-5 cells, one crossing per column. In the rippled run the perturbation decays.
* **Numerical health.** In every run: no free-energy increases; no non-finite values; n stays within [0, 1]; gas exactly 0; dt ≈ 5.2e-4, set by the unchanged Cahn–Hilliard bound; unbroken 0.1 output grid.
* **Resumed run.** The Δf = 0.1 run was interrupted and resumed bit-exactly from its t = 80 checkpoint.

## 6. Interpretation and limitations

* **Observation.** A bulk free-energy difference drives the planar GB into the higher-energy grain, at a velocity that is constant on average and proportional to the free-energy dissipation.
* **Not linear in Δf: lattice pinning.** Below a threshold between Δf = 0.025 and 0.0375 the boundary is pinned. Above it, it moves in a stop-start way, once per grid cell, at 57–85 % of the continuum prediction, approaching it as Δf grows. The diffuse GB is only ~2.7 cells wide (10–90 %), so its energy depends on its position relative to the grid cells. The once-per-cell oscillations in velocity and GB width are direct evidence.
* **Consequence.** This is a discretisation artifact, not continuum physics: GB mobilities measured on this grid are underestimated, and small driving forces are suppressed entirely.
* **Remedies (not applied, as they change the prescribed setup):**
  - a finer grid relative to the interface width;
  - larger gradient coefficients (wider GB, which also changes γ_GB);
  - time-averaging only over whole numbers of cell crossings, which was done here.
* **Other limitations:**
  - one bGB value;
  - no gas, so solute drag or segregation effects are absent by design;
  - z walls are no-flux; the run was stopped well before the GB approached them.

## 7. Commands

```
g++ -O3 -march=native -fopenmp -std=c++17 -DGRID_NX=32 -DGRID_NY=32 -DGRID_NZ=64 \
    src/phasefield_gb_migration3d.cpp -o src/pf_mig_32
python3 make_ic.py out=ic/gb_z20.bin                 # planar GB at z = 20
python3 make_ic.py out=ic/gb_z20_amp1.bin amp=1.0    # rippled GB
bash run_all.sh      # control, 0.025, 0.05, rippled 0.05 (resumable, locked)
bash run_scan.sh     # 0.0375, 0.075, 0.1
python3 analyse.py   # results/
# single run, e.g.:
src/pf_mig_32 initIn=ic/gb_z20.bin bGB=2.0 dfBias=0.05 gbOnly=1 tEnd=250 outDt=0.1 snapDt=10 outDir=runs/df0.05
```
