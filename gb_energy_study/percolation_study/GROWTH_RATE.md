# Is coverage vs time a straight line, and is the rate / acceleration the same for two angles?

The two angles compared are θ = 104° (bGB 3.0) and θ = 148° (bGB 1.2). Every other parameter and the initial condition are identical. The data are the same N = 33 and N = 25 runs as in REPORT.md: `runs/*/bubble_metrics.csv`, coverage every 0.1. No new simulation was needed.

**Script:** `tools/growth_rate.py`. **Figure:** `results/growth_rate/growth_rate.png`. **Tables:** `results/growth_rate/*.csv`.

**Fit window.** The fit runs from t = 1 to the last output before the first merger, i.e. the stage where every bubble grows on its own. After the first merger, coverage jumps each time a neck fills, so a single line cannot describe it. The first time unit is excluded because the initial sharp interfaces relax then. In that unit coverage jumps from 0.127 to 0.138 within 0.1 time units in both cases (0.1375 vs 0.1379), so the jump is not growth.

**Local rate.** dC/dt is the slope of a straight line fitted over a sliding window 1.0 wide (11 outputs). At 148° it wobbles by about ±15 % because coverage counts whole columns. Only trends larger than that wobble are interpreted.

## 1. Straight-line fits (t = 1 → first merger)

| N | θ | window | slope (coverage per time) | R² | max \|data − line\| | local dC/dt at start → end of window |
|---|---|---|---|---|---|---|
| 33 | 104° | 1–7.9 | 0.0305 | 0.9960 | 0.0091 | 0.036 → 0.026 (−29 %) |
| 33 | 148° | 1–27.5 | 0.0088 | 0.9913 | 0.0161 | 0.015 → 0.006 (−61 %) |
| 25 | 104° | 1–7.8 | 0.0272 | 0.9986 | 0.0046 | 0.029 → 0.023 (−19 %) |
| 25 | 148° | 1–26.6 | 0.0088 | 0.9963 | 0.0098 | 0.011 → 0.006 (−49 %) |

**It is not a straight line.** R² is high, but the difference between the data and the line is not random. It is a systematic arch (figure, second column): above the line in the middle, below it at both ends. That is the signature of a curve that bends downward, i.e. a growth rate that keeps falling.

The deviations reach 0.005–0.016 in coverage. The smallest coverage step the code can resolve is 1/128² = 6.1e-5, so the deviations are 75–260 times larger than that. Shorter stretches look almost straight (148°, t = 3–12: R² = 0.9991), but over the whole independent-growth stage both angles decelerate.

## 2. Is the rate the same? No.

Fitting both angles over the **same** time window (both still before any merger):

| N | window | rate 104° | rate 148° | ratio 104°/148° |
|---|---|---|---|---|
| 33 | 1–7.8 | 0.0305 | 0.0109 | 2.79 |
| 33 | 2–7.8 | 0.0295 | 0.0105 | 2.80 |
| 33 | 3–7.8 | 0.0282 | 0.0102 | 2.77 |
| 25 | 1–7.8 | 0.0272 | 0.0097 | 2.81 |
| 25 | 2–7.8 | 0.0267 | 0.0095 | 2.80 |
| 25 | 3–7.8 | 0.0261 | 0.0095 | 2.76 |

The lower angle (higher GB energy) covers the boundary **2.8 times faster**. The ratio barely changes with the window (2.76–2.81) and is the same for N = 33 and N = 25.

## 3. Is the acceleration the same? No.

The acceleration was taken from a quadratic fit C = a + bt + ct², as 2c, in units of coverage per time².

| N | θ | own window (t = 1 → merger) | same window t = 1–7.8 | same window t = 3–7.8 |
|---|---|---|---|---|
| 33 | 104° | −0.0021 | −0.0021 | −0.0010 |
| 33 | 148° | −0.00024 | −0.0006 | +0.0003 |
| 25 | 104° | −0.0011 | −0.0011 | −0.0011 |
| 25 | 148° | −0.00016 | −0.0001 | +0.0005 |

* At 104° the acceleration is negative in every window: the growth clearly slows down.
* At 148°, over the same short window, the acceleration is close to zero and its sign flips with the window. Its slowdown is too small to measure over 7 time units. Over its own longer window (t = 1 → 27) it does slow down: the average rate is 0.0103 on t = 3–12 and 0.0075 on t = 12–26.5 (−27 %).
* Per unit time, the 104° case slows down 7–9 times faster (own windows: N = 33, 0.0021 vs 0.00024; N = 25, 0.0011 vs 0.00016).

## 4. The same comparison at equal coverage instead of equal time

Fourth column of the figure, and `rate_at_equal_coverage.csv`:

* At equal coverage, the 104° rate is 3.0–3.3 times the 148° rate (mean over coverage 0.20–0.30: 3.29 for N = 33, 3.04 for N = 25).
* Multiplying the 148° rate curve by this single number places it on top of the 104° curve, within the 148° wobble. This holds over coverage 0.20–0.32 (N = 33) and 0.15–0.30 (N = 25).
* Both rates fall by about the same fraction as coverage rises. From coverage 0.18 to 0.30, the 104° rate drops 19 % and the 148° rate 20 % (N = 33); for N = 25 the drops are 26 % and 28 %.
* The exception is the stretch just before the 148° first merger (coverage above 0.33, N = 33). There the 148° rate drops further and the ratio rises to 4.1–4.3.

**What this means.** For most of the independent-growth stage, the two angles follow approximately the same rate-vs-coverage curve, and the 104° system simply moves through it about 3 times faster. If that were exact, the deceleration per unit time would differ by the square of the speed-up, about 9–11 times. The measured 7–9 times is close to that, and it is why the accelerations in section 3 are so different even though the curve shapes are similar.

## 5. Possible reason for the bending (hypothesis, not tested here)

The gas is conserved and has no source, so every bubble that grows lowers the supersaturation available to all of them. Porosity keeps rising throughout these windows. For an isolated lens growing by diffusion at constant supersaturation, the radius grows as √t, so its footprint and the coverage would grow linearly in t. Depletion of the supersaturation would then make the rate fall. This is consistent with the rate falling by the same fraction per unit of coverage at both angles, but it was not checked separately (for example by running with a constant-supersaturation reservoir).

## 6. Answer

1. **Coverage vs time is not a straight line** for either angle. Both curves bend downward (the growth rate falls) through the independent-growth stage. Only short stretches are close to linear.
2. **The rate is not the same.** θ = 104° grows 2.8 times faster than θ = 148° over the same time window, and about 3.0–3.3 times faster at equal coverage. This holds for both N = 33 and N = 25.
3. **The acceleration is not the same.** Both are negative. Per unit time, θ = 104° slows down 7–9 times faster, about what is expected if both follow the same rate-vs-coverage curve and 104° runs through it about 3 times faster.
