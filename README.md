# Diabetes Disease Progression: CTMC Calibration and Simulation in R

This repository contains an **R implementation of a continuous-time Markov chain (CTMC) model of diabetes progression**, including calibration of transition intensities using **equilibrium-constrained Maximum Caliber (MaxCal)**, deterministic population trajectories, stochastic individual-level simulation, survival analysis, and sensitivity checks.

The principal script is `DiabetesCtmcToOdeV31_Analysis.R`. It is companion analysis code for the manuscript referenced in the source as `DiabetesCtmcToOdeV31.tex`.

> **What this code is for:** Reproducing the model-based numerical analyses, diagnostic calculations, tables, and figures from the supplied, preprocessed NHANES composition targets and literature-derived transition inputs. It is a research implementation, **not** a tool for individual diagnosis, prognosis, or treatment decisions.

## 1. At a glance

The script performs the following tasks, in order:

1. Reads preprocessed diabetes-state proportions and literature transition information from three CSV files.
2. Converts published mean waiting times or cumulative event probabilities into continuous-time transition rates.
3. Checks whether the literature rates satisfy the chosen stationary population composition.
4. Calibrates positive transition rates using a constrained MaxCal (relative-entropy-rate) objective.
5. Compares the literature rates, an earlier exact-completion benchmark, and the MaxCal-calibrated rates.
6. Computes CTMC population trajectories, phase-type survival quantities, and a 2,000-person Gillespie stochastic simulation cohort.
7. Performs sensitivity, elasticity, and numerical validation analyses.
8. Writes CSV results, publication-style PNG figures, and text summaries of manuscript quantities.

The implementation uses **base R only**; no add-on R packages are required.

## 2. Repository layout

Arrange the repository as follows (the three input CSVs are **required** and must be supplied separately):

```text
.
├── README.md
├── R/
│   └── DiabetesCtmcToOdeV31_Analysis.R
├── Data/
│   ├── DiabetesCtmcToOdeV31_nhanes_average_target.csv
│   ├── DiabetesCtmcToOdeV31_literature_edge_inputs.csv
│   ├── DiabetesCtmcToOdeV31_nhanes_cycle_ratios_digitised.csv
│   └── DiabetesCtmcToOdeV31_derived_*         # generated when run
└── Figures/
    └── DiabetesCtmcToOdeV31_*.png             # generated when run
```

The `Data/` and `Figures/` directories are created if they do not exist, but **the three input CSV files are not created by the script**. The script's file paths are relative to the *current working directory*, so run it from the repository root.

If you keep the `.R` file in the repository root instead of `R/`, use `Rscript DiabetesCtmcToOdeV31_Analysis.R` in the instructions below.

## 3. Requirements and quick start

- Install a recent version of [R](https://cran.r-project.org/), including the `Rscript` command-line executable.
- No `install.packages()` step is needed.
- The figure-writing helper explicitly requests the `cairo` PNG device; check support with `capabilities("cairo")` in R if figure generation fails.
- Place the three input CSV files at the exact `Data/` paths shown above.

From a terminal opened at the **repository root**, run:

```bash
Rscript R/DiabetesCtmcToOdeV31_Analysis.R
```

Alternatively, in R or RStudio, first set the working directory to the repository root and execute:

```r
source("R/DiabetesCtmcToOdeV31_Analysis.R")
```

The script prints progress messages including `Done.`, `Review checks written.`, and `V31 4-s.f. numbers written.` as its later sections complete. Inspect the generated `Data/` files and `Figures/` images for results.

**Important:** Uploading only the `.R` script is not sufficient to reproduce the analysis. A user also needs the three specified CSV input files (subject to any data-sharing restrictions).

## 4. Input data

| Input file in `Data/` | Expected fields | Use |
| --- | --- | --- |
| `DiabetesCtmcToOdeV31_nhanes_average_target.csv` | `state`, `average_proportion` | Target six-state stationary proportions for calibration. |
| `DiabetesCtmcToOdeV31_literature_edge_inputs.csv` | `edge`, `type`, `horizon_years`, `probability` | Literature transition information used to construct reference intensities. |
| `DiabetesCtmcToOdeV31_nhanes_cycle_ratios_digitised.csv` | `cycle_label`, `P`, `T`, `Mi`, `Ma`, `B` | Cycle-specific digitised state-composition series used for the NHANES figure. |

**Input conventions:**

- The **target CSV should list states in the order** `H, P, T, Mi, Ma, B` because the script uses that ordering in matrix calculations; proportions should be nonnegative and have a positive total. The script normalises them to sum to one.
- The literature CSV should contain the transition edges `PT, TP, TMi, TMa, MiT, MiB, MiD, MaB, MaD, BD`, in that order (the `HP` reservoir entry rate is calibrated rather than read from this file).
- In the literature CSV, `type = mean_time` means the reference rate is `1 / horizon_years`; probability-based rows use `-log(1 - probability) / horizon_years`. Horizons are in **years**. Use valid event probabilities and positive horizons.
- The cycle-specific series should have a numerical `cycle_label` (plotted horizontally) and the five disease-state proportion columns shown above. The source script describes these as **digitised from an earlier figure**, not reconstructed from NHANES participant-level records.

The supplied R file does not contain the source CSVs or their original study documentation. For provenance, eligibility definitions, and sampling details, consult the relevant project data documentation and manuscript.

## 5. Model states and transition structure

The model distinguishes the following states:

| Symbol | State |
| --- | --- |
| `H` | Healthy reservoir, used in the **calibration** chain |
| `P` | Prediabetes |
| `T` | Type 2 diabetes without modelled complications |
| `Mi` | Microvascular complications only |
| `Ma` | Macrovascular complications only |
| `B` | Both microvascular and macrovascular complications |
| `D` | Death (absorbing in the **simulation** chain) |

Two related Markov chains are used:

- **Regenerative calibration chain:** states `(H, P, T, Mi, Ma, B)`. Its death-related edges (`MiD`, `MaD`, `BD`) are redirected to `H`, so that a stationary target distribution can be imposed. The `H -> P` rate is determined during calibration.
- **Absorbing progression chain:** states `(P, T, Mi, Ma, B, D)`. Death is absorbing; there is no transition back to a healthy state. This is the chain used for disease trajectories and survival calculations.

Allowed disease-progression transitions are `P -> T`, `T -> P`, `T -> Mi`, `T -> Ma`, `Mi -> T`, `Mi -> B`, `Mi -> D`, `Ma -> B`, `Ma -> D`, and `B -> D`. These are the edges explicitly implemented in the script; other transitions are not included.

## 6. Calibration and mathematical approach

Let `q_e` denote a transition intensity on allowed edge `e`, `q_e^(0)` its literature-derived reference intensity, and `pi_i` the desired stationary population proportion in state `i`. The script minimises the equilibrium-constrained MaxCal objective

$$
J(q\mid q^{(0)}) =
\sum_{e\in E_{\mathrm{ref}}} \pi_{\mathrm{from}(e)}
\left[q_e\log\left(\frac{q_e}{q_e^{(0)}}\right)-q_e+q_e^{(0)}\right],
$$

subject to

$$
\pi^{\mathsf T}Q_{\mathrm{cal}}(q)=0, \qquad q_e\geq 0,
$$

with any scenario-specific fixed rates held at their prescribed values. The reservoir edge `H -> P` is unreferenced in the objective. Interior optimality conditions lead to exponentially tilted literature rates on referenced free edges,

$$
q_e=q_e^{(0)}\exp(u_{\mathrm{to}(e)}-u_{\mathrm{from}(e)}),
$$

where the potentials and reservoir inflow rate are solved by a damped Newton procedure with a positivity safeguard.

The script retains **three rate sets** for comparison:

- `Literature`: literature-derived reference rates.
- `V16exact` (displayed in figures as **Exact completion**): a previously specified positive calibration benchmark hard-coded in the script.
- `MaxCal`: the newly solved relative-entropy-rate calibration.

For a row-vector of cohort state proportions `x(t)`, the deterministic CTMC equation is

$$
\frac{d}{dt}x(t)=x(t)Q,\qquad x(t)=x(0)e^{Qt}.
$$

The script evaluates the matrix exponential using a base-R scaling-and-squaring Taylor routine (`expmTS`), then independently simulates individual state transitions with Gillespie's stochastic simulation algorithm (`ssa_cohort`). No external differential-equation solver is used.

## 7. Simulation settings and derived analyses

**Deterministic trajectories:** Time from **0 to 50 years**, sampled in **0.1-year** increments. The script compares survival, uncomplicated T2D prevalence, and microvascular and macrovascular burden across the rate sets.

**Initial cohorts:**

| Script label | Initial proportions `(P, T, Mi, Ma, B, D)` |
| --- | --- |
| `AllPrediabetes` | `(1, 0, 0, 0, 0, 0)` |
| `DD2baseline` | `(0, 0.65, 0.12, 0.17, 0.06, 0)` |

**Stochastic simulation:** The MaxCal generator is used to simulate **2,000 independent individuals** in each initial cohort. A fixed seed (`42`) supports reproducibility; a further check uses **20 replicate cohorts** with seed `4242`. The figures display exact-semigroup trajectories, stochastic cohort proportions, and **pointwise normal-approximation 95% bands**. These are not simultaneous confidence bands, and the normal approximation can perform poorly for very rare states.

**Survival analysis:** For the transient-state generator `T`, phase-type calculations use `(-T)^{-1}` and `exp(Tt)` to compute expected time to death, expected occupation times, median survival, and model mortality probabilities at **5, 10, 20, and 30 years**.

**Sensitivity analysis:** Scenarios perturb disease-state target proportions by **±10%**, adjust selected mortality or progression priors, fix subsets of transition rates, and apply a joint perturbation. The script also calculates rate elasticities to target proportions by centred finite differences and checks them against differentiation of a KKT system.

**Validation checks:** The code evaluates stationarity and optimality residuals, rate positivity, cycle relationships, a null-space numerical comparison, a renewal occupation identity, rate-scale invariance, covariance evolution, Euler discretisation error, and stochastic-simulation diagnostics. These are **numerical diagnostics** of the implemented model, not external validation on independent patient records.

## 8. Output files

All outputs are written relative to the repository root. Running the script again **overwrites** output files with the same names.

### Tables and numeric summaries (`Data/`)

| Generated file | Contents |
| --- | --- |
| `DiabetesCtmcToOdeV31_derived_calibrated_rates.csv` | Reference, exact-completion, and MaxCal edge rates and their changes. |
| `DiabetesCtmcToOdeV31_derived_holding_times.csv` | Mean holding times by disease state and rate set. |
| `DiabetesCtmcToOdeV31_derived_survival_summaries.csv` | Survival, mortality, and state occupation summaries by generator and starting cohort. |
| `DiabetesCtmcToOdeV31_derived_trajectories_MaxCal_AllPrediabetes.csv` | Exact and simulated time trajectories for the all-prediabetes cohort. |
| `DiabetesCtmcToOdeV31_derived_trajectories_MaxCal_DD2baseline.csv` | Exact and simulated time trajectories for the second baseline cohort. |
| `DiabetesCtmcToOdeV31_derived_ssa_diagnostics.csv` | SSA-versus-exact diagnostic measures, including the Euler comparison. |
| `DiabetesCtmcToOdeV31_derived_sensitivity.csv` | Recalibrated transition rates and survival measures for sensitivity scenarios. |
| `DiabetesCtmcToOdeV31_derived_elasticities.csv` | Finite-difference transition-rate elasticity matrix. |
| `DiabetesCtmcToOdeV31_derived_elasticities_kkt.csv` | KKT-based elasticity check. |
| `DiabetesCtmcToOdeV31_derived_numbers_for_tex.txt` | Detailed manuscript-related intermediate and final quantities. |
| `DiabetesCtmcToOdeV31_derived_numbers_4sf.txt` | Compiled manuscript numerical quantities rounded to four significant figures. |
| `DiabetesCtmcToOdeV31_derived_review_checks.txt` | Additional numerical and theoretical-consistency checks. |

### Figures (`Figures/`)

| Generated PNG | What it shows |
| --- | --- |
| `DiabetesCtmcToOdeV31_ode_ssa_all_prediabetes.png` | Deterministic versus Gillespie trajectories and aggregated complication burden, starting from prediabetes. |
| `DiabetesCtmcToOdeV31_ode_ssa_dd2_baseline.png` | The corresponding plots for the second initial cohort. |
| `DiabetesCtmcToOdeV31_literature_vs_calibrated.png` | Disease-progression and survival comparisons across reference, exact-completion, and MaxCal generators. |
| `DiabetesCtmcToOdeV31_elasticities.png` | Heat map of calibrated-rate elasticities with respect to target proportions. |
| `DiabetesCtmcToOdeV31_nhanes_state_ratios.png` | Digitised NHANES cycle-specific state composition among non-healthy participants. |

## 9. Troubleshooting

- **`cannot open file 'Data/...'`:** The required CSV is missing or R is running from the wrong working directory. Verify the exact filename and launch from the repository root.
- **`object ... not found`, `undefined columns selected`, or dimension errors:** Check the CSV header names and the expected state/edge entries in Section 4. In particular, retain the required target-state ordering.
- **PNG or Cairo device error:** Check `capabilities("cairo")` and the graphics support of your R installation. On a system without Cairo, the `png_dev()` helper may need to use a supported PNG device instead of `type = "cairo"`.
- **Linear-system or numerical-optimisation error:** The input rates and target proportions may be incompatible with the model's constraints, or the Newton solver may fail for the chosen inputs. Inspect whether inputs are valid, finite, and on the expected scale. The script does not automatically repair arbitrary input datasets.
- **Different simulated results after modifying the code:** Check that both random seeds and the cohort size are unchanged. Small discrepancies may also reflect R/graphics environment differences.

## 10. Reproducibility and scope

- The script uses **base R**, fixed seeds for stochastic simulation, and explicit output filenames.
- The script header identifies version **V31** as updates to figure presentation; it states that its core computations are those of **V22**.
- Literature transition information and NHANES-derived targets are **inputs**, not raw data automatically downloaded or independently estimated by this code.
- The approach assumes the specified states, allowed transitions, time-homogeneous intensities, and chosen stationary target. Interpret outputs within those modelling assumptions.
- Review the provenance, licensing, and permitted redistribution of the input datasets before making them publicly available.

## Citation and licence

For academic use, please cite the **associated research manuscript** when its full bibliographic details are available. The script identifies the companion manuscript file as `DiabetesCtmcToOdeV31.tex`, but it does not supply a published citation or DOI.

**Licence:** No licence is declared in the supplied script. Repository maintainers should add a `LICENSE` file if they intend to specify permissions for reuse and redistribution.
