# Hybrid breeding-programme simulations in wheat and maize

Forward-in-time simulations of a two-pool hybrid breeding programme, run for two crops under
five genetic architectures with [AlphaSimR](https://cran.r-project.org/package=AlphaSimR).

Supporting code for:

> Crossa, J. et al. *Apparent additivity in crop hybrids arises from the projection of
> dominance and epistasis.*

The simulations are not used to estimate projection ratios. They compare long-term breeding
response when the underlying architecture is additive, epistatic or dominant, under otherwise
identical breeding logic. Parameters correspond to Supplementary Table S1; the heritabilities
to Supplementary Table S2.

## What is simulated

Ten simulations: two crops × five architectures. Each is a two-pool hybrid programme with
phenotypic selection — GCA on testcrosses at YT1–YT2, hybrid performance at YT3–YT5 — over a
20-year burn-in plus 20 evaluation years, replicated 5 times. Each year, 10 of the 50 parents
per pool are recycled from YT4, the elite testers are rotated, and the mean and variance of
genetic value are recorded among the YT3 inbreds and among their hybrids. No genomic
prediction and no genotype-by-year interaction.

The "pure line" and "hybrid" curves of Figure 2 are the `meanG_inbred` and `meanG_hybrid`
columns of the same run, not two separate simulations.

## The two crops

Wheat and maize use the same scripts. Three settings differ:

| Setting | Wheat | Maize |
|---|---|---|
| `species` in `runMacs()` (`CreateParents.R`) | `"WHEAT"` | `"MAIZE"` |
| `nChr` (`GlobalParameters.R`) | 21 | 10 |
| `initMeanG` (`GlobalParameters.R`) | 4.0 t ha⁻¹ | 10.0 t ha⁻¹ |

Everything else is identical, including 300 QTL and 400 markers **per chromosome**,
`initVarG = 0.5` and `VarE = 1.0`. Because chromosome number differs, genome-wide totals
differ (6,300 QTL in wheat, 3,000 in maize). Founder additive and error variances are the
same, so stage heritabilities match and absolute gains are comparable; founder means differ,
so gains as a percentage of the starting mean are not.

## The five architectures

Each `T*` folder differs only in the trait definition.

| Folder | Label | Call | Parameters |
|---|---|---|---|
| `T1/` | Additive only | `addTraitA` | — |
| `T2/` | Low epistatic | `addTraitAE` | `relAA = 0.5` |
| `T3/` | High epistatic | `addTraitAE` | `relAA = 1.0` |
| `T4/` | Low dominance | `addTraitAD` | `meanDD = 0.5`, `varDD = 0.2` |
| `T5/` | High dominance | `addTraitAD` | `meanDD = 1.0`, `varDD = 0.2` |

Founder additive variance is fixed at 0.5 across all five (`useVarA = TRUE`), so they differ
only in the amount and type of non-additive variance. In T1 the hybrid mean equals the inbred
mean, since a purely additive trait has no midparent deviation. In T4 and T5 directional
dominance over thousands of QTL accumulates substantial heterosis, so the hybrid curves sit
well above the inbred ones; that offset is a property of the `addTraitAD` parameterisation and
is not calibrated to observed heterosis, so the architectures are for relative comparison of
response.

## Requirements

R ≥ 4.1, `AlphaSimR` ≥ 1.4, and `dplyr` for the plotting script.

```r
install.packages(c("AlphaSimR", "dplyr"))
```

No input data: founder genomes are simulated with MaCS inside `CreateParents.R`.

## Layout

```
wheat/   T1/ ... T5/        maize/   T1/ ... T5/
```

Each architecture folder is self-contained — all scripts, and the results of running them.

## Running

```bash
cd wheat/T1 && Rscript run.R
```

`run.R` sources everything else in order, so this one call is the whole workflow. Repeat for
each of the ten folders. Two notes:

- **Line 1 of `run.R` is a `setwd()`** to the absolute path of the machine the simulations were
  run on. Delete it or point it at the folder on your own machine; every other path is
  relative.
- **No random seed is set.** Repeated runs are not numerically identical; published results are
  averages over 5 replicates. To reproduce a specific run, add `set.seed(1000 + REP)` inside
  the `for(REP in 1:nReps)` loop. The `.rds` files in this repository are the exact outputs
  behind the published figures and tables.

Roughly 1–3 hours and under 8 GB per folder, dominated by the MaCS founder simulation. The ten
runs are independent and can be run in parallel.

## Workflow, in execution order

| Script | When | What it does |
|---|---|---|
| `GlobalParameters.R` | once | All programme parameters. **Edit this file** to change the programme. |
| `CreateParents.R` | per replicate | 100 founder genomes via MaCS, split into two pools of 50; defines the trait and error variance; picks initial elite parents and testers. |
| `FillPipeline.R` | per replicate | Fills every stage with distinct cohorts, so year 1 starts from a steady-state pipeline. |
| `UpdateParents.R` | each year | Replaces the 10 oldest parents per pool with the best 10 YT4 inbreds. |
| `UpdateTesters.R` | each year | Replaces the oldest elite parent with the parent of the best YT5 hybrid; refreshes YT1/YT2 testers. |
| `AdvanceYear.R` | each year | Advances every stage by one year, backwards from YT5 to YT1 so data are not overwritten. |
| `StoreTrainPop.R` | burn-in years | Accumulates a training population. Retained from the pipeline this code was adapted from; no genomic prediction is run here, so it affects no reported result. |
| `ANALYZERESULTS.R` | once | Writes `Results.png`: inbred and hybrid gain and variance by year. |

`UpdateResults.R` duplicates summaries that `run.R` computes inline; it is not sourced.

Stage sizes and replication, per pool: YT1 4,000 inbreds × 1 tester, r = 1; YT2 400 × 3,
r = 2; YT3 40 × 5 elite testers, r = 4; YT4 20, r = 8; YT5 4, r = 100. The error variance of
an entry mean is σ²/(r × t) at YT1–YT2 and σ²/r at YT3–YT5, with σ² = `VarE` = 1.0. The
`#h2 = ...` comments in `GlobalParameters.R` are inherited from the programme this was adapted
from and refer to its variance scale, not this one.

## Outputs

`HybridPheno.rds` — a list of 5 data frames, one per replicate, 40 rows each.
`do.call(rbind, readRDS("HybridPheno.rds"))` gives one 200-row frame with: `year` (1–20
burn-in, 21–40 evaluation), `rep`, `scenario`, `meanG_inbred`/`varG_inbred` (the 40 YT3
inbreds, averaged over pools), `meanG_hybrid`/`varG_hybrid` (their hybrids), and two
diagnostics not reported in the manuscript, `acc_sel` and `cor`.

`Results.png` — the four diagnostic panels.

Figure 2 of the manuscript is `meanG_inbred` and `meanG_hybrid` by year, averaged over the five
replicates. Supplementary Table S2 is the mean of `varG_inbred` and `varG_hybrid` over the 20
evaluation years, converted to heritability as H² = V_G/(V_G + σ²/r) with σ² = 1.0 and r = 4
at YT3.

## Attribution

Adapted from the phenotypic hybrid maize programme of Jon Bančič, Philip Greenspoon, Chris
Gaynor and Gregor Gorjanc, <https://github.com/HighlanderLab/jbancic_alphasimr_plants>;
original scripts by Chris Gaynor. Modifications set the crop-specific genome and demography,
define the five architectures, and record inbred and hybrid gain side by side.

Please cite the manuscript above and:

> Gaynor, R. C., Gorjanc, G. & Hickey, J. M. (2021). AlphaSimR: an R package for breeding
> program simulations. *G3* **11**, jkaa017.

Questions about the code: Paolo Vitale, CIMMYT.
