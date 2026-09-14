# CopyNumBench

This repository benchmarks nearly all computational tools that infer cell-specific copy
numbers (CNs) from single-cell whole-genome sequencing (scWGS) data.  It also contains the
code that generates every figure of the manuscript.

**Contents**

1. [Repository layout](#1-repository-layout)
2. [Installation](#2-installation)
3. [Input data](#3-input-data)
4. [Running the germline benchmark](#4-running-the-germline-benchmark)
5. [Generating the manuscript figures](#5-generating-the-manuscript-figures)
6. [Reproducing everything with one command](#6-reproducing-everything-with-one-command)
7. [Statistical tests](#7-statistical-tests)
8. [Ploidy-inference benchmark (optional)](#8-ploidy-inference-benchmark-optional)
9. [Data model](#9-data-model)
10. [Methods](#10-methods)
11. [Results](#11-results)
12. [License and other notes](#12-license-and-other-notes)

---

## 1. Repository layout

| path | content |
|---|---|
| `main.py` | generates the Snakemake workflow for the germline benchmark and for the real-tumor mode |
| `common.py`, `data2from1.py`, `data3from2.py`, `data4from2and3.py`, `data_tumor.py` | pipeline steps (alignment, haplotype splitting, CN simulation, CNV calling) |
| `ploidy_eval.py`, `ploidy_tools.py` | ploidy benchmarking (experimental ploidy files, ploidy-inference tools) |
| `cnv_gather_results.py`, `bench_results/` | result tables, figures and statistical tests |
| `cnv_clustermap.py` | clustered CN heatmaps (main Fig. 4 and SI Figs. S23-S30) |
| `cnv_heatmap_montage.py` | merges the per-caller heatmaps into the multi-page SI figure |
| `entire_pipeline.sh` | runs all pipeline and figure steps, then copies the display items into the manuscript directory |
| `scDNAaccessions.tsv`, `scDNAaccessions.S04.tsv` | sample tables for the full and the fast test run |

The benchmark writes its intermediate and final files next to the repository:

| directory | content |
|---|---|
| `../data/` | germline benchmark: simulated data and CN-calling results (`data1` ... `data4`) |
| `../real_tumor_data/` | real-tumor runs (GIAB HG008, ACT breast cancers, ...) |
| `../manuscript_data/`, `../manuscript_data/raw_figs/` | submission package and the figure files the SI LaTeX and the Word files use |

---

## 2. Installation

```bash
bash -evx install_step1_by_conda.sh
bash -evx install_step2_by_download_and_setup.sh
pushd data1to2code && bash -evx install_soft1to2.sh && popd
pushd data3to4code && bash -evx install_soft3to4.sh && popd
```

The installation scripts (including the database download) run for about one day in total
in China.  They create the `copy-num-bench-scwgs` conda environment used by every command
below; activate it first:

```bash
conda activate copy-num-bench-scwgs
```

## 3. Input data

Download the FASTQ files described in `scDNAaccessions.tsv` (full run) or
`scDNAaccessions.S04.tsv` (much faster test run) into

```
../data/1from0.datdir/
```

This is the only download step; it is not part of the pipeline or of `entire_pipeline.sh`.
The full run takes about one month, the test run about one day, on a cluster with 200 CPUs.

---

## 4. Running the germline benchmark

```bash
# 1. generate the workflow (add --SraRunTable scDNAaccessions.S04.tsv for the test run)
python main.py > Snakefile

# 2. run it
snakemake --cores ${NUM_CPUS}

# 3. collect the per-cell results ...
BENCHMARK_RESULT_FILE_PREFIX=bench_results/bench-results-26-04-12-updated
python cnv_gather_results.py -i ../data/*/4from3_*.datdir/*.perf.json \
    -o ${BENCHMARK_RESULT_FILE_PREFIX}

# ... and plot them (Fig. 2, SI Figs. S1-S22, statistics and the LaTeX table)
cat ${BENCHMARK_RESULT_FILE_PREFIX}.long.tsv | \
    python bench_results/scWGS-performances-eval.py -t 0 \
        -o ${BENCHMARK_RESULT_FILE_PREFIX}.plots
```

`-t 0` writes all figures (`<prefix>.plots-all.pdf`, 22 pages, = SI Figs. S1-S22) plus the
main-text figures (`<prefix>.plots_final_grid.*`, `<prefix>.plots_main.*`,
`<prefix>.plots_multirow_main.*`).  `-t 2` redraws only the main figures.

Each metric is evaluated on the cell type it applies to, and every metric label states the
scope: `intCN_accuracy` and `CN_genome_cov_frac` use all simulated cells; the PCCs
(`intCN_PCC`, `nonintCN_PCC`) and the breakpoint metrics use the aneuploid cells only (the
normal/diploid simulations have a constant CN = 2 ground truth and therefore no expected CN
transitions); and `intCN_modal_frac`, the base-pair-weighted fraction of the CNV-call-covered
genome assigned to the modal observed CN state, uses the diploid (normal) cells only.
`intCN_modal_frac` is computed from the observed calls alone, so it is identical under Hap_0
and Hap_1 and is shown once (legend key `N/A`).

The columns (callers) are ordered by each tool's exact publication date and labelled with the
publication year, e.g. `Ginkgo\n2015`; the metric rows are grouped and sorted by importance
(`intCN_accuracy`, `intCN_PCC`, `nonintCN_PCC`, `breakpoint_f1score`, `breakpoint_precision`,
`breakpoint_recall`, `intCN_modal_frac`, `CN_genome_cov_frac`) — the modal-CN fraction and the
covered-genome fraction are the second-last and last rows.  The same year-suffixed method names
and chronological method order are used by the ploidy figures
(`bench_results/scWGS-ploidy-performances-eval.py`), the CNV heatmap titles
(`cnv_clustermap.py`) and the scRNA-seq companion repository
(`plot_cnv_heatmaps.py`).

---

## 5. Generating the manuscript figures

This is the complete, minimal command sequence that produces every display item of the
manuscript, starting from the (already downloaded) FASTQ files.  FASTQ downloads are
deliberately not included -- see [Input data](#3-input-data).

| display item | file(s) | produced by |
|---|---|---|
| Fig. 1 | not code-generated (hand-drawn scheme; panels c/d use FACS/DAPI and karyotype images) | - |
| Fig. 2 | `${BENCHMARK_RESULT_FILE_PREFIX}.plots_final_grid.png` / `.pdf` | step 1 |
| SI Figs. S1-S22 | `${BENCHMARK_RESULT_FILE_PREFIX}.plots-all.pdf` (22 pages) | step 1 |
| Fig. 3 | `${PLOIDY_PREFIX}_main_four.png` / `.pdf` (plus `_main_<group>` panels) | steps 1, 2a, 3 |
| Fig. 4 | `4from2_2_SRP047086_GIAB-HG008_..._4_step2_ginkgo_clustermap.png` / `.pdf` | step 2b |
| SI Figs. S23-S30 | `bench_results/cnv-heatmap.pdf` (8 caller pages) | steps 2b, 3 |
| Fig. 5 | `heatmaps/swarm_grid_metric_by_method.pdf` / `.pdf.png` | step 4 |
| SI Figs. S31-S36 | `heatmaps/swarmHeatmap_*_tumor_only.pdf`, `heatmaps/swarmHeatmap_Fraction_*.pdf`, `heatmaps/swarmHeatmap_Tumor_normal_classification_*.pdf` | step 4 |
| SI Table S1 | `heatmaps/dataset_summary.tsv` | step 4 |

Step 5 copies these files into `raw_figs/` under the exact names the SI LaTeX expects, and
additionally as `Fig2_...`-`Fig5_...` PNG/PDF pairs that can be pasted into the Word files
(the main-text figures are pictures inside the `.docx` files).  It also copies the pairwise
booktabs LaTeX tables written by `bench_results/scWGS-performances-eval.py` and
`bench_results/scWGS-ploidy-performances-eval.py`
(`${BENCHMARK_RESULT_FILE_PREFIX}.plots.stats.pairwise.tex` and
`${PLOIDY_PREFIX}.pooled.stats.pairwise.tex`), which `cnb-g-9-supp-FigsAndTables-k.tex` uses
for Supplementary Tables S2 and S3, plus the scRNA-seq caller pairwise table
`${SCRNA}/heatmaps/stats.pairwise.tex` whenever the companion scRNA-seq repository
provides it: that file is owned by that repository (which does not produce it yet, so the
copy is normally a no-op that only warns - see the note on a missing table below).  Each
comparison in the two scWGS tables reports, in this order and nothing else: the effective
sample size $n$ (donors), the multiplicity-adjusted $p$, the effect size $r$ and the 95% CI
of $r$; the scRNA-seq table is expected to follow the same four-statistic schema.

```
# Run from the repository root; the paths below follow the layout of the sections above:
# this repository next to the results directories (../data, ../real_tumor_data) and to
# ../manuscript_data, the submission package with the SI .tex and the Word files.
# SCRNA points at the companion scRNA-seq repository (copy-num-bench-scwgs-to-scrna);
# set it first if that repository is not cloned next to this one.
REPO=${REPO:-$PWD}
DATA=${DATA:-${REPO}/../data}
REAL=${REAL:-${REPO}/../real_tumor_data}
SCRNA=${SCRNA:-${REPO}/../copy-num-bench-scwgs-to-scrna}
MANUSCRIPT=${MANUSCRIPT:-${REPO}/../manuscript_data}
DEST=${DEST:-${MANUSCRIPT}/raw_figs}
BENCHMARK_RESULT_FILE_PREFIX=bench_results/bench-results-26-07-15-updated
PLOIDY_PREFIX=${BENCHMARK_RESULT_FILE_PREFIX}.ploidy-plots
HG008=${REAL}/SRP047086_GIAB-HG008
HG008_STEM=4from2_2_SRP047086_GIAB-HG008_3_cell-culture_ILLUMINA_nan_4_step2

# ---- step 1: germline (near-haploid) benchmark ------------------------------------
# Produces the per-cell performance files (Fig. 2, SI S1-S22) AND the ploidy-evaluation
# summaries of the simulated cell lines (left panels of Fig. 3).
python main.py --ploidy-tools scabsolute aneufinder chisel copynumber flcna ginkgo \
    hmmcopy sccnv scyn secnv > Snakefile
snakemake --cores 200
python cnv_gather_results.py -i "${DATA}"/*/4from3_*.datdir/*.perf.json \
    -o ${BENCHMARK_RESULT_FILE_PREFIX}
cat ${BENCHMARK_RESULT_FILE_PREFIX}.long.tsv | \
    python bench_results/scWGS-performances-eval.py -t 0 -o ${BENCHMARK_RESULT_FILE_PREFIX}.plots

# ---- step 2: real-tumor arms -------------------------------------------------------
# (a) ACT breast cancers: right panel of Fig. 3.  CHISEL is excluded because it needs
#     phased genotypes, which the ACT samples lack; the Fig. 3 ACT panel shows it as
#     not applicable.
python main.py --tumor-fastq --SraRunTable real_tumor/PRJNA629885_SraRunTable.csv_ref.tsv \
    --ploidy-file ploidy.PRJNA629885.tsv --ploidy-tools scabsolute --excluded-tools chisel \
    > real_tumor/PRJNA629885_SraRunTable_v02_with_scabsolute.snake
snakemake -s real_tumor/PRJNA629885_SraRunTable_v02_with_scabsolute.snake \
    --cores 200 --rerun-incomplete

# (b) GIAB HG008: Fig. 4 and SI S23-S30 (one heatmap per caller, written next to the BEDs)
python main.py --tumor-fastq \
    --SraRunTable real_tumor_metadata/GIAB_PRJNA200694_SraRunTable_bioskryb_ILM.tsv \
    > real_tumor/GIAB_PRJNA200694_SraRunTable_bioskryb_ILM.snake
snakemake -s real_tumor/GIAB_PRJNA200694_SraRunTable_bioskryb_ILM.snake \
    --cores 200 --rerun-incomplete

# ---- step 3: figures built from the results of steps 1-2 ---------------------------
# Fig. 3: both arms in one balloon-plot figure (one -i per arm)
python bench_results/scWGS-ploidy-performances-eval.py \
    -i "${DATA}"/*/*ploidy*eval*summary.json \
       "${REAL}"/SRP259526_*/*ploidy*eval*summary.json \
    -o ${PLOIDY_PREFIX}

# SI S23-S30: the eight single-page HG008 heatmaps merged into one multi-page PDF
python cnv_heatmap_montage.py -i "${HG008}/${HG008_STEM}"*_clustermap.png \
    -o bench_results/cnv-heatmap.pdf

# ---- step 4: scRNA-seq co-sequencing benchmark (Fig. 5, SI S31-S36, Table S1) ------
cd "${SCRNA}"
for yaml in coseq_configs/config_*.yaml; do
    snakemake --configfile config_template.yaml "${yaml}" \
        -s snakemake_pipeline/new_workflow.snake --cores 80 --rerun-incomplete
done
python plot_cnv_heatmaps.py               # -> heatmaps/swarmHeatmap_*.pdf,
                                          #    heatmaps/swarm_grid_metric_by_method.pdf[.png],
                                          #    heatmaps/dataset_summary.tsv

# ---- step 5: copy the display items into the manuscript directory ------------------
cd "${REPO}"
mkdir -p "${DEST}"
cp ${BENCHMARK_RESULT_FILE_PREFIX}.plots-all.pdf "${DEST}/"
cp bench_results/cnv-heatmap.pdf "${DEST}/"
cp "${SCRNA}/heatmaps/dataset_summary.tsv" "${DEST}/"
for f in swarmHeatmap_Tumor_normal_classification_ROC_AUC_scRNA_aneuploidy_score_vs_scWGS_aneuploidy_status.pdf \
         swarmHeatmap_Pearson_Correlation_Coefficient_tumor_only.pdf \
         swarmHeatmap_CopyNumber_gain_ROC_AUC_tumor_only.pdf \
         swarmHeatmap_CopyNumber_loss_ROC_AUC_tumor_only.pdf \
         swarmHeatmap_Fraction_of_the_exome_with_inferred_copy_numbers.pdf \
         swarmHeatmap_Fraction_of_the_cells_with_inferred_copy_numbers.pdf ; do
    cp "${SCRNA}/heatmaps/${f}" "${DEST}/${f}"
done
# images of the main figures for the Word files (PNG) and for vector editing (PDF)
cp ${BENCHMARK_RESULT_FILE_PREFIX}.plots_final_grid.png "${DEST}/Fig2_scWGS_caller_grid.png"
cp ${BENCHMARK_RESULT_FILE_PREFIX}.plots_final_grid.pdf "${DEST}/Fig2_scWGS_caller_grid.pdf"
cp ${PLOIDY_PREFIX}_main_four.png                      "${DEST}/Fig3_ploidy_balloons.png"
cp ${PLOIDY_PREFIX}_main_four.pdf                      "${DEST}/Fig3_ploidy_balloons.pdf"
cp "${HG008}/${HG008_STEM}"*_ginkgo_clustermap.png "${DEST}/Fig4_HG008_ginkgo_heatmap.png"
cp "${HG008}/${HG008_STEM}"*_ginkgo_clustermap.pdf "${DEST}/Fig4_HG008_ginkgo_heatmap.pdf"
cp "${SCRNA}/heatmaps/swarm_grid_metric_by_method.pdf.png" "${DEST}/Fig5_scRNA_swarm_grid.png"
cp "${SCRNA}/heatmaps/swarm_grid_metric_by_method.pdf"     "${DEST}/Fig5_scRNA_swarm_grid.pdf"
# pairwise booktabs LaTeX tables used by the SI (Supplementary Tables S2-S3)
cp "${BENCHMARK_RESULT_FILE_PREFIX}.plots.stats.pairwise.tex" "${DEST}/"
cp "${PLOIDY_PREFIX}.pooled.stats.pairwise.tex" "${DEST}/"
# scRNA-seq caller pairwise table (n, p, r, 95% CI of r), when the companion
# repository provides it; the file is absent while that repository does not write it
if [ -f "${SCRNA}/heatmaps/stats.pairwise.tex" ]; then
    cp "${SCRNA}/heatmaps/stats.pairwise.tex" "${DEST}/"
fi
```

Notes:

* Step 3 only plots; it never re-runs a caller, so it can be repeated after changes to the
  plotting code (`cnv_clustermap.py`, `cnv_heatmap_montage.py`, `bench_results/*.py`).
* If a `cnv_clustermap.py` change should reach an existing HG008 run, re-run the per-caller
  scripts directly, or run `entire_pipeline.sh` (it re-renders them by default,
  `RERUN_HEATMAPS=1`):
  `for f in "${HG008}"/2into4_*_4_step2_*.logdir/*_tumor_clustermap.sh; do bash -evx "$f"; done`
* The main-text figures live as pictures inside the Word files, so Fig. 2-Fig. 5 have to be
  pasted there from the `Fig*` files above; the SI (`cnb-g-3-suppall-k.tex` plus its two
  `\input` files) compiles directly from `raw_figs/`.
* Fig. 1 is drawn by hand; its only code-adjacent panel is the karyotype credit line, which
  is a legend-text edit.

The plotting scripts and their key options:

| script | figures | key options |
|---|---|---|
| `bench_results/scWGS-performances-eval.py` | Fig. 2, SI S1-S22 | `-t 0/1/2` (all / testing / main only), `-o` output prefix, `--stats-*` |
| `bench_results/scWGS-ploidy-performances-eval.py` | Fig. 3 | one `-i` glob per arm, `-o` output prefix, `--stats-reference` |
| `cnv_clustermap.py` | Fig. 4, SI S23-S30 | `--tool/--donor/--sample-type/--avg-spot-len` (title and provenance), `--cluster-annotation K` (row colour bar + cluster mapping TSV), `--max-inline-labels N` (hide row labels above N rows; 0 = always draw), `--bin-size`, `--fai` |
| `cnv_heatmap_montage.py` | SI S23-S30 | `-i` per-caller heatmaps, `-o` multi-page PDF, `--order`, `--rasterize-dpi` |
| `plot_cnv_heatmaps.py` (companion repository) | Fig. 5, SI S31-S36, Table S1 | `--input_glob`, `--outdir`, `--summary_out`, `--metrics`, `--no_tumor_only` |

---

## 6. Reproducing everything with one command

`entire_pipeline.sh` runs the steps of the previous section in order and copies the results
into `${DEST}`.  The default is
`../manuscript_data/raw_figs_<commit id>-<clean|dirty>` (the commit id of this repository,
e.g. `raw_figs_4f7495b-dirty`), so every run lands in a directory that identifies the code
version that produced its display items; a `DEST` given on the command line gets the same
suffix appended (the directory is created if missing):

The copy step also places the pairwise booktabs LaTeX tables used by
`cnb-g-9-supp-FigsAndTables-k.tex` into `${DEST}`:
`${PREFIX}.plots.stats.pairwise.tex` and
`${PLOIDY_PREFIX}.pooled.stats.pairwise.tex`, plus the scRNA-seq caller table
`${SCRNA}/heatmaps/stats.pairwise.tex` when the companion repository has produced it (a
missing table only warns, so the copy step never fails; that repository does not write the
file yet).  A versioned directory therefore contains the display items together with the
LaTeX source of Supplementary Tables S2 and S3.  Each row of the two scWGS tables carries
exactly four statistics, in this order: the effective sample size $n$, the adjusted $p$, the
effect size $r$ and the 95% CI of $r$.

Both code repositories' commit ids, `-clean`/`-dirty` state, commit messages and full
uncommitted diffs are printed once at the very start and once at the very end of the run
(so code edits made while the pipeline was running show up in the second report).  Both
reports are written to `${GIT_SNAPSHOT_TXT}` and, at the very end, copied into `${DEST}`
under the same file name (`entire_pipeline.git-snapshot.txt` by default; skipped together
with the rest of the copy step when `SKIP_COPY=1`).

```bash
# Modify the absolute file-path names in entire_pipeline.sh

bash entire_pipeline.sh                     # pipelines + figures + copy

# figures only, from results that are already on disk:
SKIP_GERMLINE_RUN=1 SKIP_ACT=1 SKIP_HG008_RUN=1 SKIP_SCRNA_RUN=1 bash entire_pipeline.sh

# keep the existing HG008 heatmaps (RERUN_HEATMAPS=1, the default, re-renders them):
RERUN_HEATMAPS=0 bash entire_pipeline.sh

# refresh only Fig. 2 (skip the other figure steps):
SKIP_GERMLINE_RUN=1 SKIP_ACT=1 SKIP_HG008_RUN=1 SKIP_SCRNA_RUN=1 \
    SKIP_FIG3=1 SKIP_FIG4=1 SKIP_FIG5=1 bash entire_pipeline.sh

# recompile the SI after copying:
RECOMPILE_SI=1 bash entire_pipeline.sh
```

Switches (0 = run, 1 = skip):

| switch | effect |
|---|---|
| `SKIP_GERMLINE_RUN`, `SKIP_ACT`, `SKIP_HG008_RUN`, `SKIP_SCRNA_RUN` | skip the heavy pipeline runs (reuse existing results) |
| `SKIP_FIG2`, `SKIP_FIG3`, `SKIP_FIG4`, `SKIP_FIG5` | skip the figure steps (Fig. 2+SI S1-S22, Fig. 3, Fig. 4+SI S23-S30, Fig. 5+SI S31-S36+Table S1) |
| `SKIP_COPY` | skip the copy step (including the git snapshot report) |
| `RERUN_HEATMAPS` | re-run the per-caller HG008 `cnv_clustermap.py` scripts even when their outputs exist (default 1; set 0 to keep them) |
| `RECOMPILE_SI` | recompile the SI LaTeX in `${MANUSCRIPT}` after copying (default 0) |

Further variables: `CORES` (default 200), `SCRNA_CORES` (default 80), `CONDA_ENV` (default
`copy-num-bench-scwgs`; `CONDA_ENV=""` keeps the current environment), `REPO`, `DATA`,
`REAL`, `SCRNA`, `MANUSCRIPT`, `DEST`, `PREFIX`, `PLOIDY_PREFIX`, `STATS_TEX`,
`PLOIDY_STATS_TEX`, `HEATMAP_PDF`, `HG008_DONOR`, `GIT_SNAPSHOT_TXT` (default
`${REPO}/bench_results/entire_pipeline.git-snapshot.txt`).  The code-version suffix
(`<commit id>-<clean|dirty>` of this repository) is always appended to `DEST`.

---

## 7. Statistical tests

The tests run as part of the evaluation commands above: `scWGS-performances-eval.py` runs
them right after loading the long TSV, and `scWGS-ploidy-performances-eval.py` runs them
after the balloon plots.  `bench_results/stat_tests.py` is the shared module and also works
standalone:

```bash
# Fig. 2 (CNV-caller benchmark): statistics on the long TSV
cat ${BENCHMARK_RESULT_FILE_PREFIX}.long.tsv | python bench_results/scWGS-performances-eval.py \
    --stats-only -o ${BENCHMARK_RESULT_FILE_PREFIX}.plots            # stats only, no figures
cat ${BENCHMARK_RESULT_FILE_PREFIX}.long.tsv | python bench_results/stat_tests.py \
    -o ${BENCHMARK_RESULT_FILE_PREFIX}.stats --reference ginkgo     # same thing, standalone
cat ${BENCHMARK_RESULT_FILE_PREFIX}.long.tsv | python bench_results/stat_tests.py \
    -o ${BENCHMARK_RESULT_FILE_PREFIX}.stats.fine --reference ginkgo \
    --cluster-key accession_1,accession_2,cellLine                 # finer sensitivity analysis
python bench_results/test_stat_tests.py                            # self-test + independence demo

# Fig. 3 (ploidy benchmark): statistics on the balloon-plot table
python bench_results/scWGS-ploidy-performances-eval.py -i '*_ploidy_*eval_summary.json' \
    -o ${PLOIDY_PREFIX} --stats-reference 'ginkgo|10'
python bench_results/stat_tests.py -i ${PLOIDY_PREFIX}_pct_within_long.tsv \
    -o ${PLOIDY_PREFIX} --reference 'ginkgo|10'                      # standalone
```

### 7.1 Which tests, and why

Every caller is evaluated on the same simulated cells, so per-cell performances are paired
(blocked) by cell: this pairing handles the correlation **across callers within a cell**.
The metrics are bounded and non-normal, so all tests are nonparametric and two-sided
(two-tailed).

**Independence assumption, and what was done about it.**  A block design further assumes
that the **blocks themselves are mutually independent** -- i.e. that the many per-cell
results produced by the *same* caller are independent draws.  That is not credible for this
benchmark: the ~1,989 simulated cells are all downsamplings of the haplotype-normalized BAMs
of only **nine donors** scored against **three COSMIC templates** (COLO-829, HCC1395, HeLa),
with ITH deletions nested across CNA percentages, so per-cell results of the same caller --
and the paired differences between callers -- are positively correlated within **donors**
(each donor contributes many simulated cells, and donor-specific coverage fluctuations and
artifacts are preserved in every one of them).  Treating the ~1,989 cells as independent
observations is **pseudoreplication**: with intra-cluster correlation ICC and average
cluster size m, the variance is inflated by the design effect `DE = 1 + (m-1)*ICC`; for nine
donors with ~221 cells per donor and a modest ICC of 0.3, DE ~ 67, so P values come out
orders of magnitude too small, Holm no longer controls the family-wise error, Kendall's W is
inflated, and i.i.d. cell-level bootstrap CIs lose coverage.  `python
bench_results/test_stat_tests.py` ends with a measured demonstration
(`demo_independence_failure`): under a true null with clustered data, the naive per-cell
Wilcoxon rejects far more often than the nominal 5%, while the donor/cluster-level test
stays at the nominal level.

Therefore, **Fig. 2 inference runs at the donor level**: each of the nine human donors
contributes ONE effective observation per caller and ground-truth scenario (the median of
that donor's per-cell evaluations), while per-cell quantities are kept as descriptive
statistics and flagged naive comparisons:

* **Fig. 2 default cluster key: `donor`** -- all simulated cells of one donor share that
  donor's haploid material, so every donor is one effective sample and per-cell results are
  aggregated to per-donor medians before testing.  `--cluster-key none` reverts to the naive
  per-cell tests (discouraged).  Finer keys (e.g. `donor,cellLine` or
  `accession_1,accession_2,cellLine`) are available as sensitivity analyses, but they split
  one donor's cells into several clusters and should not be treated as independent-donor
  samples.
* **Fig. 3**: within each plot group the first usable column of
  `donor -> cellLine -> dataset` defines the clusters -- germline-derived datasets of one
  donor share that donor's material; real-tumor samples (no donor metadata) are their own
  units.  Datasets with missing donor labels collapse into one shared `(missing)` cluster
  (the conservative choice).

### 7.2 Tests and outputs (all two-sided)

* **Omnibus, per (ground-truth scenario, metric):** Friedman test across the k callers on
  **per-donor caller medians** (rows `level = cluster` in `*.stats.friedman.tsv`; the cluster
  is the donor in Fig. 2), with the per-cell Friedman kept as a flagged naive row
  (`level = cell (naive)`); Kendall's W and per-caller mean ranks at both levels.
* **Post-hoc pairwise:** two-sided Wilcoxon signed-rank tests on the **per-donor medians of
  the paired differences** (reference caller vs. every other; `--stats-all-pairs` for all 36
  pairs), plus the exact two-sided sign test (`pvalue_sign_test`) as an anchor for very small
  donor counts; Holm-Bonferroni correction within each (scenario, metric) family applied to
  the donor-level P values.  For the ploidy benchmark the same tests pair the per-sample
  percentage of cells within the window by dataset, within each plot group (COLO-829 /
  HCC1395 / HeLa / ACT), aggregated to the chosen cluster level.
* **Effect sizes:** donor/cluster-level matched-pairs rank-biserial correlation r (positive =
  reference better) with a **95% percentile-bootstrap CI** obtained by resampling the
  independent units (`ci95_r_low`/`ci95_r_high`), the paired common-language effect size
  P(ref > other) + 0.5 P(tie) on the cell population, and the median per-cell difference
  with a **cluster bootstrap 95% CI**
  (donors/clusters resampled with replacement, all cells of a drawn donor kept together;
  10,000 resamples requested by default, capped at 5,000 for runtime; seeded, so all
  intervals are exactly reproducible).
* **Diagnostics of the independence assumption, per comparison** (`*.stats.pairwise.tsv`):
  `icc_within_cluster_d` (ICC(1,1) of the paired differences within donor groups, one-way
  ANOVA method of moments), `design_effect` = 1 + (m0-1)*max(ICC, 0), `n_effective_cells` =
  n/DE, and `p_inflation_ratio` = donor-level P / naive per-cell P (1 = no dependence
  detected; large = the naive P overstates significance).  The naive per-cell P value itself
  is kept in `pvalue_cell_naive` for comparison.
* **Hap_0 vs. Hap_1 agreement (the CNP/aneuploidy fix):** Spearman's rho between the
  per-caller median-performance rankings under the two ground-truth scenarios
  (`*.stats.concordance.tsv`), plus two-sided Wilcoxon signed-rank tests of the within-caller
  scenario difference at the donor level (rows with `caller_b = '<scenario>'` in
  `*.stats.pairwise.tsv`, Holm-corrected within each metric family).
* `stat_tests.mcnemar_exact()` is provided for per-cell paired binary outcomes (e.g.
  within/outside the ploidy window for two tools evaluated on the same cells) -- aggregate
  per donor/cluster or restrict to one dataset before using it, since it too assumes
  independent cells.

### 7.3 Reported values and options

Exact two-sided P values, effect sizes and confidence intervals are reported in full in
`*.stats.pairwise.tsv` (columns `pvalue_two_sided`, `pvalue_sign_test`, `pvalue_holm`,
`rank_biserial_r`, `ci95_r_low/high`, `cl_effect_paired`, `ci95_median_diff_low/high`,
`icc_within_cluster_d`,
`design_effect`, `n_effective_cells`, `p_inflation_ratio`); `n_pairs` is the number of units
the primary test runs on (donors by default in Fig. 2; `n_cells_paired`/`n_clusters` give the
cell and donor/cluster counts).  Asymptotic P values that underflow double precision
(P < 2.2e-16) are noted in the `note` column, as are donor counts too small for the exact
Wilcoxon to reach P < 0.05 (< 6 donors).  The exact n for every analysis, the chosen cluster
key, the donor count and size distribution, and the full independence record are written to
`*.stats.json`, together with the test settings, the bootstrap seed and the software
versions.  After every Fig. 2 stats run the script also writes `<output>.stats.pairwise.tex`,
a booktabs LaTeX table with one row per (metric, caller $b$); each of the two ground-truth
scenario groups (**Hap_0**, **Hap_1**) owns four sub-columns that carry, in this order and
nothing else, the effective sample size $n$ (donors), the Holm-adjusted $p$, the effect size
$r$ and the 95% bootstrap CI of $r$ -- so the manuscript table always shows the donor-level
values.

Options: `--no-stats` disables the tests; `--stats-reference`, `--stats-boot`,
`--stats-seed`, `--stats-alpha` tune them; `--stats-pair-key` overrides the columns that
identify one simulated cell (default `accession_1,accession_2,cellLine,overall_ploidy,
CNA_percent`); `--stats-cluster-key` / `--cluster-key` overrides the cluster columns
(`none` = naive per-cell level, discouraged); `--stats-no-cluster` is an alias of `none`;
`--cluster-agg` selects median (default) or mean aggregation within donors/clusters.
Nothing else in the pipeline changes: no run rule, no Snakefile difference, same figures.

---

## 8. Ploidy-inference benchmark (optional)

Tools such as [scAbsolute](https://doi.org/10.1186/s13059-024-03204-y) report a ploidy
estimate instead of a per-cell copy-number profile, so they are benchmarked by an opt-in
module of their own (`ploidy_tools.py`) rather than by the CNV-caller pipeline above.
`--ploidy-tools` also accepts any of the CNV callers, so that ploidy inference is benchmarked
across the whole tool set and not only for the tools that report a ploidy: a caller states no
ploidy, so its ploidy is taken to be the length-weighted mean of the integer copy numbers it
called (scAbsolute Eq. 1), read from the `*intcns.bed` files it has already written.
That costs no extra run and leaves the CNV benchmark untouched -- a caller named here is
registered nowhere and gains one evaluation script, no run rule -- but it does mean every
tool is scored on the same quantity, by the same code, into the same tables.
Nothing changes unless `--ploidy-tools` is passed: without it, `main.py` generates exactly
the same Snakefile as before, rule for rule.

```bash
pushd data3to4code && bash -evx install_scabsolute.sh && popd
python main.py --tumor-fastq --SraRunTable ${TUMOR_RUN_TABLE} \
    --ploidy-file ploidy.PRJNA629885.tsv --ploidy-tools scabsolute [--ploidy-facs] > Snakefile
snakemake --cores ${NUM_CPUS}
```

To benchmark ploidy inference for every tool at once, name the callers too (each also has to
be in `--tools`, which it is by default):

```bash
python main.py --tumor-fastq --SraRunTable ${TUMOR_RUN_TABLE} --ploidy-file ploidy.PRJNA629885.tsv \
    --ploidy-tools scabsolute aneufinder chisel copynumber flcna ginkgo hmmcopy sccnv scyn secnv > Snakefile
```

Each ploidy tool writes its per-cell estimates to a ploidy-calls TSV (columns `cell` and
`ploidy`, plus any tool-specific extras) and is then scored against the experimental
(FACS/DAPI) ploidy of `--ploidy-file` with the same scAbsolute metrics that `ploidy_eval.py`
applies to the CNV callers: the percentage of cells outside the `--ploidy-window` around the
experimental estimate, the mean absolute ploidy distance, and the 2x/0.5x scaling-error
diagnostics.
The resulting `*_ploidy_tool_eval_percell.tsv`, `_persample.tsv` and `_summary.json` share the
columns of the caller-side `*_ploidy_eval_*` files and can simply be concatenated, so tools
that infer ploidy and callers that imply it end up in one ranking.
The per-sample table additionally carries `sample_ploidy` and its error columns: the single
per-sample point estimate, derived from the per-cell values the way scAbsolute's own
`scripts/estimatePloidy.R` derives it.
That per-sample estimate is what a caller gains by being named in `--ploidy-tools`; the
`*_ploidy_eval_*` files that `--ploidy-file` produces for it anyway are per-cell only.
`--ploidy-facs` stays scAbsolute-only, since re-running Ginkgo with a ploidy that Ginkgo
itself implied would be circular.

`--ploidy-facs` also measures what the estimate is worth downstream: the calls are converted
into a Ginkgo FACS file and Ginkgo is re-run with them as `ginkgo_facs_<tool>`, which is
evaluated like any other caller and is therefore directly comparable with the untouched
`ginkgo` run.

Both steps also run standalone, on any ploidy-calls TSV -- or, for a caller, on the
`*intcns.bed` files it wrote:

```bash
python ploidy_tools.py eval -i '*_ploidy_calls.tsv' -o ${OUT_PREFIX} \
    --ploidy-file ploidy.PRJNA629885.tsv --metadata-tsv ${TUMOR_RUN_TABLE} --plot
python ploidy_tools.py eval -i '*_intcns.bed' -o ${OUT_PREFIX} --tool ginkgo --chroms autosomes \
    --ploidy-file ploidy.PRJNA629885.tsv --metadata-tsv ${TUMOR_RUN_TABLE} --plot
python ploidy_tools.py facs -i '*_ploidy_calls.tsv' -o ${FACS_FILE}
```

---

## 9. Data model

For a set of cells derived from the same donor (e.g. a human subject):

| data | content |
|---|---|
| `data1` | SRA raw FASTQ files |
| `data2` | reference-aligned BAM files with two germline haplotypes (PB1, PB2 and FPN are separated from each other for oocytes) |
| `data3` | reference-aligned BAM files with copy numbers derived from data downsampling/reverse-downsampling in pre-defined genome intervals |
| `data4` | single-cell copy-number (CN) calling-result TSV and CN benchmark-result TSV files |

`common.py` describes the files generated by `main.py` in more detail.

---

## 10. Methods

In brief, our benchmarking strategy satisfies two seemingly conflicting requirements.

1. Use real sequencing data (not in silico simulated data, such as data simulated from the
   hg19/GRCh37 human reference genome).  The data generated by real wet-lab sequencing
   assays capture the biological and technical variations that are intrinsic to the
   data-generation process.
2. Use data with known ground-truth copy-number profiles for each cell.

## 11. Results

<img width="900" height="300" alt="image" src="https://github.com/user-attachments/assets/2a429b0f-032f-4192-bf46-8bac16b93229" />

Our results show that Ginkgo performs the best for calling copy-number variations/variants
(CNVs) from single-cell whole-genome sequencing data.
The metrics `intCN_accuracy` and `intCN_PCC` (i.e., Pearson correlation coefficient of intCN)
measure the observed (i.e., called) versus expected (i.e., ground-truth) integer copy numbers
(CNs).
The metric `nonintCN_PCC` is the Pearson correlation coefficient of observed non-integer CN
versus expected integer CN.
The metric `CN_genome_cov_frac` refers to the fraction of the human reference genome hg19 that is
covered by the observed CN profile.
The metric `breakpoint_f1score` is the F1-score of detecting CN changes (i.e., breakpoints):
a called CN transition (breakpoint) is precise (i.e., true positive) if it is matched
one-to-one, closest pairs first, with a ground-truth CN transition within 200 000 base pairs,
and a ground-truth CN transition is recalled (i.e., true positive) if it is matched that way
with a called CN transition.
Fig. S1 shows performances as a function of each performance-related factor (e.g., ploidy
estimation accuracy and average sequencing depth).

---

## 12. License and other notes

MIT license.

Please be aware that the commercial (i.e., for-profit and non-academic) use of
`cosmic-v97/cell_lines_copy_number.csv` may require licensing from
cancer.sanger.ac.uk/cosmic.

If you would like to request any additional information (e.g., details on the benchmarking
strategy), please send an email to: cndfeifei AT aliyun DOT com
