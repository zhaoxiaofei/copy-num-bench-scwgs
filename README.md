This code repository evaluates nearly all computational tools to infer cell-specific copy numbers (CNs) from single-cell whole-genome sequencing data (scWGS).

### How to setup

```
bash -evx install_step1_by_conda.sh
bash -evx install_step2_by_download_and_setup.sh
pushd data1to2code && bash -evx install_soft1to2.sh && popd
pushd data3to4code && bash -evx install_soft3to4.sh && popd
```

The above installation scripts (which includes database download) run for about one day in total in China.

Then, download the FASTQ files described in scDNAaccessions.tsv (for performing a full run) or scDNAaccessions.S04.tsv (for performing a test run that is much faster than the full run) into the directory ../data/1from0.datdir/

The full/test run takes about one month/day to finish running on a cluster with 200 CPUs.

### How to benchmark the tools

```
python main.py > Snakefile # or python main.py --SraRunTable scDNAaccessions.S04.tsv > Snakefile # to perform the test run
snakemake --cores ${NUM_CPUS}
# Wait for the above snakemake command to finish
# Set BENCHMARK_RESULT_FILE_PREFIX to be the prefix of the files storing performance-evaluation results
#   for example, BENCHMARK_RESULT_FILE_PREFIX=bench_results/bench-results-26-04-12-updated
python cnv_gather_results.py -i ../data/*/4from3_*.datdir/*.perf.json -o ${BENCHMARK_RESULT_FILE_PREFIX}
cat ${BENCHMARK_RESULT_FILE_PREFIX}.long.tsv | python bench_results/scWGS-performances-eval.py -t 0 -o ${BENCHMARK_RESULT_FILE_PREFIX}.plots
```

### How to run the statistical tests

The statistical tests are part of the evaluation commands above: `scWGS-performances-eval.py` runs them by default right after loading the long TSV, and `scWGS-ploidy-performances-eval.py` runs them after the balloon plots. `bench_results/stat_tests.py` is the shared module and also works standalone:

```
# Fig. 2 (CNV-caller benchmark): statistics on the long TSV
cat ${BENCHMARK_RESULT_FILE_PREFIX}.long.tsv | python bench_results/scWGS-performances-eval.py \
    --stats-only -o ${BENCHMARK_RESULT_FILE_PREFIX}.plots            # stats only, no figures
cat ${BENCHMARK_RESULT_FILE_PREFIX}.long.tsv | python bench_results/stat_tests.py \
    -o ${BENCHMARK_RESULT_FILE_PREFIX}.stats --reference ginkgo     # same thing, standalone
cat ${BENCHMARK_RESULT_FILE_PREFIX}.long.tsv | python bench_results/stat_tests.py \
    -o ${BENCHMARK_RESULT_FILE_PREFIX}.stats.donor --reference ginkgo \
    --cluster-key donor,cellLine                                   # donor-level sensitivity analysis
python bench_results/test_stat_tests.py                            # self-test + independence demo

# Fig. 3 (ploidy benchmark): statistics on the balloon-plot table
python bench_results/scWGS-ploidy-performances-eval.py -i '*_ploidy_*eval_summary.json' \
    -o ${PLOIDY_PREFIX} --stats-reference 'ginkgo|10'
python bench_results/stat_tests.py -i ${PLOIDY_PREFIX}_pct_within_long.tsv \
    -o ${PLOIDY_PREFIX} --reference 'ginkgo|10'                      # standalone
```

Which tests, and why. Every caller is evaluated on the same simulated cells, so per-cell performances are paired (blocked) by cell: this pairing handles the correlation **across callers within a cell**. The metrics are bounded and non-normal, so all tests are nonparametric and two-sided (two-tailed):

**Independence assumption, and what was done about it.** A block design further assumes that the **blocks themselves are mutually independent** — i.e. that the many per-cell results produced by the *same* caller are independent draws. That is not credible for this benchmark: the ~1,989 simulated cells are all downsamplings of the haplotype-normalized BAMs of only **nine donors** scored against **three COSMIC templates** (COLO-829, HCC1395, HeLa), with ITH deletions nested across CNA percentages, so per-cell results of the same caller — and the paired differences between callers — are positively correlated within `(accession_1, accession_2, cellLine)` groups (and, coarser, within donors). Treating the cells as ~1,989 independent observations is **pseudoreplication**: with intra-cluster correlation ICC and average cluster size m, the variance is inflated by the design effect `DE = 1 + (m-1)*ICC`, so a per-cell test with ~45 clusters of ~44 cells and ICC = 0.3 runs with DE ~ 14 (effective n ~ 143, not 1,989): P values come out orders of magnitude too small, Holm no longer controls the family-wise error, Kendall's W is inflated, and i.i.d. cell-level bootstrap CIs lose coverage. `python bench_results/test_stat_tests.py` ends with a measured demonstration (`demo_independence_failure`): under a true null with that exact structure, the naive per-cell Wilcoxon rejects in **~58% of runs at the nominal 5%**, while the cluster-level test stays at **~5%**.

Therefore, **inference runs at the level of the independent experimental unit (the cluster)**, while per-cell quantities are kept as descriptive statistics and flagged naive comparisons:

* Fig. 2 default cluster key: `accession_1,accession_2,cellLine` — the shared biological material (the two haplotype BAMs) + truth template. Note that the looser relation "shares a donor" is not transitive and cannot serve as a single key; the coarser `--cluster-key donor,cellLine` is provided as a donor-level sensitivity analysis. `--cluster-key none` reverts to the naive per-cell tests (discouraged).
* Fig. 3: within each plot group the first usable column of `donor -> cellLine -> dataset` defines the clusters — germline-derived datasets of one donor share that donor's material; real-tumor samples (no donor metadata) are their own units. Datasets with missing donor labels collapse into one shared `(missing)` cluster (the conservative choice).

Tests and outputs (all two-sided):

* **Omnibus, per (ground-truth scenario, metric):** Friedman test across the k callers on **per-cluster caller medians** (rows `level = cluster` in `*.stats.friedman.tsv`), with the per-cell Friedman kept as a flagged naive row (`level = cell (naive)`); Kendall's W and per-caller mean ranks at both levels.
* **Post-hoc pairwise:** two-sided Wilcoxon signed-rank tests on the **per-cluster medians of the paired differences** (reference caller vs. every other; `--stats-all-pairs` for all 36 pairs), plus the exact two-sided sign test (`pvalue_sign_test`) as an anchor for very small cluster counts; Holm-Bonferroni correction within each (scenario, metric) family applied to the cluster-level P values. For the ploidy benchmark the same tests pair the per-sample percentage of cells within the window by dataset, within each plot group (COLO-829 / HCC1395 / HeLa / ACT), aggregated to the chosen cluster level.
* **Effect sizes:** cluster-level matched-pairs rank-biserial correlation r (positive = reference better), paired common-language effect size P(ref > other) + 0.5 P(tie) on the cell population, and the median per-cell difference with a **cluster bootstrap 95% CI** (clusters resampled with replacement, all cells of a drawn cluster kept together; 10,000 resamples requested by default, capped at 5,000 for runtime; seeded, so all intervals are exactly reproducible).
* **Diagnostics of the independence assumption, per comparison** (`*.stats.pairwise.tsv`): `icc_within_cluster_d` (ICC(1,1) of the paired differences within clusters, one-way ANOVA method of moments), `design_effect` = 1 + (m0-1)*max(ICC, 0), `n_effective_cells` = n/DE, and `p_inflation_ratio` = cluster-level P / naive per-cell P (1 = no dependence detected; large = the naive P overstates significance). The naive per-cell P value itself is kept in `pvalue_cell_naive` for comparison.
* **Hap_0 vs. Hap_1 agreement (the CNP/aneuploidy fix):** Spearman's rho between the per-caller median-performance rankings under the two ground-truth scenarios (`*.stats.concordance.tsv`), plus two-sided Wilcoxon signed-rank tests of the within-caller scenario difference at the cluster level (rows with `caller_b = '<scenario>'` in `*.stats.pairwise.tsv`, Holm-corrected within each metric family).
* `stat_tests.mcnemar_exact()` is provided for per-cell paired binary outcomes (e.g. within/outside the ploidy window for two tools evaluated on the same cells) — aggregate per cluster or restrict to one dataset before using it, since it too assumes independent cells.

Exact two-sided P values, effect sizes and confidence intervals are reported in full in `*.stats.pairwise.tsv` (columns `pvalue_two_sided`, `pvalue_sign_test`, `pvalue_holm`, `rank_biserial_r`, `cl_effect_paired`, `ci95_median_diff_low/high`, `icc_within_cluster_d`, `design_effect`, `n_effective_cells`, `p_inflation_ratio`); `n_pairs` is the number of units the primary test runs on (clusters by default; `n_cells_paired`/`n_clusters` give the cell and cluster counts). Asymptotic P values that underflow double precision (P < 2.2e-16) are noted in the `note` column, as are cluster counts too small for the exact Wilcoxon to reach P < 0.05 (< 6 clusters). The exact n for every analysis, the chosen cluster key, the cluster count and size distribution, and the full independence record are written to `*.stats.json`, together with the test settings, the bootstrap seed and the software versions.

Options: `--no-stats` disables the tests; `--stats-reference`, `--stats-boot`, `--stats-seed`, `--stats-alpha` tune them; `--stats-pair-key` overrides the columns that identify one simulated cell (default `accession_1,accession_2,cellLine,overall_ploidy,CNA_percent`); `--stats-cluster-key` / `--cluster-key` overrides the cluster columns (`none` = naive per-cell level, discouraged); `--stats-no-cluster` is an alias of `none`; `--cluster-agg` selects median (default) or mean aggregation within clusters. Nothing else in the pipeline changes: no run rule, no Snakefile difference, same figures.


### How to benchmark ploidy-inference tools (optional)

Tools such as [scAbsolute](https://doi.org/10.1186/s13059-024-03204-y) report a ploidy estimate instead of a per-cell copy-number profile, so they are benchmarked by an opt-in module of their own (`ploidy_tools.py`) rather than by the CNV-caller pipeline above.
`--ploidy-tools` also accepts any of the CNV callers, so that ploidy inference is benchmarked across the whole tool set and not only for the tools that report a ploidy: a caller states no ploidy, so its ploidy is taken to be the length-weighted mean of the integer copy numbers it called (scAbsolute Eq. 1), read from the `*intcns.bed` files it has already written.
That costs no extra run and leaves the CNV benchmark untouched -- a caller named here is registered nowhere and gains one evaluation script, no run rule -- but it does mean every tool is scored on the same quantity, by the same code, into the same tables.
Nothing changes unless `--ploidy-tools` is passed: without it, `main.py` generates exactly the same Snakefile as before, rule for rule.

```
pushd data3to4code && bash -evx install_scabsolute.sh && popd
python main.py --tumor-fastq --SraRunTable ${TUMOR_RUN_TABLE} \
    --ploidy-file ploidy.PRJNA629885.tsv --ploidy-tools scabsolute [--ploidy-facs] > Snakefile
snakemake --cores ${NUM_CPUS}
```

To benchmark ploidy inference for every tool at once, name the callers too (each also has to be in `--tools`, which it is by default):

```
python main.py --tumor-fastq --SraRunTable ${TUMOR_RUN_TABLE} --ploidy-file ploidy.PRJNA629885.tsv \
    --ploidy-tools scabsolute aneufinder chisel copynumber flcna ginkgo hmmcopy sccnv scyn secnv > Snakefile
```

Each ploidy tool writes its per-cell estimates to a ploidy-calls TSV (columns `cell` and `ploidy`, plus any tool-specific extras) and is then scored against the experimental (FACS/DAPI) ploidy of `--ploidy-file` with the same scAbsolute metrics that `ploidy_eval.py` applies to the CNV callers: the percentage of cells outside the `--ploidy-window` around the experimental estimate, the mean absolute ploidy distance, and the 2x/0.5x scaling-error diagnostics.
The resulting `*_ploidy_tool_eval_percell.tsv`, `_persample.tsv` and `_summary.json` share the columns of the caller-side `*_ploidy_eval_*` files and can simply be concatenated, so tools that infer ploidy and callers that imply it end up in one ranking.
The per-sample table additionally carries `sample_ploidy` and its error columns: the single per-sample point estimate, derived from the per-cell values the way scAbsolute's own `scripts/estimatePloidy.R` derives it.
That per-sample estimate is what a caller gains by being named in `--ploidy-tools`; the `*_ploidy_eval_*` files that `--ploidy-file` produces for it anyway are per-cell only.
`--ploidy-facs` stays scAbsolute-only, since re-running Ginkgo with a ploidy that Ginkgo itself implied would be circular.

`--ploidy-facs` also measures what the estimate is worth downstream: the calls are converted into a Ginkgo FACS file and Ginkgo is re-run with them as `ginkgo_facs_<tool>`, which is evaluated like any other caller and is therefore directly comparable with the untouched `ginkgo` run.

Both steps also run standalone, on any ploidy-calls TSV -- or, for a caller, on the `*intcns.bed` files it wrote:

```
python ploidy_tools.py eval -i '*_ploidy_calls.tsv' -o ${OUT_PREFIX} \
    --ploidy-file ploidy.PRJNA629885.tsv --metadata-tsv ${TUMOR_RUN_TABLE} --plot
python ploidy_tools.py eval -i '*_intcns.bed' -o ${OUT_PREFIX} --tool ginkgo --chroms autosomes \
    --ploidy-file ploidy.PRJNA629885.tsv --metadata-tsv ${TUMOR_RUN_TABLE} --plot
python ploidy_tools.py facs -i '*_ploidy_calls.tsv' -o ${FACS_FILE}
```

### Detail about the data

for a set of cells derived from the same donor (e.g., human subject):
    data1 are SRA raw FastQ files
    data2 are reference-aligned BAM files with two germline haplotypes (PB1, PB2 and FPN are separated from each other for oocytes)
    data3 are reference-aligned BAM files with copy numbers derived from data downsampling/reverse-downsampling in pre-defined genome intervals
    data4 are single-cell copy number (CN) calling-result TSV and CN benchmark-result TSV files

The file common.py describes the files generated by main.py in more detail

### Methods

In brief, our benchmarking strategy satisfies two seemingly conflicting requirements. 

1 - Use real sequencing data (not in silico simulated data, such as the ones simulated from the hg19/GRCh37 human reference genome). The data generated by real wet-lab sequencing assays captures the biological and technical variations that are intrinsic to the the data-generation processs.

2 - Use data with known ground truth copy number profiles for each cell. 

### Results

<img width="900" height="300" alt="image" src="https://github.com/user-attachments/assets/2a429b0f-032f-4192-bf46-8bac16b93229" />

Our results show that Ginkgo performs the best for calling copy-number variations/variants (CNVs) from single-cell whole-genome sequencing data. 
The metrics `accuracy` and `PCC_intCN` (i.e., Pearson correlation coefficient of intCN) measure the observed (i.e., called) versus expected (i.e., ground-truth) integer copy numbers (CNs). 
The metric `PCC_nonintCN` is the Pearson correlation coefficient of observed non-integer CN versus expected integer CN. 
The metric `frac_cov_genome` refers to the fraction of the human reference genome hg19 that is covered by the observed CN profile. 
The metrics `breakpoint_f1score` is the F1-score of detecting CN changes (i.e., breakpoints): 
An observed breakpoint is precise (i.e., true positive) if at least one expected breakpoint is within 200 000 base pairs of the observed breakpoint, 
and an expected breakpoint is recalled (i.e., true positive) if at least one observed breakpoint is within 200 000 base pairs of the expected breakpoint. 
Fig. S1 shows performances as a function of each performance-related factor (e.g., ploidy estimation accuracy and average sequencing depth). 

### LICENSE

MIT license

### Other things

Please be aware that the commercial (i.e., for profit and non-academic) use of `cosmic-v97/cell_lines_copy_number.csv` may require licensing from cancer.sanger.ac.uk/cosmic

If you would like to request any additional information (e.g., details on bechmarking strategy), please send an email to: cndfeifei AT aliyun DOT com
