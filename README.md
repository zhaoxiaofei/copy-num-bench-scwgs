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
