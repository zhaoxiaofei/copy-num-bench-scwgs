#!/usr/bin/env bash
#
# entire_pipeline.sh - regenerate every code-generated display item of the manuscript
# and copy the results into the manuscript directory (default:
# /stor/zxf/cnv/manuscript_data/raw_figs).
#
# The FASTQ files must already be in ../data/1from0.datdir/ - FASTQ downloads are
# deliberately NOT part of this script (see README.md, "How to setup").
#
# What it produces (see README.md, "How to generate every manuscript figure"):
#   Fig. 2 + SI S1-S22   ${PREFIX}.plots_final_grid.{png,pdf}, ${PREFIX}.plots-all.pdf
#   Fig. 3               ${PLOIDY_PREFIX}_main_four.{png,pdf} (+ per-group panels)
#   Fig. 4 + SI S23-S30  HG008 per-caller heatmaps + ${HEATMAP_PDF}
#   Fig. 5 + SI S31-S36  ${SCRNA}/heatmaps/swarmHeatmap_*.pdf, swarm_grid_metric_by_method.*
#   Table S1             ${SCRNA}/heatmaps/dataset_summary.tsv
#   Fig. 1               not code-generated (hand-drawn scheme)
#
# Usage
#   bash entire_pipeline.sh                     # everything: pipelines + figures + copy
#   SKIP_GERMLINE=1 SKIP_ACT=1 SKIP_HG008_RUN=1 SKIP_SCRNA_RUN=1 \
#       bash entire_pipeline.sh                 # figures only, from existing results
#   RERUN_HEATMAPS=1 bash entire_pipeline.sh    # also re-render the HG008 heatmaps
#   CORES=80 DEST=/path/to/raw_figs bash entire_pipeline.sh
#
# Switches (0 = run, 1 = skip).  Heavy pipeline runs:
#   SKIP_GERMLINE_RUN, SKIP_ACT, SKIP_HG008_RUN, SKIP_SCRNA_RUN
# Figure and copy steps:
#   SKIP_FIG2 (collect + plot Fig. 2 / SI S1-S22), SKIP_FIG3 (ploidy balloon plots),
#   SKIP_FIG4 (HG008 heatmaps + montage), SKIP_FIG5 (scRNA heatmaps + summary table),
#   SKIP_COPY (copy into ${DEST})
# Extra options: RERUN_HEATMAPS (re-run the per-caller HG008 clustermap scripts even if
# their outputs exist) and RECOMPILE_SI (recompile the SI LaTeX in ${MANUSCRIPT}).
set -euo pipefail

CORES=${CORES:-200}
SCRNA_CORES=${SCRNA_CORES:-80}
# Default matches the environment created by install_step1_by_conda.sh; a missing
# environment only warns, so an already-active environment keeps working.
CONDA_ENV=${CONDA_ENV:-copy-num-bench-scwgs}
REPO=${REPO:-/stor/zxf/cnv/copy-num-bench-scwgs}
DATA=${DATA:-/stor/zxf/cnv/data}
REAL=${REAL:-/stor/zxf/cnv/real_tumor_data}
SCRNA=${SCRNA:-/nfs/wxz/zxf/cnv/copy-num-bench-scwgs-to-scrna}
MANUSCRIPT=${MANUSCRIPT:-/stor/zxf/cnv/manuscript_data}
DEST=${DEST:-${MANUSCRIPT}/raw_figs}
PREFIX=${PREFIX:-${REPO}/bench_results/bench-results-26-07-15-updated}
PLOIDY_PREFIX=${PLOIDY_PREFIX:-${PREFIX}.ploidy-plots}
HEATMAP_PDF=${HEATMAP_PDF:-${REPO}/bench_results/cnv-heatmap.pdf}
HG008_DONOR=${HG008_DONOR:-SRP047086_GIAB-HG008}
HG008_DIR=${HG008_DIR:-${REAL}/${HG008_DONOR}}
HG008_STEM=${HG008_STEM:-4from2_2_${HG008_DONOR}_3_cell-culture_ILLUMINA_nan_4_step2}

# Heavy pipeline runs.
SKIP_GERMLINE_RUN=${SKIP_GERMLINE_RUN:-0}
SKIP_ACT=${SKIP_ACT:-0}
SKIP_HG008_RUN=${SKIP_HG008_RUN:-0}
SKIP_SCRNA_RUN=${SKIP_SCRNA_RUN:-0}
# Figure / copy steps.
SKIP_FIG2=${SKIP_FIG2:-0}
SKIP_FIG3=${SKIP_FIG3:-0}
SKIP_FIG4=${SKIP_FIG4:-0}
SKIP_FIG5=${SKIP_FIG5:-0}
SKIP_COPY=${SKIP_COPY:-0}
RERUN_HEATMAPS=${RERUN_HEATMAPS:-1}
RECOMPILE_SI=${RECOMPILE_SI:-0}

log() { printf '\n=== [%s] %s ===\n' "$(date '+%F %T')" "$*"; }

# Activate the conda environment unless it is already active (CONDA_ENV="" skips this).
if [[ -n "${CONDA_ENV}" ]] && command -v conda >/dev/null 2>&1 \
        && [[ "${CONDA_DEFAULT_ENV:-}" != "${CONDA_ENV}" ]]; then
    if conda env list 2>/dev/null | awk '{print $1}' | grep -qxF "${CONDA_ENV}"; then
        eval "$(conda shell.bash hook 2>/dev/null)"
        conda activate "${CONDA_ENV}" && log "activated conda environment ${CONDA_ENV}" \
            || printf 'note: could not activate %s; using the current environment\n' \
               "${CONDA_ENV}" >&2
    else
        printf 'note: conda environment %s not found; using the current environment\n' \
            "${CONDA_ENV}" >&2
    fi
fi

# ---------------------------------------------------------------------------------
# Step 1: germline (near-haploid) benchmark -> Fig. 2, SI S1-S22, ploidy summaries
# ---------------------------------------------------------------------------------
if [[ "${SKIP_GERMLINE_RUN}" != 1 ]]; then
    log "step 1/6: germline benchmark (main.py + snakemake)"
    cd "${REPO}"
    python main.py --ploidy-tools scabsolute aneufinder chisel copynumber flcna ginkgo \
        hmmcopy sccnv scyn secnv > Snakefile
    snakemake --cores "${CORES}" --rerun-incomplete --keep-going
else
    log "step 1/6: germline pipeline skipped (SKIP_GERMLINE_RUN=1)"
fi

if [[ "${SKIP_FIG2}" != 1 ]]; then
    log "step 1/6: collecting the per-cell results and plotting Fig. 2 / SI S1-S22"
    cd "${REPO}"
    python cnv_gather_results.py -i "${DATA}"/*/4from3_*.datdir/*.perf.json -o "${PREFIX}"
    python bench_results/scWGS-performances-eval.py -t 0 -o "${PREFIX}.plots" \
        < "${PREFIX}.long.tsv"
else
    log "step 1/6: Fig. 2 / SI S1-S22 skipped (SKIP_FIG2=1)"
fi

# ---------------------------------------------------------------------------------
# Step 2: ACT breast cancers -> right panel of Fig. 3
# ---------------------------------------------------------------------------------
if [[ "${SKIP_ACT}" != 1 ]]; then
    log "step 2/6: ACT ploidy arm (main.py --tumor-fastq + snakemake)"
    cd "${REPO}"
    python main.py --tumor-fastq \
        --SraRunTable real_tumor/PRJNA629885_SraRunTable.csv_ref.tsv \
        --ploidy-file ploidy.PRJNA629885.tsv --ploidy-tools scabsolute --excluded-tools chisel \
        > real_tumor/PRJNA629885_SraRunTable_v02_with_scabsolute.snake
    snakemake -s real_tumor/PRJNA629885_SraRunTable_v02_with_scabsolute.snake \
        --cores "${CORES}" --rerun-incomplete --keep-going
else
    log "step 2/6: ACT arm skipped (SKIP_ACT=1)"
fi

# ---------------------------------------------------------------------------------
# Step 3: GIAB HG008 -> Fig. 4 and SI S23-S30
# ---------------------------------------------------------------------------------
if [[ "${SKIP_HG008_RUN}" != 1 ]]; then
    log "step 3/6: GIAB HG008 pipeline (main.py --tumor-fastq + snakemake)"
    cd "${REPO}"
    python main.py --tumor-fastq \
        --SraRunTable real_tumor_metadata/GIAB_PRJNA200694_SraRunTable_bioskryb_ILM.tsv \
        > real_tumor/GIAB_PRJNA200694_SraRunTable_bioskryb_ILM.snake
    snakemake -s real_tumor/GIAB_PRJNA200694_SraRunTable_bioskryb_ILM.snake \
        --cores "${CORES}" --rerun-incomplete --keep-going
else
    log "step 3/6: HG008 pipeline skipped (SKIP_HG008_RUN=1)"
fi

# ---------------------------------------------------------------------------------
# Step 4: Fig. 3 (ploidy balloon plots, both arms) - only plots, never re-runs a caller
# ---------------------------------------------------------------------------------
if [[ "${SKIP_FIG3}" != 1 ]]; then
    log "step 4/6: Fig. 3 (ploidy balloon plots from both arms)"
    cd "${REPO}"
    python bench_results/scWGS-ploidy-performances-eval.py \
        -i "${DATA}"/*/*ploidy*eval*summary.json \
           "${REAL}"/SRP259526_*/*ploidy*eval*summary.json \
        -o "${PLOIDY_PREFIX}"
else
    log "step 4/6: Fig. 3 skipped (SKIP_FIG3=1)"
fi

# ---------------------------------------------------------------------------------
# Step 5: Fig. 4 + SI S23-S30 (per-caller HG008 heatmaps, merged into one PDF)
# ---------------------------------------------------------------------------------
if [[ "${SKIP_FIG4}" != 1 ]]; then
    if [[ "${RERUN_HEATMAPS}" == 1 ]]; then
        log "step 5/6: re-rendering the per-caller HG008 heatmaps"
        for f in "${HG008_DIR}"/2into4_*_4_step2_*.logdir/*_tumor_clustermap.sh; do
            echo bash -evx "$f" ' || true'
        done | parallel
    fi
    log "step 5/6: merging the per-caller heatmaps into ${HEATMAP_PDF}"
    cd "${REPO}"
    python cnv_heatmap_montage.py -i "${HG008_DIR}/${HG008_STEM}"*_clustermap.pdf \
        -o "${HEATMAP_PDF}"
else
    log "step 5/6: HG008 figures skipped (SKIP_FIG4=1)"
fi

# ---------------------------------------------------------------------------------
# Step 6: Fig. 5, SI S31-S36 and Table S1 (scRNA-seq co-sequencing benchmark)
# ---------------------------------------------------------------------------------
if [[ "${SKIP_SCRNA_RUN}" != 1 ]]; then
    log "step 6/6: scRNA-seq benchmark, per dataset (companion repository)"
    cd "${SCRNA}"
    for yaml in coseq_configs/config_*.yaml; do
        snakemake --configfile config_template.yaml "${yaml}" \
            -s snakemake_pipeline/new_workflow.snake \
            --cores "${SCRNA_CORES}" --rerun-incomplete
    done
else
    log "step 6/6: scRNA-seq pipeline skipped (SKIP_SCRNA_RUN=1)"
fi

if [[ "${SKIP_FIG5}" != 1 ]]; then
    log "step 6/6: scRNA-seq figures and summary table"
    cd "${SCRNA}"
    python plot_cnv_heatmaps.py
else
    log "step 6/6: scRNA-seq figures skipped (SKIP_FIG5=1)"
fi

# ---------------------------------------------------------------------------------
# Step: copy the display items into the manuscript directory
# ---------------------------------------------------------------------------------
if [[ "${SKIP_COPY}" != 1 ]]; then
    log "copying the display items into ${DEST}"
    mkdir -p "${DEST}"

    # --- files that the SI LaTeX files include directly (exact names) --------------
    cp -v "${PREFIX}.plots-all.pdf" "${DEST}/bench-results-26-07-15-updated.plots-all.pdf"
    cp -v "${HEATMAP_PDF}" "${DEST}/cnv-heatmap.pdf"
    cp -v "${SCRNA}/heatmaps/dataset_summary.tsv" "${DEST}/dataset_summary.tsv"
    for f in \
        swarmHeatmap_Tumor_normal_classification_ROC_AUC_scRNA_aneuploidy_score_vs_scWGS_aneuploidy_status.pdf \
        swarmHeatmap_Pearson_Correlation_Coefficient_tumor_only.pdf \
        swarmHeatmap_CopyNumber_gain_ROC_AUC_tumor_only.pdf \
        swarmHeatmap_CopyNumber_loss_ROC_AUC_tumor_only.pdf \
        swarmHeatmap_Fraction_of_the_exome_with_inferred_copy_numbers.pdf \
        swarmHeatmap_Fraction_of_the_cells_with_inferred_copy_numbers.pdf ; do
        cp -v "${SCRNA}/heatmaps/${f}" "${DEST}/${f}"
    done

    # --- images of the main figures: PNG for Word, PDF for vector editing ----------
    cp -v "${PREFIX}.plots_final_grid.png" "${DEST}/Fig2_scWGS_caller_grid.png"
    cp -v "${PREFIX}.plots_final_grid.pdf" "${DEST}/Fig2_scWGS_caller_grid.pdf"
    cp -v "${PLOIDY_PREFIX}_main_four.png" "${DEST}/Fig3_ploidy_balloons.png"
    cp -v "${PLOIDY_PREFIX}_main_four.pdf" "${DEST}/Fig3_ploidy_balloons.pdf"
    cp -v "${HG008_DIR}/${HG008_STEM}"*_ginkgo_clustermap.png \
          "${DEST}/Fig4_HG008_ginkgo_heatmap.png"
    cp -v "${HG008_DIR}/${HG008_STEM}"*_ginkgo_clustermap.pdf \
          "${DEST}/Fig4_HG008_ginkgo_heatmap.pdf"
    cp -v "${SCRNA}/heatmaps/swarm_grid_metric_by_method.pdf.png" \
          "${DEST}/Fig5_scRNA_swarm_grid.png"
    cp -v "${SCRNA}/heatmaps/swarm_grid_metric_by_method.pdf" \
          "${DEST}/Fig5_scRNA_swarm_grid.pdf"
else
    log "copy step skipped (SKIP_COPY=1)"
fi

# ---------------------------------------------------------------------------------
# Optional: recompile the SI LaTeX with the freshly copied raw_figs/
# ---------------------------------------------------------------------------------
if [[ "${RECOMPILE_SI}" == 1 ]]; then
    if [[ -f "${MANUSCRIPT}/cnb-g-3-suppall-k.tex" ]]; then
        log "recompiling the SI LaTeX"
        cd "${MANUSCRIPT}"
        pdflatex -interaction=nonstopmode cnb-g-3-suppall-k.tex
        biber cnb-g-3-suppall-k
        pdflatex -interaction=nonstopmode cnb-g-3-suppall-k.tex
        pdflatex -interaction=nonstopmode cnb-g-3-suppall-k.tex
    else
        log "no cnb-g-3-suppall-k.tex in ${MANUSCRIPT}: compile the SI there manually"
    fi
fi

log "done - display items are in ${DEST}"
