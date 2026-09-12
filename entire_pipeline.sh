#!/usr/bin/env bash
#
# entire_pipeline.sh - regenerate every code-generated display item of the manuscript
# and copy the results into the manuscript directory.  The output directory is the
# default /stor/zxf/cnv/manuscript_data/raw_figs plus a code-version suffix
# (raw_figs_<commit id>-<clean|dirty> of this repository, e.g.
# raw_figs_4f7495b-dirty), so a run never overwrites the output of a different code
# version.
#
# The input FASTQ files must already be present: the germline-benchmark FASTQs in
# ../data/1from0.datdir/ (full run: scDNAaccessions.tsv; fast test run:
# scDNAaccessions.S04.tsv) and the ACT / GIAB HG008 tumor FASTQs expected by those
# arms' SraRunTables.  FASTQ downloads are deliberately NOT part of this script
# (see README.md, "Input data").
#
# What it produces (see README.md, "Generating the manuscript figures"):
#   Fig. 2 + SI S1-S22   ${PREFIX}.plots_final_grid.{png,pdf}, ${PREFIX}.plots-all.pdf
#   Fig. 3               ${PLOIDY_PREFIX}_main_four.{png,pdf} (+ per-group panels)
#   Fig. 4 + SI S23-S30  HG008 per-caller heatmaps + ${HEATMAP_PDF}
#   Fig. 5 + SI S31-S36  ${SCRNA}/heatmaps/swarmHeatmap_*.pdf, swarm_grid_metric_by_method.*
#   Table S1             ${SCRNA}/heatmaps/dataset_summary.tsv
#   Provenance           ${GIT_SNAPSHOT_TXT} (copied into ${DEST} at the very end)
#   Fig. 1               not code-generated (hand-drawn scheme)
#
# Usage
#   bash entire_pipeline.sh                     # everything: pipelines + figures + copy
#   SKIP_GERMLINE_RUN=1 SKIP_ACT=1 SKIP_HG008_RUN=1 SKIP_SCRNA_RUN=1 \
#       bash entire_pipeline.sh                 # figures only, from existing results
#   RERUN_HEATMAPS=0 bash entire_pipeline.sh    # keep the existing HG008 heatmaps (default: 1)
#   RECOMPILE_SI=1 bash entire_pipeline.sh      # also recompile the SI LaTeX after copying
#   CORES=80 DEST=/path/to/raw_figs bash entire_pipeline.sh   # the suffix is still appended
#
# Switches (0 = run, 1 = skip).  Heavy pipeline runs:
#   SKIP_GERMLINE_RUN, SKIP_ACT, SKIP_HG008_RUN, SKIP_SCRNA_RUN
# Figure and copy steps:
#   SKIP_FIG2 (collect + plot Fig. 2 / SI S1-S22), SKIP_FIG3 (ploidy balloon plots),
#   SKIP_FIG4 (HG008 heatmaps + montage), SKIP_FIG5 (scRNA heatmaps + summary table),
#   SKIP_COPY (copy into ${DEST}, including the git snapshot report)
# Extra options: RERUN_HEATMAPS=1 (default) re-runs the per-caller HG008 clustermap
# scripts unconditionally - individual failures are tolerated (e.g. SCYN), so they do
# not fail the pipeline; RERUN_HEATMAPS=0 keeps the existing heatmaps and only merges
# the montage.  RECOMPILE_SI=1 (default 0) recompiles the SI LaTeX in ${MANUSCRIPT}.
#
# Both code repositories (${REPO} and ${SCRNA}) are reported with their commit id,
# a -clean/-dirty suffix, the commit message and - when dirty - the full uncommitted
# diff, once at the very start and once at the very end of the run, so code edits made
# while the pipeline was running show up as a changed commit id or as a dirty working
# tree in the final report.  Both reports are also written to ${GIT_SNAPSHOT_TXT} and,
# at the very end, copied into ${DEST} under the same file name
# (entire_pipeline.git-snapshot.txt by default); the copy is skipped together with the
# rest of the copy step when SKIP_COPY=1.
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
# Version tag of this repository's code: <short commit id>-<clean|dirty>.  It is
# appended to the output directory name, so every run writes into a directory that
# identifies the code version that produced its display items (git_snapshot() reports
# the same tag, plus the scRNA-seq repository's tag, at the start and at the end).
REPO_COMMIT=$(git -C "${REPO}" rev-parse --short HEAD 2>/dev/null || true)
REPO_STATE=clean
if [[ -n "$(git -C "${REPO}" status --porcelain --untracked-files=no 2>/dev/null || true)" ]]; then
    REPO_STATE=dirty
fi
REPO_VERSION="${REPO_COMMIT:-unknown}-${REPO_STATE}"
DEST_BASE=${DEST:-${MANUSCRIPT}/raw_figs}_$(date +"%Y-%m-%d-%H-%M-%S")
case "${DEST_BASE}" in
    *"_${REPO_VERSION}") DEST="${DEST_BASE}" ;;
    *)                   DEST="${DEST_BASE}_${REPO_VERSION}" ;;
esac
PREFIX=${PREFIX:-${REPO}/bench_results/bench-results-26-07-15-updated}
PLOIDY_PREFIX=${PLOIDY_PREFIX:-${PREFIX}.ploidy-plots}
HEATMAP_PDF=${HEATMAP_PDF:-${REPO}/bench_results/cnv-heatmap.pdf}
GIT_SNAPSHOT_TXT=${GIT_SNAPSHOT_TXT:-${REPO}/bench_results/entire_pipeline.git-snapshot.txt}
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

# ---------------------------------------------------------------------------------
# Git provenance: both code repositories are reported with their commit id, a
# -clean/-dirty suffix, the checked-out commit's message and - when the working tree
# is dirty - the number of modified tracked files and the full uncommitted diff
# against HEAD (staged + unstaged).  The -dirty suffix reflects tracked modifications
# only; untracked paths are just counted in a note, because the result directories
# keep many of them.  git_snapshot() is called once at the very start and once at the
# very end, and appends each report to ${GIT_SNAPSHOT_TXT}, which is copied into
# ${DEST} at the very end.
# ---------------------------------------------------------------------------------
# Print the snapshot to stdout.
git_snapshot_text() {
    local phase="$1" repo label commit subject changed untracked
    log "git snapshot (${phase})"
    for repo in "${REPO}" "${SCRNA}"; do
        label="$(basename "${repo}")"
        if [[ ! -d "${repo}/.git" ]]; then
            printf '  %s: not a git repository (%s)\n' "${label}" "${repo}"
            continue
        fi
        commit="$(git -C "${repo}" rev-parse --short HEAD 2>/dev/null || true)"
        subject="$(git -C "${repo}" log -1 --format='%s' 2>/dev/null || true)"
        changed="$(git -C "${repo}" status --porcelain --untracked-files=no 2>/dev/null || true)"
        # --directory keeps this cheap on the huge result directories of the scRNA repo
        untracked="$(git -C "${repo}" ls-files --others --exclude-standard --directory 2>/dev/null | wc -l | tr -d ' ')"
        : "${commit:=unknown}"
        : "${subject:=unknown}"
        : "${untracked:=0}"
        if [[ -n "${changed}" ]]; then
            printf '  %s: %s-dirty\n' "${label}" "${commit}"
            printf '    message: %s\n' "${subject}"
            printf '    modified tracked files: %s\n' "$(printf '%s\n' "${changed}" | wc -l)"
            printf '    ----- begin uncommitted diff (git diff HEAD) -----\n'
            git -C "${repo}" diff HEAD 2>/dev/null || true
            printf '    ----- end uncommitted diff -----\n'
        else
            printf '  %s: %s-clean\n' "${label}" "${commit}"
            printf '    message: %s\n' "${subject}"
        fi
        if [[ "${untracked}" -gt 0 ]]; then
            printf '    note: %s untracked path(s) (not counted by the dirty marker)\n' "${untracked}"
        fi
    done
}

# Print the snapshot and append it to ${GIT_SNAPSHOT_TXT} (the report that is copied
# into ${DEST} at the very end).  A failure to write the report only warns.
git_snapshot() {
    local phase="$1" report
    report="$(git_snapshot_text "${phase}")"
    printf '%s\n' "${report}"
    printf '%s\n' "${report}" >> "${GIT_SNAPSHOT_TXT}" 2>/dev/null \
        || printf 'warning: could not append the git snapshot report to %s\n' \
                  "${GIT_SNAPSHOT_TXT}" >&2
}

# Activate the conda environment unless it is already active (CONDA_ENV="" or no
# conda in PATH skips this; a missing environment only warns).
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

mkdir -p "$(dirname "${GIT_SNAPSHOT_TXT}")" 2>/dev/null || true
: > "${GIT_SNAPSHOT_TXT}" 2>/dev/null \
    || printf 'warning: could not create the git snapshot report %s\n' "${GIT_SNAPSHOT_TXT}" >&2
git_snapshot "start"

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
# Step 5: Fig. 4 + SI S23-S30 (re-render the per-caller HG008 heatmaps when
# RERUN_HEATMAPS=1, then merge them into one PDF)
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
git_snapshot "end"

# The git report (start + end snapshots) joins the other display items at the very end.
if [[ "${SKIP_COPY}" != 1 ]]; then
    if [[ -f "${GIT_SNAPSHOT_TXT}" ]]; then
        mkdir -p "${DEST}"
        cp -v "${GIT_SNAPSHOT_TXT}" "${DEST}/$(basename "${GIT_SNAPSHOT_TXT}")"
    else
        printf 'warning: %s does not exist; the git snapshot report is not copied\n' \
               "${GIT_SNAPSHOT_TXT}" >&2
    fi
else
    log "git snapshot report not copied (SKIP_COPY=1); it is in ${GIT_SNAPSHOT_TXT}"
fi
