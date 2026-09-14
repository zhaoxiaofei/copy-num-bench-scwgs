#!/usr/bin/env python3
"""Combine the single-page per-caller CNV heatmap PDFs into one multi-page PDF.

``cnv_clustermap.py`` writes one heatmap per caller
(``<prefix>_<caller>_clustermap.pdf`` plus a ``.png`` sibling).  The manuscript shows
the HG008 heatmaps as ONE multi-page supplementary figure (Supplementary Figs. S23-S30,
``raw_figs/cnv-heatmap.pdf``) and re-uses the Ginkgo page as main Fig. 4.  That montage
used to be made by hand with ImageMagick; this script rebuilds it deterministically
from the per-caller PDFs, so the figure can be regenerated from pipeline output alone.

The merge itself uses pypdf when it is installed, otherwise the poppler tool
``pdfunite``, and finally Ghostscript, so no Python package is required.
The per-caller heatmaps are vector PDFs of a very large grid (one is ~110 MB for the
50 kb bins of the HG008 run), so a lossless merge is ~800 MB.  There are two ways to
get a small supplementary PDF instead: point the script at the ``.png`` siblings that
``cnv_clustermap.py`` writes next to each PDF (fast, no extra tools), or pass
``--rasterize-dpi 300`` to re-render the PDF pages.  Both reproduce what the previous
ImageMagick montage did.

Usage
-----
    python cnv_heatmap_montage.py \
        -i '/stor/zxf/cnv/real_tumor_data/GIAB-HG008/4from2_2_*_4_step2_*_clustermap.pdf' \
        -o raw_figs/cnv-heatmap.pdf

    python cnv_heatmap_montage.py \
        -i '/stor/zxf/cnv/real_tumor_data/GIAB-HG008/4from2_2_*_4_step2_*_clustermap.pdf' \
        -o raw_figs/cnv-heatmap.pdf --rasterize-dpi 300

    python cnv_heatmap_montage.py \
        -i '/stor/zxf/cnv/real_tumor_data/GIAB-HG008/4from2_2_*_4_step2_*_clustermap.png' \
        -o raw_figs/cnv-heatmap.pdf

``--order`` fixes the page order by caller key.  By default the pages follow the
callers' exact publication dates (HMMcopy 2006, Copynumber 2012, Ginkgo 2015,
AneuFinder 2016, SCCNV 2020, CHISEL 2021, SCYN 2021, SeCNV 2022, FLCNA 2024),
matching the caller order of the main figures; pass ``--order`` to override it,
for example

    --order ginkgo chisel copynumber flcna hmmcopy sccnv secnv aneufinder

Callers that are not named in ``--order`` keep the alphabetical order of their file
names and are appended last.  The ``<prefix>_clustermap_relative.pdf`` variants are
skipped unless ``--include-relative`` is given, so the usual glob (``*_clustermap.pdf``)
does not accidentally double the page count.
"""

import argparse
import glob
import os
import shutil
import subprocess
import sys
import tempfile
from typing import List, Sequence


def eprint(msg: str) -> None:
    sys.stderr.write(str(msg) + "\n")
    sys.stderr.flush()


def expand_patterns(patterns: Sequence[str]) -> List[str]:
    """Expand shell-style globs; plain paths are kept as they are."""
    files: List[str] = []
    for pat in patterns:
        hits = sorted(glob.glob(pat)) if any(c in pat for c in "*?[") else [pat]
        files.extend(hits)
    return files


IMAGE_EXTS = (".png", ".jpg", ".jpeg", ".tif", ".tiff", ".bmp")

# Default page order: callers sorted by their exact publication date (see
# CALLER_PUBLICATION in bench_results/scWGS-performances-eval.py).  Callers not
# listed here keep the alphabetical order of their file names and are appended last.
DEFAULT_CALLER_ORDER = [
    "hmmcopy", "copynumber", "ginkgo", "aneufinder", "sccnv",
    "chisel", "scyn", "secnv", "flcna", "scabsolute",
]


def is_image(path: str) -> bool:
    return os.path.splitext(path)[1].lower() in IMAGE_EXTS


def caller_key(path: str) -> str:
    """Caller name of a ``<prefix>_<caller>_clustermap[_relative].pdf`` file."""
    base = os.path.splitext(os.path.basename(path))[0]
    if base.endswith("_clustermap_relative"):
        stem, suffix = base[: -len("_clustermap_relative")], " (relative)"
    elif base.endswith("_clustermap"):
        stem, suffix = base[: -len("_clustermap")], ""
    else:
        return base
    return stem.rsplit("_", 1)[-1] + suffix


def sort_files(files: Sequence[str], order: Sequence[str]) -> List[str]:
    """Page order: the caller keys named in ``order`` first, then the rest (A-Z)."""
    order = [str(x) for x in order]

    def key(path: str):
        caller = caller_key(path)
        if caller in order:
            return (0, order.index(caller), path)
        return (1, len(order), caller.lower())

    return sorted(files, key=key)


def merge_pdfs(files: Sequence[str], output: str) -> str:
    """Merge ``files`` into one PDF; returns the backend that was used."""
    try:  # optional, and by far the cleanest backend
        import pypdf  # type: ignore

        writer = pypdf.PdfWriter()
        for path in files:
            writer.append(path)
        with open(output, "wb") as fh:
            writer.write(fh)
        return "pypdf"
    except ImportError:
        pass

    if shutil.which("pdfunite"):
        subprocess.run(["pdfunite", *files, output], check=True)
        return "pdfunite"

    if shutil.which("gs"):
        subprocess.run(
            ["gs", "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=pdfwrite",
             F"-sOutputFile={output}", *files],
            check=True,
        )
        return "ghostscript"

    raise SystemExit(
        "No PDF merge backend found: install pypdf, poppler-utils (pdfunite) or ghostscript."
    )


def page_count(path: str):
    """Number of pages reported by pdfinfo, or None when it is unavailable."""
    if not shutil.which("pdfinfo"):
        return None
    try:
        out = subprocess.run(["pdfinfo", path], capture_output=True, text=True, check=True).stdout
    except subprocess.CalledProcessError:
        return None
    for line in out.splitlines():
        if line.startswith("Pages:"):
            try:
                return int(line.split()[1])
            except (IndexError, ValueError):
                return None
    return None


def render_pdf_pages(path: str, dpi: int, outdir: str) -> List[str]:
    """Render every page of one PDF to PNG; returns the image paths."""
    if not shutil.which("pdftoppm"):
        raise SystemExit("Rendering PDF pages needs poppler's pdftoppm.")
    prefix = os.path.join(outdir, "page")
    subprocess.run(["pdftoppm", "-r", str(dpi), "-png", path, prefix], check=True)
    return sorted(glob.glob(prefix + "*.png"))


def images_to_pdf(images: Sequence[str], output: str) -> int:
    """Write one PDF page per image: the .png siblings, or rasterised PDF pages."""
    from PIL import Image  # matplotlib depends on pillow, so it is always available

    pages = [Image.open(path).convert("RGB") for path in images]
    if not pages:
        raise SystemExit("No page to write.")
    try:
        dpi = float(pages[0].info.get("dpi", (300, 300))[0]) or 300.0
    except (TypeError, ValueError, IndexError):
        dpi = 300.0
    pages[0].save(output, save_all=True, append_images=pages[1:], resolution=dpi)
    for page in pages:
        page.close()
    return len(pages)


def build_raster_pdf(files: Sequence[str], output: str, dpi: int) -> int:
    """One PDF page per input page: images are used as they are, PDFs are rendered at ``dpi``."""
    with tempfile.TemporaryDirectory(prefix="cnv_heatmap_montage.") as tmp:
        images: List[str] = []
        for i, path in enumerate(files):
            if is_image(path):
                images.append(path)
                continue
            sub = os.path.join(tmp, F"in{i:03d}")
            os.makedirs(sub, exist_ok=True)
            images.extend(render_pdf_pages(path, dpi, sub))
        return images_to_pdf(images, output)


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Merge the per-caller CNV heatmap PDFs into one multi-page PDF "
                    "(materials for main Fig. 4 / Supplementary Figs. S23-S30).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-i", "--input", nargs="+", required=True,
                   help="Per-caller heatmap PDFs or shell-style glob patterns")
    p.add_argument("-o", "--output", required=True,
                   help="Merged multi-page PDF to write")
    p.add_argument("--order", nargs="*", default=[], metavar="CALLER",
                   help="Caller keys defining the page order; unnamed callers are "
                        "appended in alphabetical order")
    p.add_argument("--include-relative", action="store_true",
                   help="Also merge the <caller>_clustermap_relative.pdf variants")
    p.add_argument("--rasterize-dpi", type=int, default=0, metavar="DPI",
                   help="Render every page at this DPI and write a raster multi-page "
                        "PDF (e.g. 300); 0 merges the vector PDFs losslessly")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)

    files = expand_patterns(args.input)
    if not args.include_relative:
        files = [f for f in files if "_clustermap_relative" not in os.path.basename(f)]
    files = [f for f in files if os.path.isfile(f)]
    if not files:
        eprint("No input heatmap PDF found.")
        return 1

    files = sort_files(files, args.order or DEFAULT_CALLER_ORDER)
    for path in files:
        n = page_count(path)
        eprint(F"  page {caller_key(path)}: {path}"
               + (F" ({n} page(s))" if n is not None else ""))

    if args.rasterize_dpi > 0:
        n_pages = build_raster_pdf(files, args.output, args.rasterize_dpi)
        eprint(F"Rasterized {len(files)} file(s) at {args.rasterize_dpi} dpi -> {args.output} "
               F"({n_pages} pages)")
    elif all(is_image(path) for path in files):
        n_pages = images_to_pdf(files, args.output)
        eprint(F"Assembled {n_pages} image page(s) -> {args.output}")
    else:
        backend = merge_pdfs(files, args.output)
        n_out = page_count(args.output)
        eprint(F"Merged {len(files)} file(s) with {backend} -> {args.output}"
               + (F" ({n_out} pages)" if n_out is not None else ""))
        size_mb = os.path.getsize(args.output) / 1e6
        if size_mb > 50:
            eprint(F"Note: the merged file is {size_mb:.0f} MB because the heatmaps are vector "
                   "grids; re-run with --rasterize-dpi 300 for a small supplementary PDF.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
