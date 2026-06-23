#!/usr/bin/env python3
"""Build the combined manuscript PDF from its markdown source.

Pipeline (reproduces the recipe used in commits 45fc27f / 45df855 / 8f46532):
  1. pandoc  : markdown -> standalone HTML, styled by manuscript_style.css
  2. weasyprint : HTML -> PDF (resolves CSS + figure PNGs relative to docs/)
  3. PyMuPDF : recompress embedded images to JPEG q88 to control file size

Run from anywhere:  python3 scripts/build_combined_pdf.py
"""
from __future__ import annotations
import subprocess
from pathlib import Path

DOCS = Path(__file__).resolve().parent.parent / "docs"
STEM = "IVD_MANUSCRIPT_2026-06-12_combined"
MD = DOCS / f"{STEM}.md"
HTML = DOCS / f"{STEM}.html"
PDF = DOCS / f"{STEM}.pdf"
TITLE = "A Continuum-Aware Single-Cell Atlas of the Human Intervertebral Disc"
JPEG_QUALITY = 88


def build_html() -> None:
    # cwd=DOCS so the relative `manuscript_style.css` link matches the source tree.
    subprocess.run(
        [
            "pandoc", MD.name, "-s",
            "--metadata", f"title={TITLE}",
            "-c", "manuscript_style.css",
            "-o", HTML.name,
        ],
        cwd=DOCS,
        check=True,
    )
    print(f"[html] wrote {HTML}")


def build_pdf() -> None:
    from weasyprint import HTML as WeasyHTML

    # filename= sets the base URL to the HTML's directory, so the relative
    # stylesheet and manuscript_figures/*.png references resolve.
    WeasyHTML(filename=str(HTML)).write_pdf(str(PDF))
    print(f"[pdf] wrote {PDF} ({PDF.stat().st_size:,} bytes, pre-recompress)")


def recompress(quality: int = JPEG_QUALITY) -> None:
    import fitz

    doc = fitz.open(str(PDF))
    n_pages = doc.page_count
    doc.rewrite_images(quality=quality, lossy=True, lossless=True, color=True, gray=True)
    tmp = PDF.with_suffix(".recompress.pdf")
    doc.save(str(tmp), garbage=4, deflate=True, incremental=False)
    doc.close()
    tmp.replace(PDF)
    print(f"[pdf] recompressed images to JPEG q{quality}: "
          f"{PDF.stat().st_size:,} bytes, {n_pages} pages")


if __name__ == "__main__":
    build_html()
    build_pdf()
    recompress()
