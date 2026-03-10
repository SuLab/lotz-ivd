#!/usr/bin/env python3
"""Generate MANUSCRIPT.pdf with embedded figures and activated GitHub hyperlinks."""

import os
import re
import markdown
from xhtml2pdf import pisa

REPO_BASE = "https://github.com/andrewsu/lotz-ivd/blob/main"
RAW_BASE = "https://raw.githubusercontent.com/andrewsu/lotz-ivd/main"
RESULTS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(RESULTS_DIR)

def resolve_image_paths(html, base_dir):
    """Convert relative image src to absolute local file:// paths for PDF embedding."""
    def replace_src(m):
        src = m.group(1)
        if src.startswith(("http://", "https://")):
            return m.group(0)
        abs_path = os.path.abspath(os.path.join(base_dir, src))
        if os.path.exists(abs_path):
            return f'src="{abs_path}"'
        print(f"  WARNING: image not found: {abs_path}")
        return m.group(0)
    return re.sub(r'src="([^"]+)"', replace_src, html)

def resolve_hyperlinks(html):
    """Convert relative href links to full GitHub URLs."""
    def replace_href(m):
        href = m.group(1)
        # Skip anchors and already-absolute URLs
        if href.startswith(("#", "http://", "https://", "mailto:")):
            return m.group(0)
        # ../phylo_analysis/... -> repo blob link
        if href.startswith("../"):
            rel = href[3:]  # strip ../
            # Remove anchor for GitHub blob URLs, keep as separate anchor
            if "#" in rel:
                path, anchor = rel.split("#", 1)
                return f'href="{REPO_BASE}/{path}#{anchor}"'
            return f'href="{REPO_BASE}/{rel}"'
        # relative within results/
        return f'href="{REPO_BASE}/results/{href}"'
    return re.sub(r'href="([^"]+)"', replace_href, html)

def main():
    md_path = os.path.join(RESULTS_DIR, "MANUSCRIPT.md")
    pdf_path = os.path.join(RESULTS_DIR, "MANUSCRIPT.pdf")

    with open(md_path, "r") as f:
        md_text = f.read()

    # Convert markdown to HTML
    html_body = markdown.markdown(
        md_text,
        extensions=["tables", "toc", "fenced_code"]
    )

    # Add <a name="..."> anchors for xhtml2pdf internal linking
    # xhtml2pdf doesn't resolve id= attributes as link targets
    html_body = re.sub(
        r'<(h[1-6])\s+id="([^"]+)"',
        r'<\1 id="\2"><a name="\2"></a',
        html_body,
    )

    # Resolve paths
    html_body = resolve_image_paths(html_body, RESULTS_DIR)
    html_body = resolve_hyperlinks(html_body)

    html_doc = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<style>
  @page {{ size: A4; margin: 1.8cm; }}
  body {{
    font-family: Helvetica, Arial, sans-serif;
    font-size: 10pt;
    line-height: 1.45;
    color: #222;
  }}
  h1 {{
    font-size: 16pt;
    margin-top: 0.6em;
    color: #1a1a1a;
  }}
  h2 {{
    font-size: 14pt;
    margin-top: 1.2em;
    border-bottom: 1px solid #ccc;
    padding-bottom: 0.2em;
    color: #2c3e50;
  }}
  h3 {{
    font-size: 12pt;
    margin-top: 1em;
    color: #34495e;
  }}
  h4 {{
    font-size: 10pt;
    margin-top: 0.8em;
  }}
  table {{
    border-collapse: collapse;
    width: 100%;
    margin: 0.8em 0;
    font-size: 8pt;
  }}
  th, td {{
    border: 1px solid #bbb;
    padding: 3px 6px;
    text-align: left;
  }}
  th {{
    background-color: #3498db;
    color: white;
    font-weight: bold;
  }}
  tr:nth-child(even) {{
    background-color: #f8f9fa;
  }}
  blockquote {{
    border-left: 4px solid #4a86c8;
    margin: 0.8em 0;
    padding: 0.5em 1em;
    background: #f4f8fc;
    font-size: 9pt;
  }}
  code {{
    font-family: Courier;
    font-size: 8pt;
    background: #f5f5f5;
    padding: 1px 3px;
  }}
  hr {{
    border: none;
    border-top: 1px solid #ddd;
    margin: 1.2em 0;
  }}
  a {{
    color: #2a6496;
    text-decoration: underline;
  }}
  p {{
    margin: 0.4em 0;
  }}
  ul, ol {{
    margin: 0.4em 0;
  }}
  li {{
    margin: 0.2em 0;
  }}
  img {{
    max-width: 100%;
    margin: 8px 0;
  }}
  strong {{
    color: #111;
  }}
</style>
</head>
<body>
{html_body}
</body>
</html>"""

    with open(pdf_path, "wb") as f:
        status = pisa.CreatePDF(html_doc, dest=f, encoding="utf-8")

    if status.err:
        print(f"ERROR: PDF generation failed with {status.err} errors")
    else:
        size_kb = os.path.getsize(pdf_path) / 1024
        print(f"Created {pdf_path} ({size_kb:.0f} KB)")

if __name__ == "__main__":
    main()
