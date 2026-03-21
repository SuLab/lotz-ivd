#!/usr/bin/env python3
"""Module 05: Integration orchestrator for IVD atlas.

Orchestrates three parallel integration workflows:
  A. CCA (R, primary) — label-free canonical correlation analysis
  B. scANVI (Python) — semi-supervised deep generative model
  C. STACAS (R) — semi-supervised anchor-based integration

After workflows complete, generates shared outputs (inclusion summary, study
caveats, workflow comparison) and runs validation.

After the human checkpoint, use --select-workflow to promote the chosen
workflow's outputs to the canonical data/integrated/ path for Modules 06-12.

Usage:
    python3 scripts/05_integration.py                              # All workflows, all objects
    python3 scripts/05_integration.py --workflow scanvi             # Single workflow
    python3 scripts/05_integration.py --workflow cca --object NP    # Single workflow + object
    python3 scripts/05_integration.py --workflow all --object NP    # All workflows, one object
    python3 scripts/05_integration.py --validate-only               # Validation only
    python3 scripts/05_integration.py --force                       # Re-run even if outputs exist
    python3 scripts/05_integration.py --compare-only                # Run comparison only
    python3 scripts/05_integration.py --select-workflow cca         # Promote CCA to canonical path
    python3 scripts/05_integration.py --select-workflow scanvi      # Promote scANVI to canonical path
"""

import subprocess
import sys
import warnings
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
from jinja2 import Template

warnings.filterwarnings('ignore', category=FutureWarning)

# ── Paths ────────────────────────────────────────────────────────────────────
BASE = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = BASE / "scripts"
INT_BASE = BASE / "data" / "integrated"
RESULTS_DIR = BASE / "results" / "integration"
META_DIR = BASE / "metadata"

for d in [INT_BASE, RESULTS_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# Workflow subdirectories
for wf in ["cca", "scanvi", "stacas"]:
    (INT_BASE / wf).mkdir(parents=True, exist_ok=True)

OBJECTS = ["NP", "AF", "CEP", "all_cells"]


# ═════════════════════════════════════════════════════════════════════════════
# WORKFLOW RUNNERS
# ═════════════════════════════════════════════════════════════════════════════

def run_workflow_cca(object_name=None, force=False):
    """Run Workflow A: CCA integration (R)."""
    print("\n" + "=" * 60)
    print("Running Workflow A: CCA Integration (R)")
    print("=" * 60)

    cmd = ["Rscript", str(SCRIPTS_DIR / "05a_integration_cca.R")]
    if object_name:
        cmd += ["--object", object_name]
    if force:
        cmd.append("--force")

    print(f"  Command: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=False)

    if result.returncode != 0:
        print(f"  WARNING: CCA workflow exited with code {result.returncode}")
        return False
    return True


def run_workflow_scanvi(object_name=None, force=False):
    """Run Workflow B: scANVI integration (Python)."""
    print("\n" + "=" * 60)
    print("Running Workflow B: scANVI Integration (Python)")
    print("=" * 60)

    cmd = ["python3", str(SCRIPTS_DIR / "05b_integration_scanvi.py")]
    if object_name:
        cmd += ["--object", object_name]
    if force:
        cmd.append("--force")

    print(f"  Command: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=False)

    if result.returncode != 0:
        print(f"  WARNING: scANVI workflow exited with code {result.returncode}")
        return False
    return True


def run_workflow_stacas(object_name=None, force=False):
    """Run Workflow C: STACAS integration (R)."""
    print("\n" + "=" * 60)
    print("Running Workflow C: STACAS Integration (R)")
    print("=" * 60)

    cmd = ["Rscript", str(SCRIPTS_DIR / "05c_integration_stacas.R")]
    if object_name:
        cmd += ["--object", object_name]
    if force:
        cmd.append("--force")

    print(f"  Command: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=False)

    if result.returncode != 0:
        print(f"  WARNING: STACAS workflow exited with code {result.returncode}")
        return False
    return True


def run_comparison(object_name=None):
    """Run the comparison script across all workflows."""
    print("\n" + "=" * 60)
    print("Running Workflow Comparison")
    print("=" * 60)

    cmd = ["python3", str(SCRIPTS_DIR / "05d_integration_comparison.py")]
    if object_name:
        cmd += ["--object", object_name]

    print(f"  Command: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=False)

    if result.returncode != 0:
        print(f"  WARNING: Comparison script exited with code {result.returncode}")
        return False
    return True


# ═════════════════════════════════════════════════════════════════════════════
# SHARED OUTPUTS: INCLUSION SUMMARY & STUDY CAVEATS
# ═════════════════════════════════════════════════════════════════════════════

def generate_inclusion_summary():
    """Generate inclusion_summary.tsv and .html for manuscript supplement.

    Scans all workflow output directories for integrated objects and aggregates
    cell counts per study per object.
    """
    import scanpy as sc

    print("\n  Generating inclusion summary...")

    meta_path = META_DIR / "sample_metadata.tsv"
    meta = pd.read_csv(meta_path, sep='\t') if meta_path.exists() else None

    rows = []
    # Use scanvi outputs (h5ad) as the primary source since we can read them
    # directly. Fall back to checking other workflows for completeness.
    for obj_name in OBJECTS:
        adata = None

        # Try scanvi first (h5ad, directly loadable)
        scanvi_path = INT_BASE / "scanvi" / f"{obj_name}.h5ad"
        if scanvi_path.exists():
            try:
                adata = sc.read_h5ad(scanvi_path, backed='r')
            except Exception:
                pass

        if adata is None:
            # No loadable integration output for this object
            continue

        for study in sorted(adata.obs['study'].unique()):
            study_mask = adata.obs['study'] == study
            n_cells = int(study_mask.sum())
            n_samples = adata.obs.loc[study_mask, 'sample_id'].nunique() if 'sample_id' in adata.obs.columns else 0

            first_author = ""
            year = ""
            compartment = ""
            conditions = ""
            platform = ""
            if meta is not None:
                study_meta = meta[meta['study_accession'] == study]
                if len(study_meta) > 0:
                    first_author = study_meta.iloc[0].get('first_author', '')
                    year = str(study_meta.iloc[0].get('year', ''))
                    compartment = ', '.join(study_meta['compartment'].unique())
                    conditions = ', '.join(study_meta['condition_harmonized'].dropna().unique())
                    platform = study_meta.iloc[0].get('sequencing_platform', '')

            rows.append({
                'object': obj_name,
                'accession': study,
                'first_author': first_author,
                'year': year,
                'n_samples': n_samples,
                'n_cells': n_cells,
                'compartment': compartment,
                'conditions': conditions,
                'platform': platform,
            })
        adata.file.close()

    if not rows:
        print("  WARNING: No integration outputs found for inclusion summary")
        return

    df = pd.DataFrame(rows)
    df.to_csv(RESULTS_DIR / "inclusion_summary.tsv", sep='\t', index=False)

    # HTML version
    html = df.to_html(index=False, classes='summary-table')
    html_page = f"""<!DOCTYPE html>
<html><head><title>Inclusion Summary</title>
<style>
body {{ font-family: Arial; max-width: 1200px; margin: 0 auto; padding: 20px; }}
.summary-table {{ border-collapse: collapse; width: 100%; }}
.summary-table th {{ background: #3498db; color: white; padding: 8px; }}
.summary-table td {{ border: 1px solid #ddd; padding: 6px; }}
.summary-table tr:nth-child(even) {{ background: #f2f2f2; }}
</style></head><body>
<h1>Study Inclusion Summary</h1>
<p>Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}</p>
{html}
</body></html>"""
    (RESULTS_DIR / "inclusion_summary.html").write_text(html_page)
    print(f"  Saved: inclusion_summary.tsv and .html")


def generate_study_caveats():
    """Generate study_caveats.tsv for manuscript supplement."""
    caveats = [
        {"study": "GSE165722", "first_author": "Tu 2022",
         "caveat": "BD Rhapsody platform (not 10x)",
         "impact": "Different capture efficiency, gene detection",
         "mitigation": "Platform-aware batch correction via study-level batch key"},
        {"study": "GSE205535", "first_author": "Li 2022",
         "caveat": "BD Rhapsody platform; published corrigenda",
         "impact": "See above; potential data quality issues",
         "mitigation": "Monitor for outlier behavior in integration"},
        {"study": "CNP0002664", "first_author": "Han 2022",
         "caveat": "Singleron Matrix platform (not 10x)",
         "impact": "Different capture efficiency",
         "mitigation": "Same as above"},
        {"study": "GSE242443", "first_author": "Kuchynsky 2024",
         "caveat": "Culture-expanded CEP cells",
         "impact": "Culture alters cell states; may not reflect in vivo biology",
         "mitigation": "Caveat in all CEP results; compare with non-expanded CEP from GSE160756"},
        {"study": "GSE255768", "first_author": "Shi 2024",
         "caveat": "Degenerative endplate only; no healthy control",
         "impact": "Cannot do healthy vs. degenerated comparison for this study alone",
         "mitigation": "Healthy CEP baseline from GSE160756"},
        {"study": "GSE230809", "first_author": "Swahn 2024",
         "caveat": "All-male donors; age-disease confounded",
         "impact": "Cannot separate age from degeneration effects",
         "mitigation": "Note in interpretation; sex-specific effects cannot be assessed"},
        {"study": "GSE205535 NNP", "first_author": "Li 2022",
         "caveat": "11yo spinal cord injury, classified as 'healthy'",
         "impact": "Trauma confound",
         "mitigation": "Excluded from DE comparisons"},
        {"study": "GSE189916", "first_author": "Jiang 2022",
         "caveat": "Whole IVD (compartments not separated)",
         "impact": "Cannot assign cells to NP/AF/CEP",
         "mitigation": "Included only in all-cells object"},
    ]
    df = pd.DataFrame(caveats)
    df.to_csv(RESULTS_DIR / "study_caveats.tsv", sep='\t', index=False)
    print(f"  Saved: study_caveats.tsv")


# ═════════════════════════════════════════════════════════════════════════════
# VALIDATION
# ═════════════════════════════════════════════════════════════════════════════

def validate_all():
    """Run validation across all workflows."""
    print("\n" + "=" * 60)
    print("Validating integration outputs (all workflows)")
    print("=" * 60)

    messages = []
    all_pass = True

    # Check per-workflow outputs
    workflow_files = {
        "cca": (".rds", INT_BASE / "cca"),
        "scanvi": (".h5ad", INT_BASE / "scanvi"),
        "stacas": (".rds", INT_BASE / "stacas"),
    }

    at_least_one_per_object = {obj: False for obj in OBJECTS}

    for wf_name, (ext, wf_dir) in workflow_files.items():
        for obj_name in OBJECTS:
            path = wf_dir / f"{obj_name}{ext}"
            # Also check alternative extension
            alt_ext = '.h5ad' if ext == '.rds' else '.rds'
            alt_path = wf_dir / f"{obj_name}{alt_ext}"

            if path.exists() or alt_path.exists():
                messages.append(f"PASS: {wf_name}/{obj_name} output exists")
                at_least_one_per_object[obj_name] = True
            else:
                messages.append(f"INFO: {wf_name}/{obj_name} output missing")

    # Check that at least one workflow succeeded per object
    for obj_name, has_any in at_least_one_per_object.items():
        if has_any:
            messages.append(f"PASS: At least one workflow completed for {obj_name}")
        else:
            messages.append(f"FAIL: No workflow completed for {obj_name}")
            all_pass = False

    # Check shared outputs
    for fname in ["inclusion_summary.tsv", "study_caveats.tsv"]:
        if (RESULTS_DIR / fname).exists():
            messages.append(f"PASS: {fname} exists")
        else:
            messages.append(f"FAIL: {fname} missing")
            all_pass = False

    # Check comparison outputs
    for fname in ["workflow_comparison.tsv", "workflow_comparison_report.html"]:
        if (RESULTS_DIR / fname).exists():
            messages.append(f"PASS: {fname} exists")
        else:
            messages.append(f"WARNING: {fname} missing (run comparison)")

    for msg in messages:
        print(f"  {msg}")

    print(f"\n{'='*60}")
    print(f"OVERALL: {'ALL CHECKS PASSED' if all_pass else 'SOME CHECKS FAILED'}")
    print(f"{'='*60}")
    return all_pass, messages


# ═════════════════════════════════════════════════════════════════════════════
# WORKFLOW SELECTION (post-checkpoint)
# ═════════════════════════════════════════════════════════════════════════════

def select_workflow(workflow_name, object_name=None):
    """Promote a selected workflow's outputs to the canonical data/integrated/ path.

    After the human checkpoint chooses a workflow, this copies (scANVI) or
    converts (CCA/STACAS RDS→h5ad) the workflow outputs so that Modules 06-12
    can load from the standard data/integrated/{object}.h5ad path.
    """
    import shutil

    print("\n" + "=" * 60)
    print(f"Selecting workflow: {workflow_name}")
    print(f"Promoting outputs to canonical path: {INT_BASE}/")
    print("=" * 60)

    objects = [object_name] if object_name else OBJECTS

    if workflow_name == "scanvi":
        # scANVI outputs are already h5ad — copy directly
        for obj in objects:
            src = INT_BASE / "scanvi" / f"{obj}.h5ad"
            dst = INT_BASE / f"{obj}.h5ad"
            if not src.exists():
                print(f"  WARNING: {src} not found, skipping {obj}")
                continue
            print(f"  Copying scanvi/{obj}.h5ad -> {obj}.h5ad")
            shutil.copy2(src, dst)

    elif workflow_name in ("cca", "stacas"):
        # CCA/STACAS outputs are RDS — convert via R helper
        cmd = ["Rscript", str(SCRIPTS_DIR / "05_select_workflow.R"),
               "--workflow", workflow_name]
        if object_name:
            cmd += ["--object", object_name]

        print(f"  Converting RDS -> h5ad via: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=False)

        if result.returncode != 0:
            print(f"  ERROR: RDS->h5ad conversion failed (exit code {result.returncode})")
            sys.exit(1)

    # Verify outputs exist
    print("\n  Verifying canonical outputs:")
    all_ok = True
    for obj in objects:
        canonical = INT_BASE / f"{obj}.h5ad"
        if canonical.exists():
            size_mb = canonical.stat().st_size / (1024 * 1024)
            print(f"    OK: {obj}.h5ad ({size_mb:.0f} MB)")
        else:
            print(f"    MISSING: {obj}.h5ad")
            all_ok = False

    if all_ok:
        print(f"\n  Workflow '{workflow_name}' selected successfully.")
        print(f"  Modules 06-12 will now use these outputs.")
    else:
        print(f"\n  WARNING: Some outputs missing — check errors above.")

    # Write a record of which workflow was selected
    record_path = INT_BASE / "selected_workflow.txt"
    record_path.write_text(
        f"selected_workflow: {workflow_name}\n"
        f"selected_at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
        f"objects: {', '.join(objects)}\n"
    )
    print(f"  Recorded selection in {record_path}")


# ═════════════════════════════════════════════════════════════════════════════
# MAIN
# ═════════════════════════════════════════════════════════════════════════════

def main():
    args = sys.argv[1:]

    validate_only = '--validate-only' in args
    compare_only = '--compare-only' in args
    force = '--force' in args

    # Parse --select-workflow
    select_wf = None
    for i, a in enumerate(args):
        if a == '--select-workflow' and i + 1 < len(args):
            select_wf = args[i + 1]

    # Parse --workflow
    workflow = 'all'
    for i, a in enumerate(args):
        if a == '--workflow' and i + 1 < len(args):
            workflow = args[i + 1]

    # Parse --object
    object_name = None
    for i, a in enumerate(args):
        if a == '--object' and i + 1 < len(args):
            object_name = args[i + 1]

    # Validate workflow choice
    valid_workflows = {'all', 'cca', 'scanvi', 'stacas'}
    if workflow not in valid_workflows:
        print(f"ERROR: Unknown workflow '{workflow}'. Choose from: {valid_workflows}")
        sys.exit(1)

    if select_wf and select_wf not in {'cca', 'scanvi', 'stacas'}:
        print(f"ERROR: Unknown workflow '{select_wf}'. Choose from: cca, scanvi, stacas")
        sys.exit(1)

    # Validate object choice
    if object_name and object_name not in OBJECTS:
        print(f"ERROR: Unknown object '{object_name}'. Choose from: {OBJECTS}")
        sys.exit(1)

    # Handle --select-workflow (post-checkpoint action)
    if select_wf:
        select_workflow(select_wf, object_name)
        sys.exit(0)

    if validate_only:
        passed, _ = validate_all()
        sys.exit(0 if passed else 1)

    if compare_only:
        run_comparison(object_name)
        sys.exit(0)

    print("=" * 60)
    print("Module 05: Integration (Orchestrator)")
    print(f"Workflow: {workflow}")
    print(f"Object: {object_name or 'all'}")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)

    # Run selected workflow(s)
    results = {}

    if workflow in ('all', 'cca'):
        results['cca'] = run_workflow_cca(object_name, force)

    if workflow in ('all', 'scanvi'):
        results['scanvi'] = run_workflow_scanvi(object_name, force)

    if workflow in ('all', 'stacas'):
        results['stacas'] = run_workflow_stacas(object_name, force)

    # Print workflow results
    print("\n" + "=" * 60)
    print("Workflow Results:")
    for wf, success in results.items():
        status = "SUCCESS" if success else "FAILED"
        print(f"  {wf:8s}: {status}")
    print("=" * 60)

    # Generate shared outputs
    print("\nGenerating shared outputs...")
    generate_inclusion_summary()
    generate_study_caveats()

    # Run comparison across all workflows
    run_comparison(object_name)

    # Run validation
    passed, val_messages = validate_all()

    print(f"\nCompleted: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Overall: {'PASSED' if passed else 'FAILED'}")


if __name__ == "__main__":
    main()
