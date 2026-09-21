"""Verify the relocated viewer resolves its QC runtime without --qc-skill."""

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse
from spateo.preprocessing.slice_quality import (
    scan_h5ad_series,
    write_slice_quality_outputs,
)


def test_nested_viewer_default_runtime_and_complete_report(tmp_path):
    root = Path(__file__).resolve().parents[2]
    rng = np.random.default_rng(18)
    a = AnnData(
        sparse.csr_matrix(rng.poisson(3, (300, 12))),
        obs=pd.DataFrame(
            {"slice": np.repeat(["001", "002", "003"], 100)},
            index=[f"c{i}" for i in range(300)],
        ),
    )
    a.obsm["spatial"] = rng.normal(size=(300, 2))
    source = tmp_path / "input.h5ad"
    a.write_h5ad(source)
    result = scan_h5ad_series([source], slice_key="slice", spatial_key="spatial")
    run = tmp_path / "scan"
    write_slice_quality_outputs(result, run, write_display_payload=True)
    qc = root / "skills/spatial-slice-quality-qc"
    output = tmp_path / "viewer"
    proc = subprocess.run(
        [
            sys.executable,
            str(qc / "subskills/spatial-slice-quality-viewer/scripts/build_viewer.py"),
            "--input-dir",
            str(run),
            "--output-dir",
            str(output),
            "--policy",
            str(qc / "policies/joint_review_v2.json"),
            "--application-scope",
            "experimental_policy",
        ],
        text=True,
        capture_output=True,
    )
    assert proc.returncode == 0, proc.stderr
    assert (output / "index.html").is_file()
    ledger = json.loads((output / "viewer_workflow.json").read_text())
    assert Path(ledger["qc_skill"]) == qc
    assert ledger["source_h5ad_modified"] is False
    audit = pd.read_csv(
        output / "slice_quality_binary_audit.csv", dtype={"slice_id": str}
    )
    assert audit.slice_id.tolist() == ["001", "002", "003"]
    assert audit.final_call.isin(["keep", "exclude"]).all()
