"""Behavioral CLI checks against real Spateo readers; tiny synthetic inputs only."""
import json
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pytest
import spateo as st
from scipy import sparse
from scipy.io import mmwrite

import spateo_io as cli


def gem(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("geneID\tx\ty\tMIDCounts\ng1\t1\t2\t3\ng2\t1\t2\t4\ng1\t5\t7\t2\n")
    return path


def test_discovery_never_hashes_the_full_source(tmp_path, monkeypatch, capsys):
    source = gem(tmp_path / "sample.gem")
    def forbid_hash(_):
        raise AssertionError("Discovery scanned a full input for hashing")
    monkeypatch.setattr(cli, "sha256", forbid_hash)
    assert cli.main(["discover", str(source)]) == 0
    report = json.loads(capsys.readouterr().out)
    assert report["status"] == "pending"
    assert all(entry["materialization_attempts"] == 0 for entry in report["datasets"].values())


def test_lazy_exports_only_selected_entry_and_reports_remaining_scope(tmp_path, capsys):
    gem(tmp_path / "inputs/a.gem")
    gem(tmp_path / "inputs/b.gem")
    inspected = st.io.read_spatial(tmp_path / "inputs", load=False, load_images=False)
    assert len(inspected.datasets) == 2
    key = next(iter(inspected.datasets))
    out = tmp_path / "output"
    assert cli.main(["read", str(tmp_path / "inputs"), "--lazy", "--dataset-key", key,
                     "--output-dir", str(out), "--matrix-semantics", "counts"]) == 2
    report = json.loads(capsys.readouterr().out)
    assert report["status"] == "partial"
    assert report["export_status"] == "complete"
    assert list(report["exports"]) == [key]
    assert report["datasets"][key]["materialization_attempts"] == 1
    other = next(k for k in report["datasets"] if k != key)
    assert report["datasets"][other]["status"] == "deferred"
    assert report["datasets"][other]["materialization_attempts"] == 0
    loaded = ad.read_h5ad(report["exports"][key]["output"])
    assert loaded.shape == (2, 2) and int(loaded.X.sum()) == 9
    assert json.loads((out / "read_report.json").read_text()) == report


def test_lazy_unknown_selection_never_materializes(tmp_path, capsys):
    source = gem(tmp_path / "sample.gem")
    assert cli.main(["read", str(source), "--lazy", "--dataset-key", "not-present",
                     "--output-dir", str(tmp_path / "output")]) == 2
    report = json.loads(capsys.readouterr().out)
    assert report["unknown_selected_keys"] == ["not-present"]
    assert report["cli"]["input_file_sha256"] is None
    assert not report["exports"]
    assert all(entry["materialization_attempts"] == 0 for entry in report["datasets"].values())


def test_output_collision_rejected_before_read(tmp_path, monkeypatch):
    out = tmp_path / "output"
    out.mkdir()
    marker = out / "keep.txt"
    marker.write_text("existing result")
    def forbid_read(*args, **kwargs):
        raise AssertionError("Reader ran before checking output collision")
    monkeypatch.setattr(st.io, "read_spatial", forbid_read)
    with pytest.raises(FileExistsError):
        cli.main(["read", str(tmp_path / "not-even-present"), "--output-dir", str(out)])
    assert marker.read_text() == "existing result"


def test_missing_source_has_persisted_recovery(tmp_path, capsys):
    out = tmp_path / "output"
    assert cli.main(["read", str(tmp_path / "absent"), "--output-dir", str(out)]) == 2
    report = json.loads(capsys.readouterr().out)
    assert report["status"] == "failed" and not report["exports"]
    assert any(d.get("recovery") for d in report["diagnostics"])
    assert json.loads((out / "read_report.json").read_text()) == report


def test_invalid_options_preserve_failure_report(tmp_path):
    source = gem(tmp_path / "sample.gem")
    out = tmp_path / "output"
    with pytest.raises(ValueError):
        cli.main(["read", str(source), "--max-memory-bytes", "-1", "--output-dir", str(out)])
    report = json.loads((out / "read_report.json").read_text())
    assert report["status"] == "failed" and report["exception"]["type"] == "ValueError"
    assert not list(out.glob("*.h5ad"))


def test_manifest_serialization_failure_leaves_no_apparent_output(tmp_path):
    output = tmp_path / "result.h5ad"
    data = ad.AnnData(np.array([[1, 2]]))
    with pytest.raises(TypeError):
        cli.save_checked(data, output, {"status": "pass"}, {"not_serializable": object()})
    assert not list(tmp_path.iterdir())


def test_sidecar_publication_collision_rolls_back_only_own_file(tmp_path, monkeypatch):
    output = tmp_path / "result.h5ad"
    sidecar = output.with_name(output.name + ".manifest.json")
    original_link = cli.os.link
    def collide(source, dest):
        if Path(dest) == sidecar:
            sidecar.write_text("concurrent owner")
        return original_link(source, dest)
    monkeypatch.setattr(cli.os, "link", collide)
    with pytest.raises(FileExistsError):
        cli.save_checked(ad.AnnData(np.array([[1, 2]])), output, {}, {})
    assert not output.exists()
    assert sidecar.read_text() == "concurrent owner"
    assert list(tmp_path.iterdir()) == [sidecar]


def domestic(root, platform):
    root.mkdir(parents=True)
    matrix = sparse.csc_matrix([[1, 3], [2, 4]])
    if platform == "singleron":
        with h5py.File(root / "filtered_feature_bc_matrix.h5", "w") as handle:
            handle.attrs["chemistry_description"] = "Spatial3"
            group = handle.create_group("matrix")
            for name, value in [("data", matrix.data), ("indices", matrix.indices),
                                ("indptr", matrix.indptr), ("shape", matrix.shape),
                                ("barcodes", np.array(["c1", "c2"], dtype="S"))]:
                group.create_dataset(name, data=value)
            features = group.create_group("features")
            for name, value in [("id", ["g1", "g2"]), ("name", ["same", "same"]),
                                ("feature_type", ["Gene Expression"] * 2)]:
                features.create_dataset(name, data=np.array(value, dtype="S"))
        (root / "spatial").mkdir()
        (root / "spatial/positions_list.csv").write_text("c2,1,0,1,40,30\nc1,1,0,0,20,10\n")
    else:
        mmwrite(root / "matrix.mtx", matrix)
        (root / "features.tsv").write_text("g1\tsame\tGene Expression\ng2\tsame\tGene Expression\n")
        (root / "barcodes.tsv").write_text("c1\nc2\n")
        name, content = {
            "seekspace": ("cell_locations.tsv", "Cell_Barcode\tX\tY\nc2\t30\t40\nc1\t10\t20\n"),
            "bmkmanu": ("barcodes_pos.tsv", "c2\t30\t40\nc1\t10\t20\n"),
            "salus": ("spatial.txt", "c2 30 40\nc1 10 20\n"),
        }[platform]
        (root / name).write_text(content)
    return root


@pytest.mark.parametrize("platform", ["seekspace", "bmkmanu", "salus", "singleron"])
def test_domestic_auto_lazy_and_explicit_cli_roundtrip(tmp_path, capsys, platform):
    source = domestic(tmp_path / "sample", platform)
    out = tmp_path / "auto"
    assert cli.main(["read", str(source), "--lazy", "--output-dir", str(out),
                     "--matrix-semantics", "counts"]) == 0
    report = json.loads(capsys.readouterr().out)
    entry = next(iter(report["datasets"].values()))
    assert entry["technology"] == platform and entry["materialization_attempts"] == 1
    record = next(iter(report["exports"].values()))
    auto = ad.read_h5ad(record["output"])
    direct_path = tmp_path / "explicit.h5ad"
    assert cli.main(["convert", str(source), str(direct_path), "--reader", "read_" + platform,
                     "--matrix-semantics", "counts", "--require-spatial"]) == 0
    direct = ad.read_h5ad(direct_path)
    assert auto.var_names.tolist() == direct.var_names.tolist() == ["g1", "g2"]
    assert auto.var["gene_name"].tolist() == ["same", "same"]
    np.testing.assert_array_equal(auto.X.toarray(), [[1, 2], [3, 4]])
    np.testing.assert_array_equal(auto.obsm["spatial"], [[10, 20], [30, 40]])
    cli.assert_equal(auto.X, direct.X)
    cli.assert_equal(auto.obsm["spatial"], direct.obsm["spatial"])
