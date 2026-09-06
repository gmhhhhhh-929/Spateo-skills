#!/usr/bin/env python3
"""Small real-reader tests, using an existing scientific Python environment."""
import argparse
import contextlib
import importlib.metadata
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile


COMMIT = "d6aa68addc475dd0b56f69cebe7823b1f79933a9"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", required=True, type=Path)
    args = parser.parse_args()
    source = args.source_root.resolve()
    commit = subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip()
    if commit != COMMIT:
        raise ValueError("Source commit differs from this skill: " + commit)
    sys.dont_write_bytecode = True
    sys.path.insert(0, str(source))
    with tempfile.TemporaryDirectory(prefix="spateo-io-smoke-") as temporary:
        root = Path(temporary)
        os.environ.setdefault("MPLBACKEND", "Agg")
        os.environ.setdefault("MPLCONFIGDIR", str(root / "mpl"))
        import anndata as ad
        import h5py
        import numpy as np
        import pandas as pd
        import scipy.io
        from scipy.sparse import csc_matrix
        from PIL import Image
        import spateo.io as io
        from validate_anndata import validate
        if Path(io.__file__).resolve() != source / "spateo/io/__init__.py":
            raise RuntimeError("Not importing the selected source checkout")
        checks = []
        def checked(name):
            checks.append(name)
            print("PASS " + name, file=sys.stderr, flush=True)

        counts = np.array([[1, 2, 0, 4], [0, 5, 6, 0]], dtype=np.int32)  # genes × cells
        ids = ["c1", "c2", "c3", "c4"]
        xy = np.array([[10., 40.], [20., 50.], [30., 60.], [40., 70.]])
        def write_h5(path):
            matrix = csc_matrix(counts)
            with h5py.File(path, "w") as h:
                m = h.create_group("matrix")
                for key, value in (("data", matrix.data), ("indices", matrix.indices),
                                   ("indptr", matrix.indptr), ("shape", matrix.shape)):
                    m.create_dataset(key, data=value)
                m.create_dataset("barcodes", data=np.array(ids, dtype="S"))
                f = m.create_group("features")
                f.create_dataset("id", data=np.array(["g1", "g2"], dtype="S"))
                f.create_dataset("name", data=np.array(["Gene1", "Gene2"], dtype="S"))
                f.create_dataset("feature_type", data=np.array(["Gene Expression"] * 2, dtype="S"))
        with contextlib.redirect_stdout(sys.stderr):
            vis = root / "visium"
            (vis / "spatial").mkdir(parents=True)
            write_h5(vis / "filtered_feature_bc_matrix.h5")
            matrix = io.read_10x_h5(vis / "filtered_feature_bc_matrix.h5")
            np.testing.assert_array_equal(matrix.X.toarray(), counts.T)
            assert matrix.obs_names.tolist() == ids and "spatial" not in matrix.obsm
            checked("10x_h5_gene_cell_orientation_no_invented_spatial")
            positions = pd.DataFrame({"barcode": ids, "in_tissue": [1, 1, 1, 0],
                                      "array_row": [1, 2, 3, 4], "array_col": [5, 6, 7, 8],
                                      "pxl_row_in_fullres": xy[:, 1], "pxl_col_in_fullres": xy[:, 0]})
            positions.iloc[::-1].to_csv(vis / "spatial/tissue_positions.csv", index=False)
            (vis / "spatial/scalefactors_json.json").write_text(json.dumps({"tissue_hires_scalef": 0.5}))
            pixels = np.arange(36, dtype=np.uint8).reshape(6, 6)
            Image.fromarray(pixels).save(vis / "spatial/tissue_hires_image.png")
            a, match = io.read_auto_spatial(vis, technology="visium", load_images=True, return_match=True)
            np.testing.assert_array_equal(a.obsm["spatial"], xy)
            assert a.n_obs == 4 and int(a.obs.loc["c4", "in_tissue"]) == 0
            slot = next(iter(a.uns["spatial"].values()))
            np.testing.assert_array_equal(slot["images"]["hires"], pixels)
            assert slot["scalefactors"]["tissue_hires_scalef"] == 0.5
            assert a.uns["spateo_io"]["technology"] == match.technology == "visium"
            checked("visium_id_join_xy_axis_images_scalefactors_auto_provenance")
            cached = root / "roundtrip.h5ad"
            a.write_h5ad(cached)
            reread = ad.read_h5ad(cached)
            np.testing.assert_array_equal(reread.obsm["spatial"], xy)
            np.testing.assert_array_equal(reread.X.toarray(), counts.T)
            np.testing.assert_array_equal(next(iter(reread.uns["spatial"].values()))["images"]["hires"], pixels)
            assert reread.obs_names.tolist() == ids
            checked("h5ad_roundtrip_ids_counts_xy_image")
            missing = a.copy()
            missing.obsm["spatial"][0, 0] = np.nan
            assert validate(missing, coordinate_unit="pixel")["status"] == "fail"
            checked("nonfinite_coordinates_rejected")
            xyz = a.copy()
            xyz.obsm["spatial"] = np.column_stack([xy, [1., 2., 3., 4.]])
            assert validate(xyz, coordinate_unit="pixel")["status"] == "pass"
            assert validate(xyz, coordinate_unit="pixel", require_2d=True)["status"] == "fail"
            checked("native_xyz_preserved_but_2d_handoff_rejects_three_columns")
            a.obs["anno"] = ["A", "B", "A", "B"]
            a.obsm["X_pca"] = np.eye(2)[[0, 1, 0, 1]]
            assert validate(a, coordinate_unit="pixel", alignment_mode="expression-pca")["status"] == "fail"
            assert validate(a, coordinate_unit="pixel", alignment_mode="annotation-onehot", annotation_key="anno", categories=["A", "B"])["status"] == "pass"
            assert validate(a, coordinate_unit="pixel", alignment_mode="annotation-onehot", annotation_key="anno", categories=["B", "A"])["status"] == "fail"
            a.obsm["X_pca"] = np.array([[0., 0.], [1., -2.], [2., 3.], [-1., 4.]])
            assert validate(a, coordinate_unit="pixel", alignment_mode="expression-pca")["status"] == "fail"
            a.obsm["X_pca"] = np.ones((4, 30))
            assert validate(a, coordinate_unit="pixel", alignment_mode="spatial-only")["status"] == "pass"
            a.obsm["X_pca"][:] = 0
            assert validate(a, coordinate_unit="pixel", alignment_mode="spatial-only")["status"] == "fail"
            checked("onehot_expression_confusion_category_order_and_zero_cosine_rejected")
            atera = root / "atera"
            atera.mkdir()
            write_h5(atera / "cell_feature_matrix.h5")
            pd.DataFrame({"cell_id": ids[::-1], "x_centroid": xy[::-1, 0], "y_centroid": xy[::-1, 1]}).to_csv(atera / "cells.csv", index=False)
            (atera / "experiment.xenium").write_text(json.dumps({"run_name": "Atera WTA preview", "pixel_size": .25}))
            aa, am = io.read_auto_spatial(atera, load_image=False, load_boundaries=False,
                                        load_nucleus_boundaries=False, load_cell_groups=False, return_match=True)
            assert am.technology == "atera" and aa.uns["spateo_io"]["format_status"] == "preview-xenium-v4"
            np.testing.assert_array_equal(aa.obsm["spatial"], xy)
            assert next(iter(aa.uns["spatial"].values()))["scalefactors"]["tissue_hires_scalef"] == 4.
            checked("atera_preview_detection_id_order_units_not_rescaled")
            gem = root / "reads.gem"
            gem.write_text("geneID\tx\ty\tMIDCounts\tcell_id\nG1\t1\t2\t3\t1\nG2\t1\t2\t4\t1\nG1\t11\t12\t5\t2\n")
            grid = io.read_bgi_agg(gem)
            assert grid.uns["__type"] == "AGG" and int(grid.X.sum()) == 12
            assert validate(grid)["status"] == "fail"
            bins = io.read_bgi(gem, binsize=10)
            assert bins.shape == (2, 2) and int(bins.X.sum()) == 12
            np.testing.assert_array_equal(bins.obsm["spatial"], [[5, 5], [15, 15]])
            labels = io.read_auto_spatial(gem, technology="bgi", binsize=None, label_column="cell_id", add_props=False)
            assert labels.shape == (2, 2) and "spatial" not in labels.obsm
            try:
                io.read_bgi(gem, binsize=1, label_column="cell_id")
            except IOError:
                pass
            else:
                raise AssertionError("BGI accepted conflicting observation modes")
            checked("bgi_agg_vs_gene_matrix_counts_centroids_exclusive_modes")
            seq = root / "seqscope"
            seq.mkdir()
            (seq / "barcodes.tsv").write_text("b1\nb2\nb3\n")
            (seq / "features.tsv").write_text("g1\tG1\tGene Expression\ng2\tG2\tGene Expression\n")
            scipy.io.mmwrite(seq / "matrix.mtx", np.array([[1, 0, 2], [0, 3, 1]]))
            pos = root / "positions.txt"
            pos.write_text("b3 1 1 12 12\nb1 1 1 1 1\nb2 1 1 2 2\n")
            single = io.read_seqscope(seq, pos, binsize=None)
            grouped = io.read_seqscope(seq, pos, binsize=10, add_props=False)
            np.testing.assert_array_equal(single.obsm["spatial"], [[1, 1], [2, 2], [12, 12]])
            assert grouped.shape == (2, 2)
            np.testing.assert_array_equal(np.asarray(grouped.X.sum(axis=0)).ravel(), [3, 4])
            checked("seqscope_barcode_order_and_bin_mass_conservation")
            hd = root / "hd"
            square = hd / "binned_outputs/square_008um"
            (square / "spatial").mkdir(parents=True)
            (square / "filtered_feature_bc_matrix.h5").touch()
            (square / "spatial/tissue_positions.parquet").touch()
            seg = hd / "segmented_outputs"
            seg.mkdir()
            (seg / "filtered_feature_cell_matrix.h5").touch()
            (seg / "graphclust_annotated_cell_segmentations.geojson").touch()
            try:
                io.detect_spatial_technology(hd)
            except ValueError:
                pass
            else:
                raise AssertionError("Ambiguous HD bin/cellseg unexpectedly accepted")
            assert io.detect_spatial_technology(hd, technology="visium_hd_bin").path == square.resolve()
            checked("auto_hd_ambiguity_explicit_mode_and_subdirectory")
            (root / "table.csv").write_text("id,x\na,1\n")
            table = io.read_csv(filepath_or_buffer=str(root / "table.csv"))
            assert isinstance(table, pd.DataFrame)
            checked("general_csv_returns_dataframe")
            scripts = Path(__file__).resolve().parent
            env = dict(os.environ, PYTHONPATH=str(source), PYTHONDONTWRITEBYTECODE="1")
            def invoke(script, arguments, success=True):
                result = subprocess.run([sys.executable, str(scripts / script)] + arguments,
                                        env=env, text=True, capture_output=True)
                if success and result.returncode != 0:
                    raise AssertionError(result.stderr)
                if not success and result.returncode == 0:
                    raise AssertionError("Expected CLI failure: " + str(arguments))
                return result
            detected = invoke("spateo_io.py", ["detect", str(vis), "--technology", "visium"])
            assert json.loads(detected.stdout)[0]["technology"] == "visium"
            destination = root / "cli_output.h5ad"
            options = root / "reader_kwargs.json"
            options.write_text(json.dumps({"load_images": False}))
            command = ["convert", str(vis), str(destination), "--reader", "read_visium",
                       "--kwargs-json", str(options), "--matrix-semantics", "counts",
                       "--coordinate-unit", "pixel", "--require-spatial"]
            converted = invoke("spateo_io.py", command)
            assert json.loads(converted.stdout)["status"] == "saved"
            frozen = destination.read_bytes()
            duplicate = invoke("spateo_io.py", command, success=False)
            assert "Refusing to overwrite" in duplicate.stderr and destination.read_bytes() == frozen
            verified = invoke("validate_anndata.py", [str(destination), "--matrix-semantics", "counts",
                               "--coordinate-unit", "pixel", "--require-spatial"])
            assert json.loads(verified.stdout)["status"] == "pass"
            checked("cli_detect_convert_roundtrip_readonly_validate_and_no_overwrite")
            xyz_input = root / "native_xyz.h5ad"
            xyz.write_h5ad(xyz_input)
            xyz_output = root / "native_xyz_cli.h5ad"
            invoke("spateo_io.py", ["convert", str(xyz_input), str(xyz_output), "--reader", "read_h5ad"])
            np.testing.assert_array_equal(ad.read_h5ad(xyz_output).obsm["spatial"], xyz.obsm["spatial"])
            checked("cli_native_xyz_roundtrip_without_truncation")
        versions = {}
        for name in ["anndata", "numpy", "scipy", "pandas", "h5py", "matplotlib", "Pillow"]:
            versions[name] = importlib.metadata.version(name)
        print(json.dumps({"status": "pass", "source_commit": commit, "python": sys.version,
                          "supported_install_python": ">=3.10,<3.13", "packages": versions,
                          "tests": checks, "n_tests": len(checks), "reader_mocks_used": False}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
