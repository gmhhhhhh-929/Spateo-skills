import hashlib
import json
from pathlib import Path

import numpy as np
import pyvista as pv
import pytest
import spateo as st

import build_meshes as mesh_builder


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_point_cloud(path: Path, *, planar: bool = False) -> pv.PolyData:
    axes = [
        np.linspace(-4.0, 4.0, 13),
        np.linspace(-7.0, 7.0, 19),
        np.asarray([0.0]) if planar else np.linspace(-3.0, 3.0, 11),
    ]
    points = np.asarray(np.meshgrid(*axes, indexing="ij")).reshape(3, -1).T
    if planar:
        keep = (points[:, 0] / 4) ** 2 + (points[:, 1] / 7) ** 2 <= 1
    else:
        keep = (points[:, 0] / 4) ** 2 + (points[:, 1] / 7) ** 2 + (
            points[:, 2] / 3
        ) ** 2 <= 1
    points = points[keep]
    labels = np.where(points[:, 1] < -1.0, "anterior", "posterior")
    labels[np.abs(points[:, 1]) <= 1.0] = "center"
    model = pv.PolyData(points)
    model.point_data["obs_index"] = np.asarray(
        [f"cell-{index:04d}" for index in range(len(points))]
    )
    model.point_data["region"] = labels
    rgba = np.ones((len(points), 4), dtype=float)
    rgba[labels == "anterior", :3] = (0.85, 0.2, 0.35)
    rgba[labels == "center", :3] = (0.9, 0.7, 0.2)
    rgba[labels == "posterior", :3] = (0.2, 0.55, 0.8)
    model.point_data["region_rgba"] = rgba
    st.tdr.save_model(model, str(path), binary=True)
    return model


def run_builder(*arguments: str) -> dict:
    args = mesh_builder.build_parser().parse_args(list(arguments))
    return mesh_builder.build_meshes(args)


def test_density_body_and_multi_value_selection_roundtrip(tmp_path):
    input_path = tmp_path / "points.vtk"
    source = write_point_cloud(input_path)
    before_hash = sha256(input_path)
    output = tmp_path / "density"

    manifest = run_builder(
        "--input-vtk",
        str(input_path),
        "--output-dir",
        str(output),
        "--name",
        "specimen",
        "--label-key",
        "region",
        "--selection",
        "front=anterior|center",
        "--selection",
        "back=posterior",
        "--target-long-axis",
        "64",
        "--sigma",
        "1.2",
        "--smooth",
        "8",
        "--body-coverage-target",
        "0.9",
        "--selection-coverage-target",
        "0.85",
        "--min-component-voxels",
        "4",
        "--preview-max-points",
        "300",
    )

    assert sha256(input_path) == before_hash
    assert manifest["status"] == "pass"
    assert manifest["validated_against_revision"].startswith("615644f")
    assert (
        manifest["environment"]["runtime_mesh_methods_matches_validated_source"] is True
    )
    assert manifest["method"] == "density"
    assert [row["name"] for row in manifest["meshes"]] == [
        "body",
        "front",
        "back",
    ]
    assert (output / "specimen.preview.png").stat().st_size > 0
    stored = json.loads((output / "specimen.manifest.json").read_text())
    assert stored["shared_grid"]["n_voxels"] > 0
    assert stored["parameters"]["active_groups"][1] == "density"
    assert stored["parameters"]["density"]["target_long_axis"] == 64
    assert stored["parameters"]["density"]["voxel_size_mode"] == "target_long_axis"
    assert stored["parameters"]["preview"]["max_points"] == 300
    assert stored["meshes"][1]["selected_values"] == ["anterior", "center"]
    for row in stored["meshes"]:
        loaded = st.tdr.read_model(row["vtk"])
        assert isinstance(loaded, pv.PolyData)
        assert loaded.n_points > 0 and loaded.n_faces_strict > 0
        assert loaded.n_lines == loaded.n_verts == loaded.n_strips == 0
        assert loaded.n_open_edges == 0
        assert {"mesh_label", "mesh_label_rgba"}.issubset(loaded.cell_data)
        assert row["roundtrip"]["coordinates_allclose"] is True
        assert row["roundtrip"]["face_connectivity_equal"] is True
        assert row["roundtrip"]["required_cell_data_values_equal"] is True
        assert row["reconstruction"]["level_selection"]["target_met"] is True
    assert source.n_points == manifest["input"]["n_points"]


def test_native_spateo_core_exposes_geometry_parameters(tmp_path):
    input_path = tmp_path / "points.vtk"
    write_point_cloud(input_path)
    output = tmp_path / "native"

    manifest = run_builder(
        "--input-vtk",
        str(input_path),
        "--output-dir",
        str(output),
        "--name",
        "native",
        "--method",
        "spateo-marching-cube",
        "--levelset",
        "0.5",
        "--mc-scale-factor",
        "0.9",
        "--dist-sample-num",
        "30",
        "--smooth",
        "5",
        "--skip-preview",
    )

    reconstruction = manifest["meshes"][0]["reconstruction"]
    assert reconstruction["algorithm"] == "spateo_marching_cube_mesh_core"
    assert reconstruction["levelset"] == 0.5
    assert reconstruction["mc_scale_factor"] == 0.9
    assert reconstruction["smooth_iterations"] == 5
    assert reconstruction["random_seed"] == 0
    assert manifest["parameters"]["active_groups"][1] == "spateo-marching-cube"
    assert manifest["parameters"]["spateo-marching-cube"]["random_seed"] == 0
    loaded = st.tdr.read_model(manifest["meshes"][0]["vtk"])
    assert loaded.n_faces_strict > 0
    assert loaded.n_lines == 0

    repeated = run_builder(
        "--input-vtk",
        str(input_path),
        "--output-dir",
        str(tmp_path / "native-repeated"),
        "--name",
        "native",
        "--method",
        "spateo-marching-cube",
        "--levelset",
        "0.5",
        "--mc-scale-factor",
        "0.9",
        "--dist-sample-num",
        "30",
        "--smooth",
        "5",
        "--skip-preview",
    )
    repeated_mesh = st.tdr.read_model(repeated["meshes"][0]["vtk"])
    assert np.array_equal(loaded.faces, repeated_mesh.faces)
    assert np.allclose(loaded.points, repeated_mesh.points)


def test_repair_preserves_disconnected_components():
    parts = []
    for center in [(-3, 0, 0), (3, 0, 0)]:
        sphere = pv.Sphere(
            theta_resolution=24, phi_resolution=24, center=center
        ).triangulate()
        opened = (
            sphere.extract_cells(np.arange(sphere.n_cells - 2))
            .extract_surface()
            .triangulate()
            .clean()
        )
        parts.append(opened)
    mesh = parts[0].merge(parts[1], merge_points=False).extract_surface().triangulate()

    repaired, summary = mesh_builder._repair_open_surface(mesh, enabled=True)

    assert summary["open_edges_before"] > 0
    assert repaired.n_open_edges == 0
    assert summary["topology_before"]["components"] == 2
    assert summary["topology_after"]["components"] == 2
    assert summary["components_preserved"] is True


def test_failure_discards_staging_outputs(tmp_path, monkeypatch):
    input_path = tmp_path / "points.vtk"
    write_point_cloud(input_path)
    output = tmp_path / "atomic-output"
    real_density_mesh = mesh_builder._density_mesh
    calls = 0

    def fail_second_target(*args, **kwargs):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise RuntimeError("synthetic second-target failure")
        return real_density_mesh(*args, **kwargs)

    monkeypatch.setattr(mesh_builder, "_density_mesh", fail_second_target)
    with pytest.raises(RuntimeError, match="second-target"):
        run_builder(
            "--input-vtk",
            str(input_path),
            "--output-dir",
            str(output),
            "--label-key",
            "region",
            "--selection",
            "front=anterior|center",
            "--target-long-axis",
            "48",
            "--sigma",
            "1.2",
            "--smooth",
            "0",
            "--body-coverage-target",
            "0.85",
            "--selection-coverage-target",
            "0.80",
            "--min-component-voxels",
            "4",
            "--skip-preview",
        )

    assert not output.exists()
    assert not list(tmp_path.glob(".atomic-output.staging-*"))


def test_automatic_density_level_must_meet_coverage(tmp_path):
    input_path = tmp_path / "points.vtk"
    write_point_cloud(input_path)
    output = tmp_path / "unqualified"

    with pytest.raises(ValueError, match="did not reach its density coverage target"):
        run_builder(
            "--input-vtk",
            str(input_path),
            "--output-dir",
            str(output),
            "--target-long-axis",
            "48",
            "--sigma",
            "1.2",
            "--min-component-voxels",
            "1000000",
            "--skip-preview",
        )

    assert not output.exists()


def test_missing_selection_and_planar_point_cloud_fail(tmp_path):
    input_path = tmp_path / "points.vtk"
    write_point_cloud(input_path)
    with pytest.raises(ValueError, match="absent"):
        run_builder(
            "--input-vtk",
            str(input_path),
            "--output-dir",
            str(tmp_path / "missing"),
            "--label-key",
            "region",
            "--selection",
            "bad=not-present",
            "--skip-preview",
        )

    planar_path = tmp_path / "planar.vtk"
    write_point_cloud(planar_path, planar=True)
    with pytest.raises(ValueError, match="coordinate rank"):
        run_builder(
            "--input-vtk",
            str(planar_path),
            "--output-dir",
            str(tmp_path / "planar-output"),
            "--skip-preview",
        )


def test_output_overwrite_refusal(tmp_path):
    input_path = tmp_path / "points.vtk"
    write_point_cloud(input_path)
    output = tmp_path / "result"
    common = (
        "--input-vtk",
        str(input_path),
        "--output-dir",
        str(output),
        "--target-long-axis",
        "48",
        "--sigma",
        "1.2",
        "--smooth",
        "0",
        "--body-coverage-target",
        "0.85",
        "--min-component-voxels",
        "4",
        "--skip-preview",
    )
    run_builder(*common)
    with pytest.raises(FileExistsError, match="non-empty"):
        run_builder(*common)


def test_case_only_mesh_name_collisions_fail_before_writing(tmp_path):
    input_path = tmp_path / "points.vtk"
    write_point_cloud(input_path)
    output = tmp_path / "case-collision"

    with pytest.raises(ValueError, match="ignoring letter case"):
        run_builder(
            "--input-vtk",
            str(input_path),
            "--output-dir",
            str(output),
            "--label-key",
            "region",
            "--selection",
            "organ=anterior",
            "--selection",
            "Organ=posterior",
            "--skip-preview",
        )

    assert not output.exists()
