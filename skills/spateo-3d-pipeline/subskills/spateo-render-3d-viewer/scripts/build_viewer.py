#!/usr/bin/env python3
"""Import immutable VTK models into a portable, full-precision review page."""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import tempfile
from pathlib import Path

import numpy as np
import plotly
import pyvista as pv
import spateo as st
from plotly.offline import get_plotlyjs


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1048576), b''):
            h.update(chunk)
    return h.hexdigest()


def json_text(value):
    return json.dumps(value, ensure_ascii=False, allow_nan=False, separators=(',', ':'))


def id_hash(ids):
    return hashlib.sha256(json_text(np.asarray(ids).astype(str).tolist()).encode()).hexdigest()


def resolve(base, value):
    path = Path(value).expanduser()
    return (base / path).resolve() if not path.is_absolute() else path.resolve()


def read_model(path, sources, *, cells=False):
    if path.suffix.lower() not in ('.vtk', '.vtp'):
        raise ValueError(f'Expected .vtk or .vtp: {path}')
    sources[str(path)] = digest(path)
    model = st.tdr.read_model(str(path))
    if not isinstance(model, pv.PolyData) or not model.n_points:
        raise ValueError(f'Expected nonempty PolyData: {path}')
    if not np.isfinite(model.points).all():
        raise ValueError(f'Nonfinite coordinates: {path}')
    if cells:
        if 'obs_index' not in model.point_data:
            raise ValueError(f'Cells require obs_index: {path}')
        ids = np.asarray(model.point_data['obs_index']).astype(str)
        if ids.shape != (model.n_points,) or len(np.unique(ids)) != model.n_points:
            raise ValueError(f'Cells require unique scalar obs_index: {path}')
    else:
        # Do not include stray vertices/lines as surface primitives.
        model = model.triangulate()
        if not model.faces.size:
            raise ValueError(f'Expected polygon surface: {path}')
        model = pv.PolyData(model.points.copy(), model.faces.copy())
    return model


def unique_names(rows, context):
    if not isinstance(rows, list) or not rows:
        raise ValueError(f'{context} must be a nonempty list')
    names = [r.get('name') for r in rows]
    if any(not isinstance(n, str) or not n.strip() for n in names) or len(set(names)) != len(names):
        raise ValueError(f'{context} requires unique nonempty names')


def layers_from_dataset(spec, base, sources):
    if 'mesh_manifest' not in spec:
        return spec.get('layers')
    if 'layers' in spec:
        raise ValueError('Use layers or mesh_manifest, not both')
    path = resolve(base, spec['mesh_manifest'])
    sources[str(path)] = digest(path)
    rows = json.loads(path.read_text())['meshes']
    if not isinstance(rows, list):
        raise ValueError('mesh_manifest must be a build_meshes.py list-based manifest')
    result = []
    for row in rows:
        layer = dict(name=row['name'], role=row.get('role', 'tissue'),
                     mesh=str(resolve(path.parent, row['vtk'])), color=row['color'])
        if row.get('selected_values') and spec.get('point_cloud'):
            layer['values'] = row['selected_values']
        result.append(layer)
    return result


def select_cells(layer, spec, shared, base, sources):
    if 'cells' in layer:
        if 'values' in layer:
            raise ValueError('Use cells or values, not both')
        pc = read_model(resolve(base, layer['cells']), sources, cells=True)
        return pc.points, np.asarray(pc['obs_index']).astype(str)
    if 'values' not in layer and layer.get('role') != 'body':
        return None, None
    if shared is None:
        if 'values' in layer:
            raise ValueError('values requires a shared point_cloud')
        return None, None
    mask = np.ones(shared.n_points, dtype=bool)
    if 'values' in layer:
        values = layer['values']
        if not isinstance(values, list) or not values or any(not isinstance(v, str) for v in values):
            raise ValueError('values must be a nonempty list of exact category strings')
        key = spec.get('label_key')
        if key not in shared.point_data:
            raise ValueError(f'Missing point-cloud label_key: {key}')
        labels = np.asarray(shared.point_data[key]).astype(str)
        if labels.shape != (shared.n_points,):
            raise ValueError('label_key must contain one categorical value per point')
        missing = set(values) - set(labels)
        if missing:
            raise ValueError(f'Missing categories: {sorted(missing)}')
        mask = np.isin(labels, values)
    return np.asarray(shared.points)[mask], np.asarray(shared['obs_index']).astype(str)[mask]


def morphology(mesh, points):
    closed = bool(mesh.n_open_edges == 0 and mesh.is_manifold and mesh.volume > 0)
    inside_count = None
    selected = None
    if closed and points is not None:
        selected = pv.PolyData(points).select_enclosed_points(mesh, tolerance=1e-6, check_surface=True)
        hit = np.asarray(selected['SelectedPoints']).astype(bool)
        inside_count = int(hit.sum())
        # Density is computed explicitly below to handle zero enclosed cells.
    features = {k: float(v) for k, v in st.tdr.model_morphology(mesh).items()}
    if not closed:
        features['Volume'] = None
        features['V/SA_ratio'] = None
    elif points is not None:
        features['cell_density'] = inside_count / float(mesh.volume)
        features['Cells_inside'] = inside_count
    return features, closed


def build_payload(config, base, sources):
    datasets = config.get('datasets')
    unique_names(datasets, 'datasets')
    payload = {'title': str(config.get('title', '3D Reconstruction Review')), 'datasets': []}
    audits = []
    for spec in datasets:
        layers = layers_from_dataset(spec, base, sources)
        unique_names(layers, 'layers')
        budget = spec.get('point_budget', 20000)
        if type(budget) is not int or budget <= 0:
            raise ValueError('point_budget must be a positive integer')
        shared = read_model(resolve(base, spec['point_cloud']), sources, cells=True) if spec.get('point_cloud') else None
        dataset = {'name': spec['name'], 'units': str(spec.get('units', 'unspecified')), 'layers': []}
        audit = {'name': spec['name'], 'point_budget_per_layer': budget, 'layers': []}
        for layer in layers:
            color = layer.get('color', '#80B9DC')
            if not re.fullmatch(r'#[0-9a-fA-F]{6}', color):
                raise ValueError('color must be #RRGGBB')
            opacity = layer.get('opacity', .15 if layer.get('role') == 'body' else 1.0)
            if not isinstance(opacity, (float, int)) or not 0 <= opacity <= 1:
                raise ValueError('opacity must be between 0 and 1')
            visible = layer.get('visible', True)
            if not isinstance(visible, bool):
                raise ValueError('visible must be boolean')
            row = dict(name=layer['name'], role=layer.get('role', 'tissue'), color=color,
                       opacity=opacity, visible=visible, mesh=None, cells=None, features={})
            record = {'name': layer['name'], 'warnings': []}
            points, ids = select_cells(layer, spec, shared, base, sources)
            if layer.get('mesh'):
                path = resolve(base, layer['mesh'])
                mesh = read_model(path, sources)
                faces = mesh.faces.reshape(-1, 4)[:, 1:]
                row['mesh'] = {'points': mesh.points.tolist(), 'faces': faces.tolist()}
                row['features'], row['closed'] = morphology(mesh, points)
                record.update(mesh=str(path), vertices=mesh.n_points, faces=len(faces),
                              open_edges=mesh.n_open_edges, closed=row['closed'], features=row['features'])
                if not row['closed']:
                    record['warnings'].append('Open/nonmanifold/nonpositive-volume surface; volume metrics suppressed')
                if layer.get('role') == 'body':
                    record['warnings'].append('Body metrics describe the imported shell, not inferred anatomy')
            if points is not None:
                # Exact upper bound, stable ordered sample, independent of global RNG.
                index = np.linspace(0, len(points)-1, min(budget, len(points)), dtype=int)
                row['cells'] = {'points': points[index].tolist(), 'ids': ids[index].tolist(), 'source_count': len(points)}
                record.update(source_cells=len(points), display_cells=len(index), source_id_sha256=id_hash(ids),
                              display_id_sha256=id_hash(ids[index]), sampling='deterministic evenly spaced source indices')
            if row['mesh'] is None and row['cells'] is None:
                raise ValueError(f'Layer {layer["name"]} has neither mesh nor cells')
            dataset['layers'].append(row)
            audit['layers'].append(record)
        payload['datasets'].append(dataset)
        audits.append(audit)
    return payload, audits


def build(config_path, output_dir):
    config_path, output_dir = Path(config_path).resolve(), Path(output_dir).resolve()
    if output_dir.exists():
        raise FileExistsError(output_dir)
    sources = {str(config_path): digest(config_path)}
    config = json.loads(config_path.read_text())
    payload, datasets = build_payload(config, config_path.parent, sources)
    # Escape data rather than injecting untrusted labels into HTML/JS source.
    encoded = json_text(payload).replace('&', '\\u0026').replace('<', '\\u003c').replace('>', '\\u003e')
    template = Path(__file__).resolve().parents[1] / 'assets/viewer.html'
    html = template.read_text().replace('/*__PLOTLY__*/', get_plotlyjs()).replace('__MODEL_DATA__', encoded)
    if any(digest(Path(path)) != expected for path, expected in sources.items()):
        raise RuntimeError('An input changed during rendering')
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix='.viewer-', dir=output_dir.parent))
    try:
        (staging / 'index.html').write_text(html, encoding='utf-8')
        report = {'status': 'rendered_not_scientifically_validated', 'datasets': datasets,
                  'source_hashes': sources, 'source_files_unchanged': True,
                  'geometry_processing': 'triangulation only; no mesh simplification or coordinate rounding',
                  'versions': {'numpy': np.__version__, 'pyvista': pv.__version__, 'plotly': plotly.__version__,
                               'spateo': getattr(st, '__version__', 'unknown')},
                  'html_bytes': (staging / 'index.html').stat().st_size,
                  'html_sha256': digest(staging / 'index.html'), 'warnings': []}
        if report['html_bytes'] > 50 * 1024**2:
            report['warnings'].append('HTML exceeds 50 MiB; consider fewer layers or separate datasets')
        (staging / 'viewer_manifest.json').write_text(json.dumps(report, indent=2, allow_nan=False))
        if output_dir.exists():
            raise FileExistsError(output_dir)
        staging.rename(output_dir)
    finally:
        if staging.exists():
            shutil.rmtree(staging)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', required=True, type=Path)
    parser.add_argument('--output-dir', required=True, type=Path)
    args = parser.parse_args()
    report = build(args.config, args.output_dir)
    print(json.dumps({'html': str(args.output_dir.resolve() / 'index.html'), 'html_bytes': report['html_bytes']}))


if __name__ == '__main__':
    main()
