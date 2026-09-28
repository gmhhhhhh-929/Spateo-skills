"""Geometry diagnostics, not ground-truth anatomical or lineage validation."""
import numpy as np
from scipy.spatial import cKDTree


def pair_qc(source, target, source_labels, target_labels):
    x, y = np.asarray(source, float), np.asarray(target, float)
    a, b = np.asarray(source_labels, str), np.asarray(target_labels, str)
    da, ia = cKDTree(y).query(x)
    db, ib = cKDTree(x).query(y)
    result = {'symmetric_nn_mean': float((da.mean()+db.mean())/2),
              'nearest_label_agreement_source': float(np.mean(a == b[ia])),
              'nearest_label_agreement_target': float(np.mean(b == a[ib])),
              'interpretation': 'Descriptive geometry only; tissue growth, density and sampling affect these values.',
              'annotations': {}}
    for label in sorted(set(a) & set(b)):
        u, v = x[a == label], y[b == label]
        du, dv = cKDTree(v).query(u)[0], cKDTree(u).query(v)[0]
        result['annotations'][label] = {
            'source_n': len(u), 'target_n': len(v),
            'symmetric_nn_mean': float((du.mean()+dv.mean())/2),
            'source_to_target_p95': float(np.quantile(du, .95)),
            'target_to_source_p95': float(np.quantile(dv, .95)),
            'centroid_distance': float(np.linalg.norm(u.mean(0)-v.mean(0))),
            'source_extent': np.ptp(u, axis=0).tolist(),
            'target_extent': np.ptp(v, axis=0).tolist(),
        }
    # Only report this anatomical clue when both observed labels are present.
    # Z is a specimen axis, not a universal anatomical direction; review the
    # original reconstruction before treating its sign as an acceptance gate.
    if {'dorsal epidermal', 'ventral epidermal'} <= set(a) & set(b):
        delta = [float(points[labels == 'dorsal epidermal', 2].mean()
                       - points[labels == 'ventral epidermal', 2].mean())
                 for points, labels in ((x, a), (y, b))]
        result['dorsoventral_z'] = {
            'source_dorsal_minus_ventral': delta[0],
            'target_dorsal_minus_ventral': delta[1],
            'opposite_signs': bool(delta[0] * delta[1] < 0),
            'interpretation': 'Mean Z ordering is an orientation clue, not a universal anatomical axis or proof by itself.'
        }
    return result
