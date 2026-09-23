import numpy as np
import pyvista as pv
import pytest
from coverage_envelope import make_grid, coverage, triangle_containment, smooth_envelope, regional_candidate, repair_supported_gaps


def test_regional_audit_detects_sparse_missed_end_despite_high_global_coverage():
    rng=np.random.default_rng(42)
    core=rng.normal(0,.12,(1000,3))
    end=rng.normal([2,0,0],.03,(35,3))
    points=np.vstack([core,end]);mesh=pv.Sphere(radius=.8).triangulate()
    report,_=coverage(points,mesh,np.array([1.,0,0]),np.zeros(3))
    assert report['global']>.95
    assert report['min_axial']==0
    assert report['min_spatial_block']==0


def test_triangle_certificate_is_not_a_vertex_only_test():
    shell=pv.ParametricTorus().triangulate().compute_normals(auto_orient_normals=True,consistent_normals=True)
    tri=pv.PolyData(np.array([[1.,0,0],[-1,0,0],[0,1,0]]),np.array([3,0,1,2]))
    assert not triangle_containment(tri,shell)['passed']
    inner=pv.Sphere(radius=.2).triangulate()
    outer=pv.Sphere(radius=1.).triangulate().compute_normals(auto_orient_normals=True,consistent_normals=True)
    assert triangle_containment(inner,outer)['passed']


def test_smooth_union_does_not_double_count_overlapping_tissues():
    rng=np.random.default_rng(7);points=rng.uniform(-.3,.3,(300,3))
    tissue=pv.Sphere(radius=.7,theta_resolution=20,phi_resolution=20).triangulate()
    grid=make_grid(points,[tissue],long_axis=48)
    shell,qc=smooth_envelope(points,[tissue,tissue.copy()],grid)
    assert shell.n_open_edges==0
    assert all(c['passed'] for c in qc['attempts'][-1]['certificates'])
    assert qc['outside_tissue_voxels']==0
    assert 0<qc['union_over_body_voxel_fraction']<1
    assert qc['overlap_over_body_voxel_fraction']==qc['union_over_body_voxel_fraction']


def test_supported_regional_restoration_keeps_source_points_immutable():
    rng=np.random.default_rng(22)
    points=np.vstack([rng.normal(0,.1,(400,3)),rng.normal([.8,0,0],.04,(45,3))])
    original=points.copy();old=pv.Sphere(radius=.4).triangulate()
    grid=make_grid(points,[old],long_axis=64)
    mesh,_=regional_candidate(points,old,grid,np.array([1.,0,0]),np.zeros(3),quantile=.02,sigma=2,min_component_cells=10)
    mesh,_=repair_supported_gaps(points,mesh,grid,sigma=2,min_component_cells=10)
    before,_=coverage(points,old,np.array([1.,0,0]),np.zeros(3))
    after,_=coverage(points,mesh,np.array([1.,0,0]),np.zeros(3))
    assert after['global']>before['global']
    assert after['min_axial']>before['min_axial']
    assert mesh.n_open_edges==0
    assert np.array_equal(points,original)


def test_invalid_grid_fails_before_allocating_large_volume():
    with pytest.raises(ValueError):make_grid(np.zeros((2,2)))
    with pytest.raises(ValueError):make_grid(np.array([[0,0,0],[1,1,1]]),long_axis=1000,max_voxels=100)
