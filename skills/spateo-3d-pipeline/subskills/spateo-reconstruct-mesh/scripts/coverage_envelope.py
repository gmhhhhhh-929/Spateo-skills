"""Reusable regional density refinement and smooth enclosing volume fields.

All coordinates are physical XYZ in one frame. This module never edits source
points. The caller owns versioned paths, Spateo IO and biological review.
"""
from dataclasses import dataclass
import numpy as np
import pyvista as pv
from scipy import ndimage as ndi
from scipy.spatial import cKDTree
from skimage.measure import marching_cubes


@dataclass
class Grid:
    origin: np.ndarray
    shape: tuple
    spacing: float

    def image(self):
        return pv.ImageData(dimensions=self.shape, spacing=(self.spacing,)*3, origin=self.origin)

    def indices(self, points):
        ids=np.rint((np.asarray(points)-self.origin)/self.spacing).astype(int)
        if np.any(ids<0) or np.any(ids>=np.asarray(self.shape)):
            raise ValueError('Points outside shared grid')
        return ids


def make_grid(points, meshes=(), long_axis=220, padding=14, max_voxels=20000000):
    points=np.asarray(points,dtype=float)
    if points.ndim!=2 or points.shape[1]!=3 or not len(points) or not np.isfinite(points).all():
        raise ValueError('Expected nonempty finite XYZ')
    if long_axis<24 or padding<4:raise ValueError('Grid resolution/padding too small')
    bounds=np.vstack([points,*[m.points for m in meshes]])
    spacing=float(np.ptp(bounds,axis=0).max()/long_axis)
    if spacing<=0:raise ValueError('Zero coordinate extent')
    origin=bounds.min(0)-spacing*padding
    shape=tuple((np.ceil((bounds.max(0)-origin)/spacing).astype(int)+padding+1).tolist())
    if np.prod(shape)>max_voxels:raise ValueError('Grid exceeds configured memory bound')
    return Grid(origin,shape,spacing)


def voxelize(mesh,grid):
    if mesh.n_open_edges:raise ValueError('Closed source surface required')
    mask=np.zeros(grid.shape,dtype=bool)
    reference=grid.image()
    for part in mesh.split_bodies():
        surf=part.extract_surface().triangulate().clean()
        surf=pv.PolyData(surf.points,surf.faces)
        vox=surf.voxelize_binary_mask(reference_volume=reference)
        mask|=np.asarray(vox['mask']).reshape(grid.shape,order='F').astype(bool)
    return mask


def extract(field,grid,level=.5,iterations=30):
    v,f,_,_=marching_cubes(field.astype(np.float32),level,spacing=(grid.spacing,)*3,allow_degenerate=False)
    mesh=pv.PolyData(v+grid.origin,np.column_stack([np.full(len(f),3),f]).ravel()).clean()
    if iterations:mesh=mesh.smooth_taubin(n_iter=iterations,pass_band=.1)
    mesh=mesh.compute_normals(auto_orient_normals=True,consistent_normals=True)
    if mesh.n_open_edges:raise ValueError('Extracted surface is not closed')
    return mesh


def inside(points,mesh):
    return np.asarray(pv.PolyData(points).select_enclosed_points(mesh,tolerance=1e-7)['SelectedPoints']).astype(bool)


def coverage(points,mesh,axis,center,bins=10,min_cells=30):
    """Global + equal-width axial + spatial-block coverage; zero-support bins excluded."""
    points=np.asarray(points);hit=inside(points,mesh)
    axial=(points-center)@axis
    ids=np.minimum(bins-1,((axial-axial.min())/max(np.ptp(axial),1e-12)*bins).astype(int))
    span=np.maximum(np.ptp(points,axis=0),1e-12)
    block_shape=np.array([6,6,3])
    blocks=np.minimum(block_shape-1,((points-points.min(0))/span*block_shape).astype(int))
    blockids=np.ravel_multi_index(blocks.T,block_shape)
    def summarize(groups):
        return [{'id':int(k),'cells':int((groups==k).sum()),'coverage':float(hit[groups==k].mean())}
                for k in np.unique(groups)]
    axial_rows=summarize(ids);block_rows=summarize(blockids)
    supported=lambda rows:[r['coverage'] for r in rows if r['cells']>=min_cells]
    return {'global':float(hit.mean()),'cells':len(points),'inside_cells':int(hit.sum()),
            'axial':axial_rows,'spatial_blocks':block_rows,'minimum_region_cells':min_cells,
            'min_axial':min(supported(axial_rows),default=None),
            'min_spatial_block':min(supported(block_rows),default=None)},hit


def density(points,grid,sigma=3.2):
    ids=grid.indices(points);occ=np.zeros(grid.shape,dtype=np.float32)
    occ[tuple(ids.T)]=1
    return ndi.gaussian_filter(occ,sigma=sigma,mode='constant')


def regional_candidate(points,previous,grid,axis,center,quantile=.02,sigma=3.2,min_component_cells=20):
    """Restore sparse regions with smooth axial density thresholds, preserving prior volume.

    This is a candidate generator, not an anatomical acceptance decision. Axial
    calibration must be audited with independent transverse/spatial-block checks.
    """
    if not 0<quantile<.5 or sigma<=0 or min_component_cells<1:raise ValueError('Invalid refinement parameters')
    ids=grid.indices(points);field=density(points,grid,sigma)
    values=ndi.map_coordinates(field,((points-grid.origin)/grid.spacing).T,order=1)
    coord=(points-center)@axis;span=np.ptp(coord)
    anchors=np.linspace(coord.min(),coord.max(),36);window=max(span/9,grid.spacing*3)
    levels=[]
    for anchor in anchors:
        neighbors=values[np.abs(coord-anchor)<=window]
        # Empty axial windows occur between genuinely disconnected parts.
        levels.append(np.quantile(neighbors,quantile) if len(neighbors) else np.quantile(values,quantile))
    levels=ndi.gaussian_filter1d(np.asarray(levels),1.5)
    xyz=np.indices(grid.shape,dtype=np.float32)
    gridcoord=sum((grid.origin[d]+xyz[d]*grid.spacing-center[d])*axis[d] for d in range(3))
    thresholds=np.interp(gridcoord,anchors,levels)
    oldmask=voxelize(previous,grid)
    mask=(field>=thresholds)|oldmask
    cc,n=ndi.label(mask);support=np.bincount(cc[tuple(ids.T)],minlength=n+1);support[0]=0
    # Prior substantial geometry is protected; sparse supported islands are not
    # discarded solely because the largest tissue component is much larger.
    protected=np.unique(cc[oldmask]);protected=protected[protected>0]
    keep=np.union1d(np.flatnonzero(support>=min_component_cells),protected)
    mask=np.isin(cc,keep)
    mesh=extract(ndi.gaussian_filter(mask.astype(np.float32),1),grid)
    return mesh,{'quantile':quantile,'sigma_voxels':sigma,'axial_window':window,
                 'threshold_anchors':anchors.tolist(),'threshold_values':levels.tolist(),
                 'minimum_component_cells':min_component_cells,'component_cell_support':support.tolist(),
                 'components_pruned':int(n-len(keep)),'previous_volume_protected':True}


def triangle_containment(tissue,shell,max_depth=7,tolerance=1e-5):
    """Conservative signed-distance certificate for entire triangles, not just vertices.

    Distance is 1-Lipschitz. If a triangle centroid's inward distance exceeds
    every vertex's distance to that centroid, the whole triangle is inside.
    Recursively subdivide unresolved triangles; unresolved is not a pass.
    Requires an oriented, closed, non-self-intersecting shell; numerical VTK
    signed-distance evaluation is not a formal exact-arithmetic proof.
    """
    if shell.n_open_edges:raise ValueError('Containment requires a closed shell')
    tri=tissue.points[np.asarray(tissue.faces).reshape(-1,4)[:,1:]]
    tests=0;minimum=float('inf')
    for depth in range(max_depth+1):
        cent=tri.mean(1);radius=np.linalg.norm(tri-cent[:,None,:],axis=2).max(1)
        distance=np.asarray(pv.PolyData(cent).compute_implicit_distance(shell)['implicit_distance'])
        tests+=len(tri);minimum=min(minimum,float((-distance).min()))
        if np.any(distance>=-tolerance):
            return {'passed':False,'reason':'outside_or_contacting_centroid','tests':tests,'depth':depth}
        tri=tri[distance+radius>=-tolerance]
        if len(tri)==0:return {'passed':True,'tests':tests,'depth':depth,'minimum_tested_centroid_clearance':minimum}
        if depth==max_depth:break
        a,b,c=tri[:,0],tri[:,1],tri[:,2];ab=(a+b)/2;bc=(b+c)/2;ca=(c+a)/2
        tri=np.concatenate([np.stack(t,axis=1) for t in [(a,ab,ca),(ab,b,bc),(ca,bc,c),(ab,bc,ca)]])
    return {'passed':False,'reason':'unresolved_at_depth_limit','remaining':len(tri),'tests':tests}


def repair_supported_gaps(points,mesh,grid,sigma=3.2,min_neighbors=8,min_component_cells=20):
    """Locally lower density thresholds around supported missed cells, not globally.

    The neighborhood and patch radii are recorded scale choices. They do not
    classify excluded points as biological noise or delete them from the source.
    """
    hit=inside(points,mesh);ids=grid.indices(points)
    support=cKDTree(points).query_ball_point(points,r=3*sigma*grid.spacing,return_length=True)
    missed=(~hit)&(support>=min_neighbors)
    if not missed.any():return mesh,{'supported_missed_cells':0,'changed':False}
    field=density(points,grid,sigma)
    values=ndi.map_coordinates(field,((points-grid.origin)/grid.spacing).T,order=1)
    seed=np.full(grid.shape,np.inf,dtype=np.float32)
    np.minimum.at(seed,tuple(ids[missed].T),values[missed]*.5)
    # Only a three-voxel neighborhood of supported missed cells can be added.
    local=ndi.minimum_filter(seed,size=7,mode='constant',cval=np.inf)
    oldmask=voxelize(mesh,grid);mask=(field>=local)|oldmask
    cc,n=ndi.label(mask);counts=np.bincount(cc[tuple(ids.T)],minlength=n+1);counts[0]=0
    protected=np.unique(cc[oldmask]);protected=protected[protected>0]
    keep=np.union1d(protected,np.flatnonzero(counts>=min_component_cells))
    mask=np.isin(cc,keep)
    result=extract(ndi.gaussian_filter(mask.astype(np.float32),1.5),grid)
    return result,{'supported_missed_cells':int(missed.sum()),'changed':True,
                   'support_radius':3*sigma*grid.spacing,'minimum_neighbors':min_neighbors,
                   'threshold_fraction':.5,'patch_radius_voxels':3,'antialias_sigma_voxels':1.5,'minimum_component_cells':min_component_cells,
                   'component_cell_support':counts.tolist()}


def smooth_envelope(points,tissues,grid,sigma=2.5,clearance_voxels=1.8):
    """Signed-distance smoothing of point-supported body + tissue union, with margin.

    Body is a display envelope, not a mutually exclusive tissue segmentation.
    Cavities are filled only for the body. Tissue geometry is not clipped/scaled.
    """
    field=density(points,grid,3.2)
    vals=ndi.map_coordinates(field,((points-grid.origin)/grid.spacing).T,order=1)
    body_density=field>=float(np.quantile(vals,.005))
    union=np.zeros(grid.shape,dtype=bool);overlap=np.zeros(grid.shape,dtype=np.uint16)
    for tissue in tissues:
        mask=voxelize(tissue,grid);union|=mask;overlap+=mask
    base=ndi.binary_fill_holes(body_density|union)
    sdf=ndi.distance_transform_edt(base,sampling=grid.spacing)-ndi.distance_transform_edt(~base,sampling=grid.spacing)
    field=ndi.gaussian_filter(sdf,sigma=sigma)
    initial=float(field[base].min()-clearance_voxels*grid.spacing)
    attempts=[]
    for increment in (0,.75,1.5):
        level=initial-increment*grid.spacing
        occupied=field>=level
        if any(np.take(occupied,[0,-1],axis=d).any() for d in range(3)):
            raise ValueError('Envelope reaches grid boundary; increase padding')
        mesh=extract(field,grid,level=level,iterations=30)
        certs=[triangle_containment(t,mesh) for t in tissues]
        attempts.append({'level':level,'certificates':certs})
        if all(c['passed'] for c in certs):break
    else:raise ValueError('Smooth shell failed triangle containment')
    bodymask=voxelize(mesh,grid);bodycount=int(bodymask.sum())
    stats={'grid_origin':grid.origin.tolist(),'grid_shape':list(grid.shape),'spacing':grid.spacing,
           'sdf_sigma_voxels':sigma,'clearance_voxels':clearance_voxels,'attempts':attempts,
           'point_coverage':float(inside(points,mesh).mean()),
           'union_over_body_voxel_fraction':float(np.sum(union&bodymask)/bodycount),
           'overlap_over_body_voxel_fraction':float(np.sum((overlap>1)&bodymask)/bodycount),
           'uncovered_body_voxel_fraction':float(np.sum(bodymask&~union)/bodycount),
           'outside_tissue_voxels':int(np.sum(union&~bodymask)),
           'interpretation':'Voxel estimates; overlap is counted once. Envelope margin and cavities can contribute to unfilled body volume.'}
    return mesh,stats
