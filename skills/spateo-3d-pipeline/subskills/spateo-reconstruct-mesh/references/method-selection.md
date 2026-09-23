# Surface reconstruction method selection

This contract is reviewed against `gmhhhhhh-929/spateo-release` commit `615644f88613bea8ceb2e2df1e2391d16de55ec1` and PyVista 0.46.5.

## Why the native method is sensitive

`st.tdr.construct_surface(..., cs_method="marching_cube")` delegates to `marching_cube_mesh`. The pinned implementation computes a global voxel width as the largest sampled nearest-neighbor distance multiplied by `mc_scale_factor`. A detached point can therefore make every voxel too coarse, while a small factor can create a very large volume and many disconnected single-cell shells. With `dist_sample_num=None`, the pairwise distance matrix is quadratic in the number of points; sampling reduces memory but can miss the problematic region. This wrapper seeds that legacy NumPy sampling deterministically for each mesh name and records both the base and derived seed.

The source rasterizes unsmoothed impulses and pads only the maximum grid boundary. Its default `levelset=0` is unsuitable for a binary 0/1 field; use `0.5`. `smooth` is applied later as Laplacian mesh smoothing, so thousands of iterations can shrink the model and erase thin anatomy rather than repair the scalar field.

The full `construct_surface` then splits bodies, calls MeshFix and ACVD on each, and chooses the larger of its inside/outside point sets. On the validated 16.53 h Drosophila body, mixed connectivity from `split_bodies` caused MeshFix to raise a reshape error before any surface was returned. The native mode therefore calls the pinned Spateo marching-cubes core directly, converts it to triangle-only PolyData, and performs optional face-only repair without claiming equivalence to the full post-processing chain. Repair explicitly disables MeshFix's default removal of smaller components and records topology before and after, because paired or disconnected anatomy can be valid.

## Default density route

Spatial-transcriptomic cell centroids sample the tissue volume, not only its boundary. The default route therefore:

1. creates one padded isotropic coordinate grid shared by body and annotation outputs;
2. maps points to their nearest grid nodes and caps occupancy at one per node;
3. applies a Gaussian filter to form a smooth volume field;
4. searches descending normalized isovalues and chooses the first that meets the requested source-point voxel coverage after component filtering, failing rather than silently accepting an unmet automatic target;
5. labels density support with face-connected (6-neighbor) voxels and retains components above both an absolute voxel count and a fraction of the largest component, preserving substantial paired anatomy without treating corner contact as a connected surface;
6. extracts a Lewiner marching-cubes surface with physical voxel spacing;
7. applies Taubin smoothing, which reduces shrinkage relative to repeated Laplacian smoothing;
8. repairs remaining open edges one connected surface at a time, preserves valid smaller components, and explicitly records any numerically degenerate fragment small enough to drop under `--repair-max-degenerate-cells`, then saves and reloads through Spateo;
9. requires a configurable minimum final point-containment fraction for density results and publishes all outputs only after every target passes.

The component filter changes only the derived density volume. It does not remove points from the input VTK or observations from the source H5AD. The manifest records pre-filter component sizes, retained components, chosen isovalue, voxel coverage, final point containment, and repair activity.

## Parameter effects

| Parameter | Increasing it usually does | Main risk |
| --- | --- | --- |
| `target_long_axis` | Uses smaller voxels and preserves more detail | Cubic memory/time growth and fragmented sparse surfaces |
| `voxel_size` | Uses coarser geometry | Loss of thin structures and false bridges |
| `sigma` | Connects slice gaps and smooths occupancy | Merges nearby structures and inflates surfaces |
| coverage target | Includes more source points by lowering the isovalue | Admits low-support islands or broad halos |
| `min_component_fraction` | Removes more small components | Deletes legitimate bilateral/disconnected anatomy |
| `smooth` in density mode | Adds Taubin iterations | Excessive rounding despite low shrinkage |
| `mc_scale_factor` in native mode | Coarsens the global occupancy grid | Inflated, bridged geometry and detail loss |
| `smooth` in native mode | Adds Laplacian iterations | Volume shrinkage and loss of branches |

Choose parameters independently for a full body and each annotation class when necessary. Keep the input point cloud immutable and write a new result directory for every comparison.

## Alternatives reviewed

- PyVista `gaussian_splatting` implements a related point-to-volume operation but was added only in PyVista 0.46; the SciPy implementation keeps this skill compatible with the wider Spateo-supported PyVista range.
- PyVista `reconstruct_surface` assumes points lie on the surface of a solid and is not the default for volume-filling cell centroids.
- Open3D ball pivoting and Poisson reconstruction require normals and surface samples. Poisson may extrapolate surfaces in low-density areas even for suitable input.
- Alpha shapes avoid normals but remain strongly dependent on sampling density, scale, and outliers.

Primary references: [PyVista Gaussian splatting](https://docs.pyvista.org/api/core/_autosummary/pyvista.DataSetFilters.gaussian_splatting.html), [PyVista surface reconstruction](https://docs.pyvista.org/api/core/_autosummary/pyvista.polydatafilters.reconstruct_surface), [PyVista Taubin smoothing](https://docs.pyvista.org/api/core/_autosummary/pyvista.polydatafilters.smooth_taubin), [scikit-image marching cubes](https://scikit-image.org/docs/stable/api/skimage.measure.html#skimage.measure.marching_cubes), and [Open3D surface reconstruction](https://www.open3d.org/docs/latest/tutorial/geometry/surface_reconstruction.html).
