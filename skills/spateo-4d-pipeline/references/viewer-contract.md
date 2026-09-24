# Viewer contract

The payload is versioned `spateo-4d-viewer/v1`. It contains stage statuses, frame unit, source/target labels, full aligned and original XYZ by observation ID, optional mapped endpoints/raw and fitted vectors, scalar features, trajectory seed IDs/time arrays, actual GLM curves and CSV tables, and QC summaries. It contains no full expression matrix.

Capabilities drive tabs: registration always; mapping when mapped data exists; field only with a fitted field; trajectory only with a trajectory checkpoint; features only after feature computation. Disabled stages remain explicit in the evidence page. Imported field runs identify their mapping/field as imported, not newly computed.

Coordinate values in the payload are preserved at Python float precision. Source and target share one scene with aspectmode=data. Before/after changes the display array without changing H5AD. Scientific raw displacement satisfies endpoint minus source coordinates. Fitted field is a separate vector source.

Arrow default factor = 0.04 × source bounding-box diagonal / positive vector-magnitude q95. The user multiplier scales this factor only. Render explicit shafts and arrowheads; camera/length controls never alter vector values, field fitting, or integration. Keep direction meaningful even when magnitudes vary. The notebook uses 0.015 of its model extent; this viewer's default is an intentional display choice and is not numerical equivalence to the notebook.

Trajectory payload orientation is paths × time × XYZ from the audited native backend. Playback reveals prefixes of saved paths, not straight-line endpoint interpolation. Quantify source-envelope exits in QC. Current trajectory display uses an explicit deterministic seed cap; all scientific trajectories remain in H5AD/NPZ. Feature values and GLM curves come from their own checkpoint rather than recomputation in JavaScript.

A display-only request should use an existing payload or a verified parent; do not rebuild mapping. The legacy standalone dashboard remains separate. Mesh surface reconstruction is delegated to the 3D skill when requested; this default viewer is cell/field/trajectory based. A future imported surface must carry a verified matching frame before overlay; an original-frame surface cannot be casually reused after registration.
