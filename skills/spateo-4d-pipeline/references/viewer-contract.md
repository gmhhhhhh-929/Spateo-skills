# Viewer contract

The payload is `spateo-4d-viewer/v1`: stage states, units, source/target labels, aligned and original XYZ by observation ID, optional mapping/field vectors, features, trajectory seed IDs/time arrays, actual GLM fits, full test CSVs and QC. No full expression matrix is embedded. `registrationReview` is optional for historical payloads and contains coordinate-only per-annotation diagnostics; state whether they use full cells or exported display samples.

## Shared 3D visual language

Use the companion 3D reconstruction viewer's dark background (`#080e15`), compact header, left-side collapsible groups, independent layer visibility/opacity, reset camera and PNG export. Keep technical details in Run evidence. Preserve coordinate geometry and native aspect ratios. Surfaces, if requested, come from the 3D mesh subskill and need a verified matching frame; do not overlay untransformed original surfaces after registration.

Registration defaults to separate timepoints with the same axis ranges, manual physical aspect ratios and synchronized cameras. Offer a common-frame overlay and before/after overlays. Do not independently fit/rescale each specimen: that hides size differences. Both source and target need opacity controls, stable colors and tissue filtering. Dorsal, oblique and lateral camera presets are display choices. Label the selected scientific frame (rigid/nonrigid/imported); never silently display a nonrigid result while using rigid coordinates for downstream inference.

Capabilities drive tabs. Registration is always available; mapping, field, trajectories and features require their corresponding artifacts. Imported stages are not newly computed. A display-only change must reuse verified artifacts and must not refit science.

## Arrows, paths and colors

Use compact solid 3D cone heads and thin shafts. Head length is capped independently at `min(0.25 × displayed vector length, 0.008 × source diagonal × head-size control)`; head radius is 20% of head length. Never make cone width proportional to the entire mapping displacement. Mapping defaults to full endpoint links (vector length ×1 preserves actual endpoints), 100 links and 1.5 px semi-transparent shafts. Field display uses `0.045 × source diagonal / magnitude q95`, then the independent vector-length control (0.05–2). Default field density is 200, with brighter 2.5 px shafts and subdued source cells. Trajectories default to 120 paths and 2 px lines, sharing the compact-head limit. Head size has its own 0.25–2 control and reset. Scientific magnitudes remain in the color scale; glyph geometry, density and appearance never change stored vectors, mapping endpoints or paths. Check both default and 600-arrow stress views; thinning alone cannot fix oversized heads.

Observed source and target cells keep one categorical color each in mapping/field/trajectory views. Inferred endpoints are a separate named diamond layer, disabled when changing tabs, and never presented as observed cells. Paths use model-time color and explicit heads in the direction of the saved integration order. Playback reveals saved path prefixes, not interpolated endpoint lines. Cell visibility and path visibility are independent.

Feature mode colors source cells and paths by the same selected feature. Join paths to feature values by seed-cell ID. A constant seed-cell value is propagated along each displayed path and explicitly labeled as such: it is not evaluation of the feature at predicted locations. Missing seed values must not be invented. Keep arrows visible. Offer full-range and 2–98% color limits; clipping is color-only, never data deletion. Use explicit color stop arrays to avoid named-colormap fallbacks in Plotly JS.

## Feature–gene interpretation

Preserve actual Spateo fitted mean/95% mean-confidence bounds. The default log1p display transforms observations, mean and bounds together; fitting stays on normalized expression. Offer a linear display. Do not present the log of an expected value as an expected log-expression fit. Outliers remain present and no curve is cosmetically refitted.

For every displayed candidate, show FDR q, full-cell Spearman rho, the screening family size, and whether the spline test selected it. FDR significance is not strong correlation. The exploratory `|rho| ≥ 0.3` table filter is labeled and adjustable; nonlinear associations can have low rho. If none meet the chosen cutoff, say so. A stored zero q means numerical underflow, not a literal zero probability. Gene symbols must be backed by a verified mapping; otherwise keep IDs.

Explain the spline-versus-intercept NB2 test, mean confidence interval, limited candidate panel versus genome-wide search, spatial dependence and available biological replication. Association is not causation. CSV export retains all tests, including failures and nonsignificant results.

## Verification

Use real payloads to inspect split and overlay registration, tissue selection, independent cell opacity, solid arrowheads, trajectory playback, feature-colored paths, gene switching, expression scale and CSV export. Inspect console errors and screenshots. Check an alignment-only payload separately. Generic controls/notes must not hardcode timepoint labels. A usable viewer does not certify anatomical correctness: display observed registration concerns and keep downstream results exploratory until the scientific frame is accepted.

## Reproducible registration comparison

`viewer_payload.py --payload PAYLOAD --output NEW_HTML --comparison "Legacy default" BENCHMARK_DIR` can attach one or more completed `compare_registration.py` outputs. It verifies their coordinate-file hash, joins by unique cell ID, checks original coordinates and annotations, and shows the alternate rigid coordinates only in the registration tab. Mapping, fields, trajectories and genes stay attached to the current scientific run. Different units/original coordinate values fail the original-coordinate check. Do not imply that the result selector recomputes downstream analyses.

Unavailable GLM predictions or confidence bounds export as JSON null and plot gaps, never zero. The viewer discloses missing-value counts and candidate-panel versus genome-wide scope. Stored curves prioritize FDR-selected genes then absolute Spearman correlation; CSV retains every requested gene.

## Mapping display separation

Cell mapping offers a signed Z-separation slider, initialized to zero, and Reset spacing. For display separation `g`, source points move by `−g/2`, target points and mapped endpoint markers by `+g/2`. Mapping links use the source display point and the displayed endpoint (at vector length ×1). The offset is added after vector scaling and never enters magnitude colors or saved scientific arrays. Fitted field arrows, if selected, move with the source without acquiring a false Z component. Non-mapping tabs never receive the offset. Labels and exported view settings record it. A nonzero offset switches a top-XY view to oblique so separation is visible; the user can then choose another camera. Validate positive/negative/zero spacing, endpoint consistency and resetting without modifying the payload.
