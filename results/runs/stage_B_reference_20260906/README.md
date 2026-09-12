# Stage B: published reference geometry

These files contain data extracted from the local paper figures, not reproduced bootstrap solutions. No optimization or parameter fitting was performed.

- `purSplot_points.csv`: 498 displayed Fig.3 boundary markers.
- `chiplot_points.csv`: six Fig.4 boundary series and distinct highlighted / black markers. Select `role=boundary_sample` for contour comparisons; the `epsilon` column gives the caption mapping.
- `report.json`: axis calibration, displayed extrema, marker semantics, physical reference discrepancy, and the comparison plan.
- `extracted_reference.pdf/png`: an independent redraw of those extracted marker centers. `purSplot.png` and `chiplot.png` render the original source figures.

The x and y coordinates were recovered directly from the bounding-box centers of filled circular PDF vector paths using a linear fit to every axis tick. Duplicate highlighted paths were removed. Typical tick-fit residual is at most 0.002 PDF points. A conservative 0.01 PDF point coordinate resolution is recorded separately; neither this nor the CSV decimal precision is an estimate of the original solver accuracy.

Fig.3 displayed samples span x ≈ [−2.90195, 2.23289], y ≈ [−0.733589, 0.0793227]. The Fig.4 outer-to-inner tolerances are 0.006 (blue), 0.004 (orange), 0.002 (green), 0.001 (red-orange), 0.0006 (purple), 0.0002 (brown). Their maximum displayed x values are approximately 0.161880, 0.125595, 0.0825730, 0.0547821, 0.0410012, 0.0223626. Fig.4 shows only x>0.

The actual Fig.4 black marker is approximately (0.07131334, −0.00475418). The explicit paper parameters mπ=140 MeV and fπ=92 MeV instead give (0.07332139, −0.00488809) through Eq.(2.18). Both lie on y=−x/15 within the plot resolution, but their x values differ by 2.7387%. They must remain separately labeled, without tuning the declared inputs to the drawing. Both are inside the displayed green sample hull and outside the smaller tolerances.

Section 3.1 (TeX paragraph around line 931) explicitly describes maximizing ray distance t with x=a+t cos α, y=b+t sin α from an interior center. Support maximization uses an outward normal and is an equivalent way to recover the convex projected set, but its angle does not label the same boundary point. Compare convex sets / support functions in common normals for Fig.3 and upper/lower vertical slices for the visible x>0 portion of Fig.4. Retain computed inner/outer support gaps separately from paper scan/rendering resolution and finite-prescription differences.
