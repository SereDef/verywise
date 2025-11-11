# Locate vertices belonging to specific cortical ROIs

This function identifies which vertices in a cortical surface mesh
correspond to one or more specified regions of interest (ROIs), based on
the 36 regions of the Desikan-Killiany atlas distributed with FreeSurfer
(`aparc.annot`). If no ROIs are specified, it returns a full lookup
table of all ROIs with vertex counts, proportions, and lobe assignments.

## Usage

``` r
locate_roi(rois = NULL, n_verts = 163842, hemi = c("lh", "rh"), verbose = TRUE)
```

## Arguments

- rois:

  Optional character vector of ROI names to locate (e.g.,
  `c("insula", "precentral")`). If `NULL` (default), the function
  returns a full ROI lookup table.

- n_verts:

  Integer. Total number of surface vertices to consider (default:
  `163842` "fsaverage").

- hemi:

  String: "lh" or "rh". Used to point to the hemisphere-specific aparc
  annotation file.

- verbose:

  Logical (default: `TRUE`).

## Value

- If `rois` is `NULL`, returns a `data.frame` with:

  - `roi_id`: ROI numeric ID

  - `roi_label`: ROI name

  - `vw_count`: Number of vertices in the ROI

  - `vw_prop`: Proportion of total vertices in the ROI

  - `roi_lobe`: Anatomical lobe classification

  - `lobe_count`: Total number of vertices in the lobe

  - `lobe_prop`: Proportion of vertices in the lobe

- If `rois` is provided, returns a logical vector of length `n_verts`,
  where `TRUE` indicates the vertex belongs to one of the selected ROIs.

## Examples

``` r
# Get the full ROI lookup table
roi_table <- locate_roi()

# Get a logical vector of vertices in the insula
insula_mask <- locate_roi(rois = "insula")
#>  * selected 5229 vertices (3.2%)
```
