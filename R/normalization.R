# SCTransform for MatisseObject is registered as an S3 method in dispatch.R,
# following the Seurat ecosystem convention. It provides mode-aware defaults:
# event mode normalises the "transcript" assay; junction mode normalises the
# active default assay (usually "RNA"). Users can override with the `assay`
# argument.
#
# See dispatch.R: SCTransform.MatisseObject
NULL
