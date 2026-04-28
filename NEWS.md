# Matisse 0.1.0

## Phase 6 (April 2026): unify event input, drop junction_data

* **`CreateMatisseObject()` now takes a single `events` parameter** that
  accepts either a SUPPA2 `.ioe` path (or vector of paths) or a pre-built
  `data.frame`. The previous `ioe_files` and `event_data` parameters are
  removed — Matisse auto-dispatches on input class. The IOE-precedence
  warning that gated the two-parameter form (P14) is gone, since the
  conflict isn't possible anymore.
* **`junction_data` parameter dropped (P7).** Per-junction genomic
  coordinates are auto-parsed from junction IDs at construction and
  written to `seurat[["isoform"]][[]]` (the isoform assay's
  feature-metadata table). Supported ID formats:
  `chr-start-end-strand`, `chr:start-end:strand`, `chr_start_end_strand`.
  Names that don't match get NA coords with a warning — counts still
  work, only sashimi plots degrade.
* **Renamed event-annotation columns:** `inclusion_junctions` /
  `exclusion_junctions` → `inclusion_features` / `exclusion_features`
  (mode-agnostic — they hold junction IDs in junction mode, transcript
  IDs in transcript mode).
* `PlotSashimi` now reads junction coords from
  `seurat[["isoform"]][[]]` instead of `obj@misc[["junction_data"]]`.

## Code-review pass (April 2026): API surface tightened, several bug fixes

This release applies a multi-phase code-review pass that simplifies the
API surface, fixes four real correctness bugs, and refactors event
annotation storage so it lives natively inside the embedded Seurat
object.

### Architectural changes

* Per-event annotation now lives in `seurat[["psi"]]@meta.data` (the
  Assay5 feature-metadata table), accessed via the `[[ ]]` API.
  Previously stored in `MatisseObject@misc[["event_data"]]`. Single
  source of truth — annotation rides with the assay through any
  subset/merge automatically.
* `CreateMatisseObject()` now folds the PSI calculation into
  construction by default. The two-step "construct then `CalculatePSI`"
  flow collapses to one. Pass `defer_psi = TRUE` to skip the PSI step
  (rare; for power users who want to inspect intermediate state).
  `min_coverage` parameter added to the constructor.
* `MatisseObject@misc[["event_data_path"]]` records the source path of
  the event annotation (normalized) for provenance, or `NA_character_`
  when annotation was passed as a `data.frame` directly.

### Breaking API changes

* **Removed accessor wrappers** (use the native Seurat API instead):
  - `GetJunctionCounts(obj)`  → `Matrix::t(SeuratObject::GetAssayData(GetSeurat(obj)[["isoform"]], "counts"))`
  - `GetTranscriptCounts(obj)` → `SeuratObject::GetAssayData(GetSeurat(obj)[["isoform"]], "counts")`
  - `GetInclusionCounts(obj)`  → `Matrix::t(SeuratObject::GetAssayData(GetSeurat(obj)[["psi"]], "counts"))`
  - `GetExclusionCounts(obj)`  → `Matrix::t(SeuratObject::GetAssayData(GetSeurat(obj)[["psi"]], "exclusion"))`
  - `GetEventData(obj)`        → `GetSeurat(obj)[["psi"]][[]]` (post-CalculatePSI)
  - `GetJunctionData(obj)`     → `obj@misc[["junction_data"]]`
  - `AddIsoformMetadata(obj, df)` → `SeuratObject::AddMetaData(matisse_obj, df)` (S3 method delegates)
* **Removed `BuildSimpleEvents()`.** Matisse no longer ships a heuristic
  event caller. Junction-mode users must bring an external event
  annotation (SUPPA2 IOE, rMATS, or equivalent), as transcript-mode
  already requires `ioe_files`.
* **Single-event MatisseObjects are no longer representable.**
  SeuratObject's `Assay5` requires at least 2 features, so the `[`
  operator now errors when an event subset would leave fewer than 2
  events. For single-event values, use `GetPSI(obj)[, event_id]`
  directly.

### Bug fixes

* `MatisseMeta(obj) <- value` previously wrote columns by positional
  index, silently writing the wrong values when `value`'s rownames
  disagreed with `meta.data`'s ordering. Now delegates to
  `SeuratObject::AddMetaData` for proper barcode alignment.
* `rlang::abort()` was being called with multiple positional string
  arguments at 8 sites; this silently truncated the message and
  consumed the second string as the error class. All sites collapsed
  to a single string so the full message reaches the user.
* `cli::cli_alert_*()` had the same multi-arg truncation bug at 5
  sites (the second string is `id`, not message-continuation). Fixed.
* `SCTransform.MatisseObject` referenced `object@mode == "event"`,
  but the slot is `input.mode` and valid values are
  `"junction"` / `"transcript"`. The function crashed on every call;
  it was untested. The default assay in transcript mode also pointed
  to a non-existent `"transcript"` assay (the actual assay name is
  `"isoform"`). Both fixed; new tests cover both modes.

### Removed dead code

* `.psi_rowmeans_sparse` (utils.R) — only caller was its own test.
* `.parse_junction_list` (psi.R) — never called; `.calculate_psi_matrix`
  inlined the same logic.
* `na_fill` parameter — declared on `CalculatePSI` and threaded through
  the implementation, but never actually applied. Removed in P5.

# Matisse 0.1.0 (initial release)

* New `MatisseObject` S4 class wrapping a `Seurat` object with isoform-resolved
  layers: junction counts, PSI matrix, inclusion/exclusion counts, event and
  junction annotation tables, and per-cell isoform metadata.
* `CreateMatisseObject()`: primary constructor with input validation and
  automatic cell-barcode alignment.
* `CalculatePSI()`: computes PSI matrices from junction count data and splice
  event definitions. Supports both `MatisseObject` and bare matrix input.
* `SummarizePSI()`: per-event summary statistics (mean, median, sd, coverage).
* `FilterCells()` and `FilterEvents()`: threshold-based filtering with
  informative removal summaries.
* Visualization: `PlotUMAP()`, `PlotViolin()`, `PlotHeatmap()`,
  `PlotSashimi()`.
* `MergeMatisse()`: concatenate two `MatisseObject`s with cell-prefix deduplication.
* GitHub Actions workflows for R CMD check (multi-OS), pkgdown website
  deployment, and test coverage reporting.
