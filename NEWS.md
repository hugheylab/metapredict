# metapredict 1.1.3
- Fixed issues caught by lintr.

# metapredict 1.1.2
- Tweaked usage of `data.table` to avoid warnings.

# metapredict 1.1.1
- Revised code to not need `globalVariables()` in order to pass R CMD check without notes.

# metapredict 1.1.0
- Switched to data.table under the hood.

# metapredict 1.0.4
- Updated links in vignette to download data from the Bhattacharjee study.

# metapredict 1.0.3
- Fixed a bug in mapping probes to genes for some platforms that was caused by mixing up characters and integers.
- Enabled support for gzipped series matrix files.

# metapredict 1.0.2
- Vignettes now properly link to each other and to documentation.

# metapredict 1.0.1
- Added `pkgdown` site.
- Updated roxygen-based documentation.
- `installCustomCdfPackages` now uses latest BrainArray version.
- `metapredict` now works with updated `glmnet::predict.glmnet`.
