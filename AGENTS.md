# Repository Guidelines

## Project Structure

This repository is an **R package** (`algebraic.mle`).

- `R/`: core implementation (S3 classes and methods for MLE objects)
- `man/`: generated Rd documentation (do not edit by hand)
- `vignettes/`: package tutorials/examples (`*.Rmd`)
- `docs/`: pkgdown site output (generated; deployed via GitHub Actions)
- `fixing/`: experimental/WIP code (excluded from builds via `.Rbuildignore`)

## Build, Test, and Development Commands

Typical local workflow uses `devtools`:

- Install dependency + this package:
  - `R -q -e 'devtools::install_github("queelius/algebraic.dist")'`
  - `R -q -e 'devtools::install()'`
- Load without installing (fast iteration): `R -q -e 'devtools::load_all()'`
- Regenerate docs/NAMESPACE from roxygen: `R -q -e 'devtools::document()'`
- Validate package (R CMD check): `R -q -e 'devtools::check()'`
- Build site locally: `R -q -e 'pkgdown::build_site()'`
- Render a vignette: `R -q -e 'rmarkdown::render("vignettes/statistics.Rmd")'`

## Coding Style & Naming Conventions

- Indentation: 4 spaces; use `<-` for assignment.
- Names: `snake_case` for functions/objects; S3 methods use `generic.class` (e.g., `print.mle`).
- Documentation: write roxygen2 comments in `R/`; run `devtools::document()` to update `man/` and `NAMESPACE`.

## Testing Guidelines

There is currently **no** `tests/` directory. Treat `devtools::check()` plus vignette execution as the primary validation path.

If adding unit tests, prefer `testthat` (edition 3 is configured in `DESCRIPTION`) with `tests/testthat/test-*.R`.

## Commit & Pull Request Guidelines

- Commit subjects in this repo are generally short and imperative (e.g., “Fix…”, “Add…”, “Remove…”). Prefer descriptive messages over “updates”.
- PRs should include: a clear description, the motivation/linked issue (if any), and confirmation that `devtools::check()` is clean.
- When changing the public API, update roxygen docs and relevant vignettes/README, and ensure pkgdown still builds.

