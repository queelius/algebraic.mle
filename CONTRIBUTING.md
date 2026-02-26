# Contributing to algebraic.mle

Thanks for your interest in contributing to algebraic.mle!

## Reporting Bugs

Please open an issue at <https://github.com/queelius/algebraic.mle/issues> with:

- A minimal reproducible example
- Expected vs actual behavior
- Output of `sessionInfo()`

## Proposing Changes

1. Fork the repository and create a branch from `master`
2. Make your changes
3. Run `devtools::check()` to ensure no new issues are introduced
4. Run `devtools::test()` to verify tests pass
5. Submit a pull request

## Style

- Follow existing code conventions (S3 classes, roxygen2 documentation)
- Add tests for new functionality in `tests/testthat/`
- Update documentation with `devtools::document()` after changing roxygen comments

## Code of Conduct

Please note that this project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md).
By participating in this project you agree to abide by its terms.
