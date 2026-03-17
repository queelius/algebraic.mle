# Design: MLE Composition Algebra + Documentation Reframe

**Date:** 2026-02-27
**Status:** Approved
**Version target:** 1.2.0

## Mission Statement

The MLE is a technology. Under regularity conditions, any MLE is asymptotically
normal with known variance. `algebraic.mle` exploits this: once you have an MLE
-- however you found it -- this package gives you an algebra of composition,
transformation, and inference.

The package does not solve for MLEs (that's `compositional.mle` / `likelihood.model`).
It takes MLE objects as input and provides algebraic operations closed over them.

## New Features

### 1. `joint()` -- Block-Diagonal Composition

Combines independent MLEs with **disjoint** parameter sets into a joint MLE.

**Signature:**
```r
joint(x, ...) # generic
joint.mle(x, ...) # method: accepts 2+ mle objects
```

**Semantics:**
- `theta.hat` = concatenation of all parameter vectors
- `sigma` = block-diagonal matrix from individual vcov matrices
- `loglike` = sum of individual log-likelihoods (when all available, else NULL)
- `nobs` = NULL (different experiments; no single sample size)
- `info` = block-diagonal from individual FIMs (when all available, else NULL)
- `score` = concatenation of score vectors (when all available, else NULL)

**Validation:**
- Parameter names must be disjoint (error on overlap)
- All inputs must be `mle` objects
- At least 2 inputs required

**Result type:** `mle` (plain base class, no subclass)

**Algebraic closure:** After `joint()`, existing operations work:
- `marginal(joint_mle, indices)` -> recovers component MLEs
- `conditional(joint_mle, ...)` -> conditional inference via Schur complement
- `rmap(joint_mle, g)` -> delta method on the joint parameter space
- `as_dist(joint_mle)` -> MVN for distribution algebra

### 2. `combine()` -- Optimal Same-Parameter Combination

Combines independent MLEs for the **same** parameter via inverse-variance weighting.

**Signature:**
```r
combine(x, ...) # generic
combine.mle(x, ...) # method: accepts 2+ mle objects or a list
```

**Semantics (same as mle_weighted, improved):**
- Weighted average: `theta_combined = (sum I_i)^{-1} sum(I_i theta_i)`
- `sigma` = `ginv(sum I_i)`
- `info` = `sum I_i`
- `nobs` = sum of all nobs
- `loglike` = NULL (cannot meaningfully combine)

**Improvements over mle_weighted:**
- Accepts varargs: `combine(a, b, c)` instead of `combine(list(a, b, c))`
- Also accepts a list for programmatic use: `combine(list_of_mles)`
- Falls back to `ginv(vcov)` when FIM is not directly available
- Single MLE passthrough: `combine(x)` returns `x`

**Backwards compatibility:** `mle_weighted()` stays unchanged. `combine()` is
the recommended API going forward.

### 3. Documentation Reframe

**DESCRIPTION:** Rewrite package description with "MLE as technology" narrative.

**README.Rmd:** Restructure around:
1. What is the MLE technology? (asymptotic normality, invariance, Fisher info)
2. What can you do with it? (compose, transform, bridge to dist algebra, infer)
3. How do you create MLE objects? (mle_numerical, mle_boot, mle)

**New vignette `vignettes/mle-algebra.Rmd`:** "The Algebra of MLEs"
- Motivating example: component reliability -> system reliability via joint + rmap
- Full pipeline: independent MLEs -> joint -> transform -> inference
- Meta-analysis via combine()
- as_dist() bridge to distribution algebra

## Files Changed

| File | Action | Description |
|------|--------|-------------|
| `R/joint.R` | CREATE | joint() generic + joint.mle() |
| `R/combine.R` | CREATE | combine() generic + combine.mle() |
| `R/algebraic.mle.R` | EDIT | Export joint, combine |
| `DESCRIPTION` | EDIT | Bump to 1.2.0, rewrite Description |
| `README.Rmd` | EDIT | Restructure with new narrative |
| `vignettes/mle-algebra.Rmd` | CREATE | New algebra vignette |
| `tests/testthat/test-joint.R` | CREATE | Tests for joint() |
| `tests/testthat/test-combine.R` | CREATE | Tests for combine() |

## What We're NOT Doing

- No MLE arithmetic operators (+.mle, *.mle)
- No shared-parameter joint() (disjoint only)
- No changes to likelihood.model or compositional.mle
- No pruning of unused features
- No breaking changes to any existing API
