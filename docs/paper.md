# Summary

`algebraic.mle` is an R package that defines an algebra over maximum
likelihood estimators (MLEs). It provides a unified S3 class hierarchy —
`mle` as the base type, with specializations `mle_numerical`,
`mle_boot`, `rmap_mle`, and `mle_weighted` — together with a rich set of
generic methods for statistical inference. The package treats MLEs as
first-class objects — in the sense of @casella2002 and @lehmann1998 —
that can be inspected, transformed, combined, and compared through a
consistent interface.

Given any fitted MLE, `algebraic.mle` provides parameter extraction,
standard errors and variance-covariance matrices, asymptotic and
bootstrap confidence intervals, model comparison via AIC, likelihood
ratio statistics, bias and mean squared error estimation, predictive
intervals through Monte Carlo integration, and parameter transformations
via the delta method or Monte Carlo simulation. A typical workflow wraps
the output of [`stats::optim()`](https://rdrr.io/r/stats/optim.html)
into an `mle_numerical` object, after which the full suite of methods
becomes available:

``` r
sol <- mle_numerical(optim(
  par = c(b0 = 0, b1 = 0), fn = loglik,
  control = list(fnscale = -1), hessian = TRUE))
confint(sol)
rmap(sol, g, method = "delta")
```

The package is implemented in R \[@R2024\] and is available on CRAN. It
builds on `algebraic.dist` \[@algebraic.dist\], inheriting generics for
`params`, `sampler`, `rmap`, `expectation`, and `marginal`, so that MLEs
and probability distributions share a common algebraic vocabulary.

# Statement of Need

Researchers who fit custom parametric models by maximum likelihood
frequently need to perform post-estimation inference — computing
confidence intervals, propagating uncertainty through transformations,
comparing nested models, or generating predictive intervals that account
for parameter uncertainty. In R,
[`stats::optim()`](https://rdrr.io/r/stats/optim.html) returns a raw
list with point estimates and an optional Hessian, leaving the analyst
to manually invert the Hessian, compute standard errors, and implement
the delta method each time.

`algebraic.mle` fills this gap for researchers in statistics,
reliability engineering, and applied sciences who build custom
likelihood-based models and need a lightweight, composable abstraction
for MLE results. It targets users who work below the level of
formula-based model fitting (as in `glm` or `lme4`) but above raw
numerical optimization, providing the “missing middle” where the output
of any optimizer can be elevated to a statistically rich object with a
uniform API.

# State of the Field

Several R packages address maximum likelihood estimation, each with a
different scope. The built-in
[`stats4::mle`](https://rdrr.io/r/stats4/mle.html) class \[@stats4\]
couples optimization and inference into a single S4 framework; users
must specify models through its `minuslogl` interface, and the resulting
S4 objects are not easily composed with other estimators. The `bbmle`
package \[@bbmle\] extends `stats4` with `mle2`, adding AIC tables,
profile likelihood, and multi-model comparison, but remains tied to its
own S4 class system. The `maxLik` package \[@maxLik\] provides a
comprehensive optimization back-end supporting Newton-Raphson, BHHH, and
constrained methods, returning `maxLik` objects with summary and
covariance methods; its focus is on the optimization step rather than on
algebraic manipulation of the resulting estimator. The `mle.tools`
package \[@mle.tools\] focuses specifically on Cox-Snell bias correction
and Fisher information computation for known distribution families. The
`MLE` package \[@MLE\] offers ready-made MLE functions for over 100
named distributions but does not provide a general-purpose estimator
class. The `fitdistrplus` package \[@fitdistrplus\] is widely used for
fitting parametric distributions by MLE, moment matching, and other
methods; it returns `fitdist` objects with summary and plotting methods,
but targets distribution fitting rather than algebraic manipulation of
arbitrary estimators.

`algebraic.mle` takes a different approach. Rather than providing its
own optimizer or restricting attention to particular distribution
families, it defines a post-estimation algebra: any optimization result
can be wrapped into an `mle` object, and from that point forward the
user has access to variance estimation, confidence intervals, delta
method transformations (`rmap`), bootstrap inference (`mle_boot`),
inverse-variance weighted combination (`mle_weighted`), marginal
distributions, and predictive intervals. This separation of optimization
from inference makes the package composable with arbitrary optimizers
and model specifications, including the downstream reliability packages
`maskedcauses` and `maskedhaz` that use `algebraic.mle` for MLE results
from masked failure data models.

# Software Design

The package is organized around an S3 class hierarchy rooted at `mle`.
The base constructor accepts a parameter vector, log-likelihood value,
score, variance- covariance matrix, and observed Fisher information
matrix, producing an object on which all generic methods dispatch. Three
subclasses extend this base: `mle_numerical` wraps the return value of
[`stats::optim()`](https://rdrr.io/r/stats/optim.html), computing the
variance-covariance as the inverse of the negative Hessian; `mle_boot`
wraps a [`boot::boot`](https://rdrr.io/pkg/boot/man/boot.html) object,
estimating bias and variance from bootstrap replicates \[@efron1993\]
rather than asymptotic theory; and `rmap_mle` represents a transformed
MLE $g\left( \widehat{\theta} \right)$, with variance estimated either
by the delta method (using the Jacobian from
[`numDeriv::jacobian`](https://rdrr.io/pkg/numDeriv/man/jacobian.html))
or by Monte Carlo simulation.

A fourth constructor, `mle_weighted`, combines a list of independent
MLEs of the same parameter via inverse-variance weighting, producing a
single MLE whose Fisher information is the sum of the individual
information matrices. All subclasses satisfy the `mle` interface, so
methods like `confint`, `vcov`, `aic`, `sampler`, `pred`, and `rmap`
work uniformly. The package re-exports generics from `algebraic.dist`
(including `params`, `sampler`, `rmap`, `expectation`, and `marginal`)
and provides S3 methods for base R generics (`coef`, `logLik`, `vcov`,
`confint`, `nobs`), enabling interoperability with
[`AIC()`](https://rdrr.io/r/stats/AIC.html) and
[`BIC()`](https://rdrr.io/r/stats/AIC.html).

# Research Impact Statement

`algebraic.mle` serves as the inference layer for a family of R packages
addressing reliability engineering with masked failure data. The
packages `maskedcauses` and `maskedhaz` use `mle` and `mle_numerical`
objects to represent fitted series system models, and
`compositional.mle` \[@compositional.mle\] builds composable MLE solvers
that return `mle` objects. The package is published on CRAN and
available through the author’s
[r-universe](https://queelius.r-universe.dev). It is used in ongoing
research on model selection for masked series systems at Southern
Illinois University Edwardsville.

# AI Usage Disclosure

Claude (Anthropic) was used to assist with drafting and structuring this
JOSS paper. The author reviewed, edited, and validated all text. All
software design decisions, implementations, and mathematical content are
the author’s own work. No generative AI was used in the development of
the `algebraic.mle` package source code.

# Acknowledgements

The author thanks Southern Illinois University Edwardsville for
institutional support.

# References
