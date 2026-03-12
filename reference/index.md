# Package index

## All functions

- [`as_dist(`*`<mle_fit>`*`)`](https://queelius.github.io/algebraic.mle/reference/as_dist.mle_fit.md)
  : Convert an MLE to a distribution object.

- [`as_dist(`*`<mle_fit_boot>`*`)`](https://queelius.github.io/algebraic.mle/reference/as_dist.mle_fit_boot.md)
  : Convert a bootstrap MLE to an empirical distribution.

- [`bias()`](https://queelius.github.io/algebraic.mle/reference/bias.md)
  : Generic method for computing the bias of an estimator object.

- [`bias(`*`<mle_fit>`*`)`](https://queelius.github.io/algebraic.mle/reference/bias.mle_fit.md)
  : Computes the bias of an \`mle_fit\` object assuming the large sample
  approximation is valid and the MLE regularity conditions are
  satisfied. In this case, the bias is zero (or zero vector).

- [`bias(`*`<mle_fit_boot>`*`)`](https://queelius.github.io/algebraic.mle/reference/bias.mle_fit_boot.md)
  : Computes the estimate of the bias of a \`mle_fit_boot\` object.

- [`cdf(`*`<mle_fit>`*`)`](https://queelius.github.io/algebraic.mle/reference/cdf.mle_fit.md)
  : CDF of the asymptotic distribution of an MLE.

- [`coef(`*`<mle_fit>`*`)`](https://queelius.github.io/algebraic.mle/reference/coef.mle_fit.md)
  :

  Extract coefficients from an `mle_fit` object.

- [`combine()`](https://queelius.github.io/algebraic.mle/reference/combine.md)
  : Combine independent MLEs for the same parameter.

- [`conditional(`*`<mle_fit>`*`)`](https://queelius.github.io/algebraic.mle/reference/conditional.mle_fit.md)
  : Conditional distribution from an MLE.

- [`confint(`*`<mle_fit>`*`)`](https://queelius.github.io/algebraic.mle/reference/confint.mle_fit.md)
  : Function to compute the confidence intervals of \`mle_fit\` objects.

- [`confint(`*`<mle_fit_boot>`*`)`](https://queelius.github.io/algebraic.mle/reference/confint.mle_fit_boot.md)
  : Method for obtained the confidence interval of an \`mle_fit_boot\`
  object. Note: This impelements the \`vcov\` method defined in
  \`stats\`.

- [`confint_from_sigma()`](https://queelius.github.io/algebraic.mle/reference/confint_from_sigma.md)
  : Function to compute the confidence intervals from a
  variance-covariance matrix

- [`density(`*`<mle_fit>`*`)`](https://queelius.github.io/algebraic.mle/reference/density.mle_fit.md)
  : PDF of the asymptotic distribution of an MLE.

- [`density(`*`<mle_fit_boot>`*`)`](https://queelius.github.io/algebraic.mle/reference/density.mle_fit_boot.md)
  : PDF of the empirical distribution of bootstrap replicates.

- [`dim(`*`<mle_fit>`*`)`](https://queelius.github.io/algebraic.mle/reference/dim.mle_fit.md)
  : Dimension (number of parameters) of an MLE.

- [`dim(`*`<mle_fit_boot>`*`)`](https://queelius.github.io/algebraic.mle/reference/dim.mle_fit_boot.md)
  : Dimension (number of parameters) of a bootstrap MLE.

- [`expectation(`*`<mle_fit>`*`)`](https://queelius.github.io/algebraic.mle/reference/expectation.mle_fit.md)
  : Expectation operator applied to \`x\` of type \`mle_fit\` with
  respect to a function \`g\`. That is, \`E(g(x))\`.

- [`inv_cdf(`*`<mle_fit>`*`)`](https://queelius.github.io/algebraic.mle/reference/inv_cdf.mle_fit.md)
  : Quantile function of the asymptotic distribution of an MLE.

- [`is_mle()`](https://queelius.github.io/algebraic.mle/reference/is_mle.md)
  : Determine if an object \`x\` is an \`mle_fit\` object.

- [`is_mle_boot()`](https://queelius.github.io/algebraic.mle/reference/is_mle_boot.md)
  : Determine if an object is an \`mle_fit_boot\` object.

- [`joint()`](https://queelius.github.io/algebraic.mle/reference/joint.md)
  : Compose independent MLEs into a joint MLE.

- [`logLik(`*`<mle_fit>`*`)`](https://queelius.github.io/algebraic.mle/reference/logLik.mle_fit.md)
  :

  Log-likelihood of an `mle_fit` object.

- [`marginal(`*`<mle_fit>`*`)`](https://queelius.github.io/algebraic.mle/reference/marginal.mle_fit.md)
  : Method for obtaining the marginal distribution of an MLE that is
  based on asymptotic assumptions:

- [`mean(`*`<mle_fit>`*`)`](https://queelius.github.io/algebraic.mle/reference/mean.mle_fit.md)
  : Mean of the asymptotic distribution of an MLE.

- [`mean(`*`<mle_fit_boot>`*`)`](https://queelius.github.io/algebraic.mle/reference/mean.mle_fit_boot.md)
  : Mean of bootstrap replicates.

- [`mle()`](https://queelius.github.io/algebraic.mle/reference/mle.md) :
  Constructor for making \`mle_fit\` objects, which provides a common
  interface for maximum likelihood estimators.

- [`mle_boot()`](https://queelius.github.io/algebraic.mle/reference/mle_boot.md)
  : Bootstrapped MLE

- [`mle_numerical()`](https://queelius.github.io/algebraic.mle/reference/mle_numerical.md)
  : This function takes the output of \`optim\`, \`newton_raphson\`, or
  \`sim_anneal\` and turns it into an \`mle_fit_numerical\` (subclass of
  \`mle_fit\`) object.

- [`mse()`](https://queelius.github.io/algebraic.mle/reference/mse.md) :
  Generic method for computing the mean squared error (MSE) of an
  estimator, \`mse(x) = E\[(x-mu)^2\]\` where \`mu\` is the true
  parameter value.

- [`mse(`*`<mle_fit>`*`)`](https://queelius.github.io/algebraic.mle/reference/mse.mle_fit.md)
  : Computes the MSE of an \`mle_fit\` object.

- [`mse(`*`<mle_fit_boot>`*`)`](https://queelius.github.io/algebraic.mle/reference/mse.mle_fit_boot.md)
  : Computes the estimate of the MSE of a \`boot\` object.

- [`nobs(`*`<mle_fit>`*`)`](https://queelius.github.io/algebraic.mle/reference/nobs.mle_fit.md)
  : Method for obtaining the number of observations in the sample used
  by an \`mle_fit\`.

- [`nobs(`*`<mle_fit_boot>`*`)`](https://queelius.github.io/algebraic.mle/reference/nobs.mle_fit_boot.md)
  : Method for obtaining the number of observations in the sample used
  by an \`mle_fit_boot\`.

- [`nparams(`*`<mle_fit>`*`)`](https://queelius.github.io/algebraic.mle/reference/nparams.mle_fit.md)
  : Method for obtaining the number of parameters of an \`mle_fit\`
  object.

- [`nparams(`*`<mle_fit_boot>`*`)`](https://queelius.github.io/algebraic.mle/reference/nparams.mle_fit_boot.md)
  : Method for obtaining the number of parameters of an \`boot\` object.

- [`obs(`*`<mle_fit>`*`)`](https://queelius.github.io/algebraic.mle/reference/obs.mle_fit.md)
  : Method for obtaining the observations used by the \`mle_fit\` object
  \`x\`.

- [`obs(`*`<mle_fit_boot>`*`)`](https://queelius.github.io/algebraic.mle/reference/obs.mle_fit_boot.md)
  : Method for obtaining the observations used by the \`mle_fit_boot\`.

- [`observed_fim()`](https://queelius.github.io/algebraic.mle/reference/observed_fim.md)
  : Generic method for computing the observed FIM of an \`mle_fit\`
  object.

- [`observed_fim(`*`<mle_fit>`*`)`](https://queelius.github.io/algebraic.mle/reference/observed_fim.mle_fit.md)
  : Function for obtaining the observed FIM of an \`mle_fit\` object.

- [`orthogonal()`](https://queelius.github.io/algebraic.mle/reference/orthogonal.md)
  : Generic method for determining the orthogonal parameters of an
  estimator.

- [`orthogonal(`*`<mle_fit>`*`)`](https://queelius.github.io/algebraic.mle/reference/orthogonal.mle_fit.md)
  : Method for determining the orthogonal components of an \`mle_fit\`
  object \`x\`.

- [`params(`*`<mle_fit>`*`)`](https://queelius.github.io/algebraic.mle/reference/params.mle_fit.md)
  : Method for obtaining the parameters of an \`mle_fit\` object.

- [`params(`*`<mle_fit_boot>`*`)`](https://queelius.github.io/algebraic.mle/reference/params.mle_fit_boot.md)
  : Method for obtaining the parameters of an \`boot\` object.

- [`pred()`](https://queelius.github.io/algebraic.mle/reference/pred.md)
  : Generic method for computing the predictive confidence interval
  given an estimator object \`x\`.

- [`pred(`*`<mle_fit>`*`)`](https://queelius.github.io/algebraic.mle/reference/pred.mle_fit.md)
  : Estimate of predictive interval of \`T\|data\` using Monte Carlo
  integration.

- [`print(`*`<mle_fit>`*`)`](https://queelius.github.io/algebraic.mle/reference/print.mle_fit.md)
  : Print method for \`mle_fit\` objects.

- [`print(`*`<summary_mle_fit>`*`)`](https://queelius.github.io/algebraic.mle/reference/print.summary_mle_fit.md)
  : Function for printing a \`summary\` object for an \`mle_fit\`
  object.

- [`rmap(`*`<mle_fit>`*`)`](https://queelius.github.io/algebraic.mle/reference/rmap.mle_fit.md)
  : Computes the distribution of \`g(x)\` where \`x\` is an \`mle_fit\`
  object.

- [`sampler(`*`<mle_fit>`*`)`](https://queelius.github.io/algebraic.mle/reference/sampler.mle_fit.md)
  : Method for sampling from an \`mle_fit\` object.

- [`sampler(`*`<mle_fit_boot>`*`)`](https://queelius.github.io/algebraic.mle/reference/sampler.mle_fit_boot.md)
  : Method for sampling from an \`mle_fit_boot\` object.

- [`score_val()`](https://queelius.github.io/algebraic.mle/reference/score_val.md)
  : Generic method for computing the score of an estimator object
  (gradient of its log-likelihood function evaluated at the MLE).

- [`score_val(`*`<mle_fit>`*`)`](https://queelius.github.io/algebraic.mle/reference/score_val.mle_fit.md)
  : Computes the score of an \`mle_fit\` object (score evaluated at the
  MLE).

- [`se()`](https://queelius.github.io/algebraic.mle/reference/se.md) :
  Generic method for obtaining the standard errors of an estimator.

- [`se(`*`<mle_fit>`*`)`](https://queelius.github.io/algebraic.mle/reference/se.mle_fit.md)
  : Function for obtaining an estimate of the standard error of the MLE
  object \`x\`.

- [`summary(`*`<mle_fit>`*`)`](https://queelius.github.io/algebraic.mle/reference/summary.mle_fit.md)
  : Function for obtaining a summary of \`object\`, which is a fitted
  \`mle_fit\` object.

- [`sup(`*`<mle_fit>`*`)`](https://queelius.github.io/algebraic.mle/reference/sup.mle_fit.md)
  : Support of the asymptotic distribution of an MLE.

- [`vcov(`*`<mle_fit>`*`)`](https://queelius.github.io/algebraic.mle/reference/vcov.mle_fit.md)
  : Computes the variance-covariance matrix of \`mle_fit\` object.

- [`vcov(`*`<mle_fit_boot>`*`)`](https://queelius.github.io/algebraic.mle/reference/vcov.mle_fit_boot.md)
  : Computes the variance-covariance matrix of \`boot\` object. Note:
  This impelements the \`vcov\` method defined in \`stats\`.
