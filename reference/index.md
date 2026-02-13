# Package index

## All functions

- [`aic()`](https://queelius.github.io/algebraic.mle/reference/aic.md) :
  Generic method for obtaining the AIC of a fitted distribution object
  fit.
- [`aic(`*`<mle>`*`)`](https://queelius.github.io/algebraic.mle/reference/aic.mle.md)
  : Method for obtaining the AIC of an \`mle\` object.
- [`algebraic.mle-package`](https://queelius.github.io/algebraic.mle/reference/algebraic.mle.md)
  [`algebraic.mle`](https://queelius.github.io/algebraic.mle/reference/algebraic.mle.md)
  : \`algebraic.mle\`: A package for algebraically operating on and
  generating maximum likelihood estimators from existing maximum
  likelihood estimators.
- [`bias()`](https://queelius.github.io/algebraic.mle/reference/bias.md)
  : Generic method for computing the bias of an estimator object.
- [`bias(`*`<mle>`*`)`](https://queelius.github.io/algebraic.mle/reference/bias.mle.md)
  : Computes the bias of an \`mle\` object assuming the large sample
  approximation is valid and the MLE regularity conditions are
  satisfied. In this case, the bias is zero (or zero vector).
- [`bias(`*`<mle_boot>`*`)`](https://queelius.github.io/algebraic.mle/reference/bias.mle_boot.md)
  : Computes the estimate of the bias of a \`mle_boot\` object.
- [`coef(`*`<mle>`*`)`](https://queelius.github.io/algebraic.mle/reference/coef.mle.md)
  : Extract coefficients from an \`mle\` object.
- [`confint(`*`<mle>`*`)`](https://queelius.github.io/algebraic.mle/reference/confint.mle.md)
  : Function to compute the confidence intervals of \`mle\` objects.
- [`confint(`*`<mle_boot>`*`)`](https://queelius.github.io/algebraic.mle/reference/confint.mle_boot.md)
  : Method for obtained the confidence interval of an \`mle_boot\`
  object. Note: This impelements the \`vcov\` method defined in
  \`stats\`.
- [`confint_from_sigma()`](https://queelius.github.io/algebraic.mle/reference/confint_from_sigma.md)
  : Function to compute the confidence intervals from a
  variance-covariance matrix
- [`expectation(`*`<mle>`*`)`](https://queelius.github.io/algebraic.mle/reference/expectation.mle.md)
  : Expectation operator applied to \`x\` of type \`mle\` with respect
  to a function \`g\`. That is, \`E(g(x))\`.
- [`is_mle()`](https://queelius.github.io/algebraic.mle/reference/is_mle.md)
  : Determine if an object \`x\` is an \`mle\` object.
- [`is_mle_boot()`](https://queelius.github.io/algebraic.mle/reference/is_mle_boot.md)
  : Determine if an object is an \`mle_boot\` object.
- [`logLik(`*`<mle>`*`)`](https://queelius.github.io/algebraic.mle/reference/logLik.mle.md)
  : Extract log-likelihood from an \`mle\` object.
- [`loglik_val()`](https://queelius.github.io/algebraic.mle/reference/loglik_val.md)
  : Generic method for obtaining the log-likelihood value of a fitted
  MLE object.
- [`loglik_val(`*`<mle>`*`)`](https://queelius.github.io/algebraic.mle/reference/loglik_val.mle.md)
  : Method for obtaining the log-likelihood of an \`mle\` object.
- [`marginal(`*`<mle>`*`)`](https://queelius.github.io/algebraic.mle/reference/marginal.mle.md)
  : Method for obtaining the marginal distribution of an MLE that is
  based on asymptotic assumptions:
- [`mle()`](https://queelius.github.io/algebraic.mle/reference/mle.md) :
  Constructor for making \`mle\` objects, which provides a common
  interface for maximum likelihood estimators.
- [`mle_boot()`](https://queelius.github.io/algebraic.mle/reference/mle_boot.md)
  : Bootstrapped MLE
- [`mle_numerical()`](https://queelius.github.io/algebraic.mle/reference/mle_numerical.md)
  : This function takes the output of \`optim\`, \`newton_raphson\`, or
  \`sim_anneal\` and turns it into an \`mle_numerical\` (subclass of
  \`mle\`) object.
- [`mle_weighted()`](https://queelius.github.io/algebraic.mle/reference/mle_weighted.md)
  : Accepts a list of \`mle\` objects for some parameter, say \`theta\`,
  and combines them into a single estimator \`mle_weighted\`.
- [`mse()`](https://queelius.github.io/algebraic.mle/reference/mse.md) :
  Generic method for computing the mean squared error (MSE) of an
  estimator, \`mse(x) = E\[(x-mu)^2\]\` where \`mu\` is the true
  parameter value.
- [`mse(`*`<mle>`*`)`](https://queelius.github.io/algebraic.mle/reference/mse.mle.md)
  : Computes the MSE of an \`mle\` object.
- [`mse(`*`<mle_boot>`*`)`](https://queelius.github.io/algebraic.mle/reference/mse.mle_boot.md)
  : Computes the estimate of the MSE of a \`boot\` object.
- [`nobs(`*`<mle>`*`)`](https://queelius.github.io/algebraic.mle/reference/nobs.mle.md)
  : Method for obtaining the number of observations in the sample used
  by an \`mle\`.
- [`nobs(`*`<mle_boot>`*`)`](https://queelius.github.io/algebraic.mle/reference/nobs.mle_boot.md)
  : Method for obtaining the number of observations in the sample used
  by an \`mle\`.
- [`nparams(`*`<mle>`*`)`](https://queelius.github.io/algebraic.mle/reference/nparams.mle.md)
  : Method for obtaining the number of parameters of an \`mle\` object.
- [`nparams(`*`<mle_boot>`*`)`](https://queelius.github.io/algebraic.mle/reference/nparams.mle_boot.md)
  : Method for obtaining the number of parameters of an \`boot\` object.
- [`obs(`*`<mle>`*`)`](https://queelius.github.io/algebraic.mle/reference/obs.mle.md)
  : Method for obtaining the observations used by the \`mle\` object
  \`x\`.
- [`obs(`*`<mle_boot>`*`)`](https://queelius.github.io/algebraic.mle/reference/obs.mle_boot.md)
  : Method for obtaining the observations used by the \`mle\`.
- [`observed_fim()`](https://queelius.github.io/algebraic.mle/reference/observed_fim.md)
  : Generic method for computing the observed FIM of an \`mle\` object.
- [`observed_fim(`*`<mle>`*`)`](https://queelius.github.io/algebraic.mle/reference/observed_fim.mle.md)
  : Function for obtaining the observed FIM of an \`mle\` object.
- [`orthogonal()`](https://queelius.github.io/algebraic.mle/reference/orthogonal.md)
  : Generic method for determining the orthogonal parameters of an
  estimator.
- [`orthogonal(`*`<mle>`*`)`](https://queelius.github.io/algebraic.mle/reference/orthogonal.mle.md)
  : Method for determining the orthogonal components of an \`mle\`
  object \`x\`.
- [`params(`*`<mle>`*`)`](https://queelius.github.io/algebraic.mle/reference/params.mle.md)
  : Method for obtaining the parameters of an \`mle\` object.
- [`params(`*`<mle_boot>`*`)`](https://queelius.github.io/algebraic.mle/reference/params.mle_boot.md)
  : Method for obtaining the parameters of an \`boot\` object.
- [`pred()`](https://queelius.github.io/algebraic.mle/reference/pred.md)
  : Generic method for computing the predictive confidence interval
  given an estimator object \`x\`.
- [`pred(`*`<mle>`*`)`](https://queelius.github.io/algebraic.mle/reference/pred.mle.md)
  : Estimate of predictive interval of \`T\|data\` using Monte Carlo
  integration.
- [`print(`*`<mle>`*`)`](https://queelius.github.io/algebraic.mle/reference/print.mle.md)
  : Print method for \`mle\` objects.
- [`print(`*`<summary_mle>`*`)`](https://queelius.github.io/algebraic.mle/reference/print.summary_mle.md)
  : Function for printing a \`summary\` object for an \`mle\` object.
- [`rmap(`*`<mle>`*`)`](https://queelius.github.io/algebraic.mle/reference/rmap.mle.md)
  : Computes the distribution of \`g(x)\` where \`x\` is an \`mle\`
  object.
- [`sampler(`*`<mle>`*`)`](https://queelius.github.io/algebraic.mle/reference/sampler.mle.md)
  : Method for sampling from an \`mle\` object.
- [`sampler(`*`<mle_boot>`*`)`](https://queelius.github.io/algebraic.mle/reference/sampler.mle_boot.md)
  : Method for sampling from an \`mle_boot\` object.
- [`score_val()`](https://queelius.github.io/algebraic.mle/reference/score_val.md)
  : Generic method for computing the score of an estimator object
  (gradient of its log-likelihood function evaluated at the MLE).
- [`score_val(`*`<mle>`*`)`](https://queelius.github.io/algebraic.mle/reference/score_val.mle.md)
  : Computes the score of an \`mle\` object (score evaluated at the
  MLE).
- [`se()`](https://queelius.github.io/algebraic.mle/reference/se.md) :
  Generic method for obtaining the standard errors of an estimator.
- [`se(`*`<mle>`*`)`](https://queelius.github.io/algebraic.mle/reference/se.mle.md)
  : Function for obtaining an estimate of the standard error of the MLE
  object \`x\`.
- [`summary(`*`<mle>`*`)`](https://queelius.github.io/algebraic.mle/reference/summary.mle.md)
  : Function for obtaining a summary of \`object\`, which is a fitted
  \`mle\` object.
- [`vcov(`*`<mle>`*`)`](https://queelius.github.io/algebraic.mle/reference/vcov.mle.md)
  : Computes the variance-covariance matrix of \`mle\` object.
- [`vcov(`*`<mle_boot>`*`)`](https://queelius.github.io/algebraic.mle/reference/vcov.mle_boot.md)
  : Computes the variance-covariance matrix of \`boot\` object. Note:
  This impelements the \`vcov\` method defined in \`stats\`.
