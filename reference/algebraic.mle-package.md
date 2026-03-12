# algebraic.mle: Algebraic Maximum Likelihood Estimators

The maximum likelihood estimator (MLE) is a technology: under regularity
conditions, any MLE is asymptotically normal with variance given by the
inverse Fisher information. This package exploits that structure by
defining an algebra over MLEs. Compose independent estimators into joint
MLEs via block-diagonal covariance ('joint'), optimally combine repeated
estimates via inverse-variance weighting ('combine'), propagate
transformations via the delta method ('rmap'), and bridge to
distribution algebra via conversion to normal or multivariate normal
objects ('as_dist'). Supports asymptotic ('mle', 'mle_numerical') and
bootstrap ('mle_boot') estimators with a unified interface for
inference: confidence intervals, standard errors, AIC, Fisher
information, and predictive intervals. For background on maximum
likelihood estimation, see Casella and Berger (2002,
ISBN:978-0534243128). For the delta method and variance estimation, see
Lehmann and Casella (1998, ISBN:978-0387985022).

## See also

Useful links:

- <https://github.com/queelius/algebraic.mle>

- <https://queelius.github.io/algebraic.mle/>

- Report bugs at <https://github.com/queelius/algebraic.mle/issues>

## Author

**Maintainer**: Alexander Towell <lex@metafunctor.com>
([ORCID](https://orcid.org/0000-0001-6443-9897))
