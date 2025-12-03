# Estimate of predictive interval of \`T\|data\` using Monte Carlo integration.

Let \`T\|x ~ f(t\|x)“ be the pdf of vector \`T\` given MLE \`x\` and \`x
~ MVN(params(x),vcov(x))“ be the estimate of the sampling distribution
of the MLE for the parameters of \`T\`. Then, \`(T,x) ~ f(t,x) = f(t\|x)
f(x) is the joint distribution of \`(T,x)\`. To find \`f(t)\` for a
fixed \`t\`, we integrate \`f(t,x)\` over \`x\` using Monte Carlo
integration to find the marginal distribution of \`T\`. That is, we:

## Usage

``` r
# S3 method for class 'mle'
pred(x, samp, alpha = 0.05, R = 50000, ...)
```

## Arguments

- x:

  an \`mle\` object.

- samp:

  The sampler for the distribution that is parameterized by the MLE
  \`x\`, i.e., \`T\|x\`.

- alpha:

  (1-alpha)-predictive interval for \`T\|x\`. Defaults to 0.05.

- R:

  number of samples to draw from the sampling distribution of \`x\`.
  Defaults to 50000.

- ...:

  additional arguments to pass into \`samp\`.

## Details

1\. Sample from MVN \`x\` 2. Compute \`f(t,x)\` for each sample 3. Take
the mean of the \`f(t,x)\` values asn an estimate of \`f(t)\`.

The \`samp\` function is used to sample from the distribution of
\`T\|x\`. It should be designed to take
