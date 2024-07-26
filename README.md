# README

This package contains a function that outputs a family object for fitting a binomial model when the outcome is continuous.

To install:
```{r}
devtools::install_github("Qiu-Hongxiang-David/BinomCont")
```

Example:
```{r}
library(BinomCont)
x <- rnorm(100)
y <- 1/(1 + exp(-x)) + rnorm(100)

glm(y~x, family = binomial_continuous())
```

More naive approaches appear not working.
