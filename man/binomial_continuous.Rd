% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/family.R
\name{binomial_continuous}
\alias{binomial_continuous}
\title{Family object for fitting binomial regression model with potentially unbounded continuous outcome}
\usage{
binomial_continuous(
  link = "logit",
  variance = "mu(1-mu)",
  family = c("gaussian", "binomial")
)
}
\arguments{
\item{link}{see \code{\link[stats:family]{family}}. Default to \code{"logit"}}

\item{variance}{see \code{\link[stats:family]{family}}. Default to \code{"mu(1-mu)"}, same as logistic regression.}

\item{family}{The family of the returned family object. Either \code{"gaussian"} or \code{"binomial"}. Default to \code{"gaussian"}. Seems not to matter for glm.}
}
\value{
a family object whose \code{family} is \code{"gaussian"}
}
\description{
A family object for fitting binomial models (i.e., generalized linear models with range contained in the open unit interval \eqn{(0,1)}) with continuous outcomes that may fall outside the unit interval \eqn{[0,1]}. Also works with \code{\link[glmnet:glmnet]{glmnet::glmnet}}.
}
\details{
This family is useful, for example, when the estimand is a conditional probability function while the outcome is a transformed psudo-outcome so that the estimator is multiply robust, or estimating a regression function with known bounds while the outcome might not respect the known bounds. Naive approaches such as \code{glm(family=binomial())}, \code{glm(family=quasibinomial())}, \code{glm(family=gaussian(link="logit"))}, \code{glm(family=quasi(link="logit",variance="constant"))} etc. might not work appropriately in such cases.

The output has \code{family="gaussian"} by default to be compatible with other learners in \code{\link[SuperLearner:SuperLearner]{SuperLearner::SuperLearner}}, because when the outcome is continuous, other learners might not perform correctly with \code{family="binomial"}.
}
\section{Warning}{


This function tweaks basic family objects and might remove some safety features. USE WITH CARE!!!
}

\examples{
set.seed(123)
expit <- binomial()$linkinv

# glm
x <- rnorm(100)
y <- expit(1 + x) + rnorm(100)
glm(y~x, family = binomial_continuous()) # or family=binomial_continuous, or family="binomial_continuous"
# Errors
try(glm(y~x, family = binomial()), outFile=stdout())
try(glm(y~x, family = quasibinomial()), outFile=stdout())
try(glm(y~x, family = gaussian(link = "logit")), outFile=stdout())
try(glm(y~x, family = quasi(link = "logit", variance = "constant")), outFile=stdout())

#glmnet
X <- matrix(rnorm(100 * 5), nrow = 100)
y <- expit(1 + X[,1]) + rnorm(100)
glmnet::glmnet(X, y, family = binomial_continuous()) # or family=binomial_continuous; cannot use family="binomial_continuous"
# Errors
try(glmnet::glmnet(X, y, family = binomial()), outFile=stdout())
try(glmnet::glmnet(X, y, family = gaussian(link = "logit")), outFile=stdout())
try(glmnet::glmnet(X, y, family = quasi(link = "logit", variance = "constant")), outFile=stdout())

# other links/variance for glm
x <- rnorm(100)
y <- expit(1 + x) + rnorm(100)
glm(y~x, family = binomial_continuous("probit")) # probit regression
glm(y~x, family = binomial_continuous(variance = "constant")) # least squares

# within SuperLearner
X <- matrix(rnorm(100 * 3), nrow = 100)
y <- expit(1 + X[,1]) + rnorm(10)
library(SuperLearner)
SuperLearner(y, data.frame(X), family=binomial_continuous(), SL.library = c("SL.glm", "SL.ipredbagg"), cvControl = list(V = 2))
# Error in SL.ipredbagg
SuperLearner(y, data.frame(X), family = binomial_continuous(family = "binomial"), SL.library = c("SL.glm", "SL.ipredbagg"), cvControl = list(V = 2))


}
