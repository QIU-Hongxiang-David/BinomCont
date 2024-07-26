clip_interval<-function(x,lower=-Inf,upper=Inf){
    out<-pmax(pmin(x,upper),lower)
    if(is.matrix(x)){
        out<-matrix(out,nrow=nrow(x),ncol=ncol(x))
    }
    out
}


#' @title Family object for fitting binomial regression model with potentially unbounded continuous outcome
#' @name binomial_continuous
#' @description
#' A family object for fitting binomial models (i.e., generalized linear models with range contained in the open unit interval \eqn{(0,1)}) with continuous outcomes that may fall outside the unit interval \eqn{[0,1]}. Also works with \code{\link[glmnet:glmnet]{glmnet::glmnet}}.
#' @details
#' This family is useful, for example, when the estimand is a conditional probability function while the outcome is a transformed psudo-outcome so that the estimator is multiply robust, or estimating a regression function with known bounds while the outcome might not respect the known bounds. Naive approaches such as \code{glm(family=binomial())}, \code{glm(family=quasibinomial())}, \code{glm(family=gaussian(link="logit"))}, \code{glm(family=quasi(link="logit",variance="constant"))} etc. might not work appropriately in such cases.
#'
#' The output has \code{family="gaussian"} by default to be compatible with other learners in \code{\link[SuperLearner:SuperLearner]{SuperLearner::SuperLearner}}, because when the outcome is continuous, other learners might not perform correctly with \code{family="binomial"}.
#'
#' @section Warning:
#'
#' This function tweaks basic family objects and might remove some safety features. USE WITH CARE!!!
#' @param link see \code{\link[stats:family]{family}}. Default to \code{"logit"}
#' @param variance see \code{\link[stats:family]{family}}. Default to \code{"mu(1-mu)"}, same as logistic regression.
#' @param family The family of the returned family object. Either \code{"gaussian"} or \code{"binomial"}. Default to \code{"gaussian"}. Seems not to matter for glm.
#' @returns a family object whose \code{family} is \code{"gaussian"}
#'
#' @examples
#' set.seed(123)
#' expit <- binomial()$linkinv
#'
#' # glm
#' x <- rnorm(100)
#' y <- expit(1 + x) + rnorm(100)
#' glm(y~x, family = binomial_continuous()) # or family=binomial_continuous, or family="binomial_continuous"
#' # Errors
#' try(glm(y~x, family = binomial()), outFile=stdout())
#' try(glm(y~x, family = quasibinomial()), outFile=stdout())
#' try(glm(y~x, family = gaussian(link = "logit")), outFile=stdout())
#' try(glm(y~x, family = quasi(link = "logit", variance = "constant")), outFile=stdout())
#'
#' #glmnet
#' X <- matrix(rnorm(100 * 5), nrow = 100)
#' y <- expit(1 + X[,1]) + rnorm(100)
#' glmnet::glmnet(X, y, family = binomial_continuous()) # or family=binomial_continuous; cannot use family="binomial_continuous"
#' # Errors
#' try(glmnet::glmnet(X, y, family = binomial()), outFile=stdout())
#' try(glmnet::glmnet(X, y, family = gaussian(link = "logit")), outFile=stdout())
#' try(glmnet::glmnet(X, y, family = quasi(link = "logit", variance = "constant")), outFile=stdout())
#'
#' # other links/variance for glm
#' x <- rnorm(100)
#' y <- expit(1 + x) + rnorm(100)
#' glm(y~x, family = binomial_continuous("probit")) # probit regression
#' glm(y~x, family = binomial_continuous(variance = "constant")) # least squares
#'
#' # within SuperLearner
#' X <- matrix(rnorm(100 * 3), nrow = 100)
#' y <- expit(1 + X[,1]) + rnorm(10)
#' library(SuperLearner)
#' SuperLearner(y, data.frame(X), family=binomial_continuous(), SL.library = c("SL.glm", "SL.ipredbagg"), cvControl = list(V = 2))
#' # Error in SL.ipredbagg
#' SuperLearner(y, data.frame(X), family = binomial_continuous(family = "binomial"), SL.library = c("SL.glm", "SL.ipredbagg"), cvControl = list(V = 2))
#'
#'
#' @export
binomial_continuous<-function(link="logit",variance="mu(1-mu)",family=c("gaussian","binomial")){
    linktemp <- substitute(link)
    if (!is.character(linktemp))
        linktemp <- deparse(linktemp)
    if (linktemp %in% c("logit", "probit", "cloglog", "identity", "inverse", "log", "1/mu^2", "sqrt")){
        stats <- make.link(linktemp)
    }else if (is.character(link)) {
        stats <- make.link(link)
        linktemp <- link
    } else {
        stats <- link
        linktemp <- stats$name %||% deparse(linktemp)
    }

    family<-match.arg(family)

    variancetemp <- substitute(variance)
    if (!is.character(variancetemp))
        variancetemp <- deparse(variancetemp)

    eps<-.Machine$double.eps

    out<-do.call(quasi,list(link=linktemp,variance=variancetemp))
    out$family<-family
    out$linkfun<-function(mu){
        .Call(stats:::C_logit_link, clip_interval(mu,eps,1-eps))
    }
    out$linkinv<-function(eta){
        clip_interval(.Call(stats:::C_logit_linkinv, eta),eps,1-eps)
    }
    out$initialize<-expression({
        n <- rep.int(1, nobs)
        mustart <- pmax(0.001, pmin(0.999, y))
    })
    out$validmu<-binomial()$validmu
    if(is.character(variancetemp) && variancetemp=="mu(1-mu)"){
        out$dev.resids<-function(y,mu,wt){
            mu<-clip_interval(mu,eps,1-eps)
            -2*wt*(y*log(mu)+(1-y)*log1p(-mu))
        }
    }
    out
}
