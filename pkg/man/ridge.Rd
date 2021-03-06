\name{ridge}
\alias{ridge}
\title{Fit a linear model with a structured ridge regularization}
\usage{
  ridge(x, y, lambda2 = NULL, struct = NULL,
    intercept = TRUE, normalize = TRUE, nlambda2 = 100,
    lambda.min = ifelse(n <= p, 0.01, 1e-04),
    lambda.max = 100, control = list(), checkargs = TRUE)
}
\arguments{
  \item{x}{matrix of features. Do NOT include intercept.
  When normalized os \code{TRUE}, coefficients will then be
  rescaled to the original scale.}

  \item{y}{response vector.}

  \item{lambda2}{sequence of decreasing
  \eqn{\ell_2}{l2}-penalty levels. If \code{NULL} (the
  default), a vector is generated with \code{nlambda2}
  entries.}

  \item{struct}{matrix structuring the coefficients,
  possibly sparsely encoded. Must be at least positive
  semidefinite (this is checked internally if the
  \code{checkarg} argument is \code{TRUE}). If \code{NULL}
  (the default), the identity matrix is used. See details
  below.}

  \item{intercept}{logical; indicates if an intercept
  should be included in the model. Default is \code{TRUE}.}

  \item{normalize}{logical; indicates if variables should
  be normalized to have unit L2 norm before fitting.
  Default is \code{TRUE}.}

  \item{nlambda2}{integer that indicates the number of
  values to put in the \code{lambda2} vector.  Ignored if
  \code{lambda2} is provided.}

  \item{lambda.min}{the minimal amount of penalty used to
  generated the vector \code{lambda2}. Ignored if
  \code{lambda2} is provided.}

  \item{lambda.max}{the maximal amount of penalty used to
  generated the vector \code{lambda2}. Ignored if
  \code{lambda2} is provided.}

  \item{control}{list of argument controlling low level
  options of the algorithm --use with care and at your own
  risk-- : \itemize{%

  \item{\code{verbose}: }{integer; activate verbose mode
  --this one is not too much risky!-- set to \code{0} for
  no output; \code{1} for warnings only, and \code{2} for
  tracing the whole progression. Default is \code{1}.
  Automatically set to \code{0} when the method is embedded
  within cross-validation or stability selection.}

  \item{\code{timer}: }{logical; use to record the timing
  of the algorithm. Default is \code{FALSE}.}

  }}

  \item{checkargs}{logical; should arguments be checked to
  (hopefully) avoid internal crashes? Default is
  \code{TRUE}. Automatically set to \code{FALSE} when calls
  are made from cross-validation or stability selection
  procedures.}
}
\value{
  an object with class \code{quadrupen}, see the
  documentation page \code{\linkS4class{quadrupen}} for
  details.
}
\description{
  Adjust a linear model with ridge regularization (possibly
  structured \eqn{\ell_2}{l2}-norm). The solution path is
  computed at a grid of values for the
  \eqn{\ell_2}{l2}-penalty. See details for the criterion
  optimized.
}
\note{
  The optimized criterion is the following:
  \if{latex}{\deqn{% \hat{\beta}_{\lambda_2} = \arg
  \min_{\beta} \frac{1}{2} (y - X \beta)^T (y - X \beta) +
  \frac{\lambda_2}{2} \beta^T S \beta, }} \if{html}{\out{
  <center> &beta;<sup>hat</sup>
  <sub>&lambda;<sub>2</sub></sub> = argmin<sub>&beta;</sub>
  1/2 RSS(&beta) + + &lambda;/2 <sub>2</sub>
  &beta;<sup>T</sup> S &beta;, </center> }}
  \if{text}{\deqn{beta.hat(lambda2) = argmin_beta 1/2
  RSS(beta) + lambda2 beta' S beta,}} where the
  \eqn{\ell_2}{l2} structuring positive semidefinite matrix
  \eqn{S}{S} is provided via the \code{struct} argument
  (possibly of class \code{Matrix}).
}
\seealso{
  See also \code{\linkS4class{quadrupen}},
  \code{\link{plot.quadrupen}} and \code{\link{crossval}}.
}
\keyword{models,}
\keyword{regression}

