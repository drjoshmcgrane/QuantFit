#' 1000 sampled 3-matrices from simulated Rasch data
#'
#' A pre-computed [ConjointChecks()] result: Rasch data were simulated
#' (2000 respondents, 20 items, normal abilities and difficulties) and the
#' double-cancellation check was run on 1000 randomly sampled 3x3
#' submatrices. Useful for exploring the `checks` class, [plot][plot.checks]
#' and [summary][summary.checks] methods without waiting for MCMC.
#'
#' @format An object of S4 class [`checks`][checks-class].
#' @source Simulated via the Rasch model; see the example code in the
#'   original ConjointChecks package documentation (seed 8675309).
#' @examples
#' data(rasch1000)
#' summary(rasch1000)
"rasch1000"
