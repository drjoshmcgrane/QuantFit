#' Simulate item responses from a specified latent-structure model
#'
#' Generates dichotomous or polytomous item response data from any of the six
#' models in the Torres Irribarra & Diakow hierarchy, for validation and power
#' studies. The dichotomous case reproduces the generators of the original paper
#' (random logits on \eqn{U(-4, 4)}, sorted across classes for class
#' monotonicity, across items for invariant item ordering, or both for double
#' monotonicity; a partial-credit/Rasch structure for the quantitative models),
#' and the polytomous case generalises them through a cumulative-threshold
#' (for UN/MON/IIO/DM) or partial-credit (for LCR/RM) formulation.
#'
#' @details
#' For the non-parametric models each item-by-class combination is given a
#' location \eqn{L_{cj}} and the category structure is built from shared,
#' increasing threshold offsets, so that
#' \eqn{P(X \ge k \mid c, j) = \mathrm{logit}^{-1}(L_{cj} - \tau_k)}:
#' \describe{
#'   \item{UN}{\eqn{L} left unsorted (class and item profiles may cross).}
#'   \item{MON}{each column of \eqn{L} sorted, so category distributions are
#'     stochastically ordered across classes.}
#'   \item{IIO}{each row of \eqn{L} sorted, so expected item scores share one
#'     ordering across classes.}
#'   \item{DM}{both, so the data are doubly monotone but not necessarily
#'     partial-credit.}
#' }
#' The quantitative models use the partial credit model
#' \eqn{P(X = x \mid \theta, j) \propto \exp\{\sum_{l \le x}(\theta - \delta_{jl})\}}
#' with discrete, separated class locations (LCR) or a continuous
#' \eqn{\theta \sim N(0, 1)} (RM). When `n_cat = 2` every model reduces to the
#' dichotomous generator.
#'
#' @param model One of `"UN"`, `"MON"`, `"IIO"`, `"DM"`, `"LCR"`, `"RM"`.
#' @param n_persons Number of respondents.
#' @param n_items Number of items.
#' @param n_classes Number of latent classes (ignored for `"RM"`).
#' @param n_cat Number of ordered response categories (2 = dichotomous).
#' @param class_probs Optional class mixing proportions (length `n_classes`);
#'   defaults to equal.
#' @param seed Optional random seed.
#'
#' @return An integer matrix of responses (`n_persons` rows, `n_items` columns),
#'   scored `0..n_cat-1`, with attributes `"model"` and `"params"` (the
#'   generating parameters and, for the mixture models, the class memberships).
#'
#' @examples
#' # doubly monotone dichotomous data, 3 classes
#' d <- simulate_responses("DM", n_persons = 300, n_items = 8, n_classes = 3,
#'                         seed = 1)
#' # partial-credit (Rasch) data with 4 categories
#' p <- simulate_responses("RM", n_persons = 300, n_items = 6, n_cat = 4,
#'                         seed = 1)
#'
#' @seealso [select_model_ll()], [quant_fit()]
#' @export
simulate_responses <- function(model = c("UN", "MON", "IIO", "DM", "LCR", "RM"),
                               n_persons = 500, n_items = 10, n_classes = 3,
                               n_cat = 2L, class_probs = NULL, seed = NULL) {
  model <- match.arg(model)
  if (!is.null(seed)) set.seed(seed)
  n_cat <- as.integer(n_cat)
  if (n_cat < 2L) stop("n_cat must be >= 2")
  m <- n_cat - 1L
  tau <- if (m == 1L) 0 else seq(-1.2, 1.2, length.out = m)  # increasing offsets

  # cumulative P(X>=k) (length m) -> category probs (length m+1)
  cum_to_cat <- function(cum) c(1, cum) - c(cum, 0)

  params <- list()

  if (model %in% c("UN", "MON", "IIO", "DM")) {
    if (is.null(class_probs)) class_probs <- rep(1 / n_classes, n_classes)
    L <- matrix(stats::runif(n_classes * n_items, -4, 4), n_classes, n_items)
    if (model %in% c("MON", "DM")) L <- apply(L, 2L, sort)          # classes ordered
    if (model %in% c("IIO", "DM")) L <- t(apply(L, 1L, sort))       # items ordered
    cls <- sample.int(n_classes, n_persons, replace = TRUE, prob = class_probs)
    resp <- matrix(0L, n_persons, n_items)
    for (j in seq_len(n_items)) {
      cumj <- plogis(outer(L[, j], tau, "-"))                       # C x m: P(X>=k|c)
      P <- cbind(1, cumj) - cbind(cumj, 0)                          # C x (m+1)
      for (c in seq_len(n_classes)) {
        who <- which(cls == c)
        if (length(who))
          resp[who, j] <- sample.int(m + 1L, length(who), replace = TRUE,
                                     prob = P[c, ]) - 1L
      }
    }
    params <- list(L = L, tau = tau, class = cls, class_probs = class_probs)

  } else {  # LCR, RM: partial credit model
    b <- stats::runif(n_items, -2, 2)                              # item locations
    delta_list <- lapply(b, function(bj) bj + tau)                 # item step params
    if (model == "LCR") {
      if (is.null(class_probs)) class_probs <- rep(1 / n_classes, n_classes)
      repeat {                                                     # separated classes
        a <- sort(stats::runif(n_classes, -3, 3))
        if (all(diff(a) > 0.5)) break
      }
      cls <- sample.int(n_classes, n_persons, replace = TRUE, prob = class_probs)
      theta <- a[cls]
      params <- list(theta_class = a, item = b, tau = tau, class = cls)
    } else {
      theta <- stats::rnorm(n_persons)
      params <- list(item = b, tau = tau)
    }
    ip <- compute_pcm_probs(theta, delta_list, use_cpp = TRUE)      # list of n x (m+1)
    resp <- matrix(0L, n_persons, n_items)
    for (j in seq_len(n_items)) {
      cdf <- t(apply(ip[[j]], 1L, cumsum))
      resp[, j] <- rowSums(stats::runif(n_persons) > cdf)
    }
  }

  colnames(resp) <- paste0("Item", seq_len(n_items))
  attr(resp, "model") <- model
  attr(resp, "params") <- params
  resp
}
