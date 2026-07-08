#' QuantFit: Latent Structure Model Selection Framework
#'
#' @description
#' Implements the model selection framework from Torres Irribarra & Diakow's paper
#' "Categorization, Ordering and Quantification: Selecting a Latent Variable Model
#' by Comparing Latent Structures."
#'
#' The package enables researchers to compare six latent structure models to determine
#' whether data supports classificatory, ordinal, or quantitative interpretations of
#' a latent variable.
#'
#' @section The Six Models:
#' \describe{
#'   \item{UN (Unconstrained)}{Standard latent class analysis with no ordering constraints.
#'     Represents a purely classificatory interpretation.}
#'   \item{MON (Class Monotonicity)}{Constrains item probabilities to be non-decreasing
#'     across classes. Implies classes can be ordered.}
#'   \item{IIO (Invariant Item Ordering)}{Constrains items to maintain the same relative
#'     difficulty ordering across all classes.}
#'   \item{DM (Double Monotonicity)}{Combines MON and IIO constraints. Implies a
#'     unidimensional ordinal scale similar to Mokken scaling.}
#'   \item{LCR (Latent Class Rasch)}{Uses Rasch parameterization with discrete ability
#'     classes. Implies interval-level measurement.}
#'   \item{RM (Rasch Model)}{Standard Rasch model with continuous latent trait.
#'     Represents the strongest quantitative interpretation.}
#' }
#'
#' @section Main Functions:
#' \describe{
#'   \item{\code{\link{fit_un}}}{Fit unconstrained latent class model}
#'   \item{\code{\link{fit_mon}}}{Fit class monotonicity model}
#'   \item{\code{\link{fit_iio}}}{Fit invariant item ordering model}
#'   \item{\code{\link{fit_dm}}}{Fit double monotonicity model}
#'   \item{\code{\link{fit_lcr}}}{Fit latent class Rasch model}
#'   \item{\code{\link{fit_rm}}}{Fit Rasch model (mirt wrapper)}
#'   \item{\code{\link{compare_models}}}{Compare all six models}
#'   \item{\code{\link{successive_comparison}}}{Stepwise model comparison}
#' }
#'
#' @section Typical Workflow:
#' \enumerate{
#'   \item Prepare binary response data (0/1 matrix or data frame)
#'   \item Use \code{compare_models()} for quick comparison of all models
#'   \item Or use \code{successive_comparison()} for guided interpretation
#'   \item Examine \code{plot_irfs()} to visualize item response functions
#'   \item Use \code{summary()} on the best-fitting model
#' }
#'
#' @section Data Requirements:
#' \itemize{
#'   \item Binary responses only (0/1)
#'   \item Complete data (no missing values by default)
#'   \item Sufficient sample size for stable estimation (N > 200 recommended)
#'   \item Multiple items (I >= 5 recommended)
#' }
#'
#' @references
#' Torres Irribarra, D., & Diakow, R. Categorization, Ordering and Quantification:
#' Selecting a Latent Variable Model by Comparing Latent Structures.
#'
#' @examples
#' \dontrun{
#' # Generate example data
#' set.seed(123)
#' n <- 500
#' n_items <- 10
#'
#' # Data with ordinal structure
#' theta <- rnorm(n)
#' delta <- seq(-1.5, 1.5, length.out = n_items)
#' data <- matrix(0, n, n_items)
#' for (i in 1:n) {
#'   for (j in 1:n_items) {
#'     prob <- plogis(theta[i] - delta[j])
#'     data[i, j] <- rbinom(1, 1, prob)
#'   }
#' }
#'
#' # Compare all models
#' comparison <- compare_models(data, n_classes = 4)
#' print(comparison)
#'
#' # Guided comparison
#' result <- successive_comparison(data, n_classes = 4)
#' print(result$conclusion)
#'
#' # Visualize best model
#' best_fit <- get_model(comparison, comparison$best_model)
#' plot_irfs(best_fit)
#' }
#'
#' @docType package
#' @name QuantFit-package
#' @aliases QuantFit
"_PACKAGE"

# Package imports
#' @import stats
#' @import utils
#' @importFrom alabama constrOptim.nl
#' @importFrom nloptr nloptr
#' @importFrom Matrix Matrix
#' @importFrom numDeriv hessian
NULL

# Package startup message
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "QuantFit: Latent Structure Model Selection Framework\n",
    "For help, type: ?QuantFit or vignette('introduction')"
  )
}
