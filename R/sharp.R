
#' Sharp null monte-carlo sensitivity analysis for continuous exposures
#' and binary outcomes.
#'
#' @param Z A length N vector of (nonnegative) observed doses.
#' @param Q A length N vector of observed binary outcomes.
#' @param index A length N vector of indices indicating matched set membership.
#' @param mc An integer for the total number of Monte-Carlo samples desired.
#' @param gamma The nonnegative sensitivity parameter; gamma = 0 means no
#' unmeasured confounding.
#' @param trans The transformation of the doses to use for the test statistic.
#' Default is the identity function.
#' @param direct The direction of the test - "upper" or "lower"; default is upper.
#' @param weights Weights for each stratum to apply for the test statistic
#' @param seed seed for random number generation.
#' @param obsT The observed value of the test statistic; default is NULL
#'
#' @return A list containing two objects:
#'
#' \item{mc}{A length mc vector containing the monte-carlo samples of the test
#' statistic.}
#' \item{p}{The monte-carlo p-value.}
#' @import gtools
#' @export
#'
#' @examples
#' # Load the data
#' data <- treat_out_match
#' # Make a threshold at log(3.5) transformation function.
#' above = function(Z) { return(Z > log(3.5)) }
#' # Conduct randomization test.
#' solution <- dose_sensitivity_mc_gen(data$treat, data$complain, data$match_ind,
#' mc = 250, gamma = 0, trans = above)
#'
dose_sensitivity_mc_gen <- function(Z, Q, index, mc, gamma, weights = NA, obsT = NULL,
                                    trans = identity, direct = "upper", seed = 1) {
  # set the seed.
  set.seed(seed)

  # error-checking
  if(!all(Q %in% c(0,1))) {
    stop("Non-binary outcomes")
  }
  if (gamma < 0) {
    stop("Negative gamma")
  }
  if (!(direct %in% c("upper","lower"))) {
    stop("invalid direction")
  }

  # Matched set indices
  match_index = unique(index)
  # number of matched sets
  nostratum <- length(unique(index))
  # initilize weights if null.
  if (any(is.na(weights))) {
    weights <- rep(1, nostratum)
  }
  # compute test statistic.
  if(is.null(obsT)) {
    individual_weights <- rep(weights, times = as.vector(table(index)))
    obsT <- sum(individual_weights * trans(Z) * Q)
  }
  # worst case distribution in each matched set
  worst_case <- vector('list',nostratum)
  for (j in 1:nostratum) {
    doses <- Z[which(index == match_index[j])]
    resp <- Q[which(index == match_index[j])]
    ns <- length(doses)
    permZ <- matrix(doses[permutations(ns,ns)],nrow=factorial(ns))
    zTr <- (trans(permZ)) %*% resp * weights[j]
    if (gamma == 0) {
      probs <- rep(1/nrow(permZ),nrow(permZ))
    } else {
      if (direct == "upper") {
        probs <- exp(gamma*(permZ%*% resp))
      } else if (direct == "lower") {
        probs <- exp(gamma*(permZ%*% (1-resp)))
      }
      sum_probs <- sum(probs)
      probs <- probs/sum_probs
    }


    worst_case[[j]] <- data.frame(probs = probs, zTr = zTr)

  }
  print("Initialized")

  mc_samps = rep(NA,mc)
  for (l in 1:mc) {
    if (l %% 1000 == 0) {
      print(l)
    }
    samp = 0
    for (j in 1:nostratum) {
      samp <- samp + sample(worst_case[[j]]$zTr, 1, prob=worst_case[[j]]$probs, replace=TRUE)
    }

    mc_samps[l] <- samp
  }

  # monte carlo p-value
  return(list(mc = mc_samps, p = (1 + sum(mc_samps >= obsT))/ (1 + mc)))

}


#' Sharp null sensitivity analysis for continuous exposures and binary outcomes
#' using normal approximation.
#'
#' @param Z A length N vector of (nonnegative) observed doses.
#' @param Q A length N vector of observed binary outcomes.
#' @param index A length N vector of indices indicating matched set membership.
#' @param gamma The nonnegative sensitivity parameter; gamma = 0 means no
#' unmeasured confounding.
#' @param trans The transformation of the doses to use for the test statistic.
#' Default is the identity function.
#' @param obsT The observed value of the test statistic; default is NULL.
#' @param direct The direction of the test - "upper" or "lower"; default is upper.
#' @param weights Weights to apply for the test statistic
#'
#' @return A list containing the following:
#'
#' \item{obsT}{The observed value of the test statistic}
#' \item{exp}{The worst-case expectation}
#' \item{var}{The worst-case variance.}
#' \item{deviate}{The normal approximation deviate.}
#' @import gtools
#' @export
#'
#' @examples
#' # Load the data
#' data <- treat_out_match
#' # Make a threshold at log(3.5) transformation function.
#' above = function(Z) { return(Z > log(3.5)) }
#' # Conduct randomization test using normal approximation.
#' solution <- normal_test_gen(data$treat, data$complain, data$match_ind,
#' gamma = 0, trans = above)
normal_test_gen <- function(Z, Q, index, gamma, trans = identity, weights = NA,
                            obsT = NULL, direct = "upper") {
  if(!all(Q %in% c(0,1))) {
    stop("Non-binary outcomes")
  }
  if (gamma < 0) {
    stop("Negative gamma")
  }
  if (!(direct %in% c("upper","lower"))) {
    stop("invalid direction")
  }

  # Matched set indices
  match_index = unique(index)

  # number of matched sets
  nostratum <- length(unique(index))

  # initialize weights if null
  if (any(is.na(weights))) {
    weights <- rep(1, nostratum)
  }
  # compute test statistic.
  if(is.null(obsT)) {
    individual_weights <- rep(weights, times = as.vector(table(index)))
    obsT <- sum(individual_weights * trans(Z) * Q)
  }

  # compute worst-case expectation and variance.
  variance = 0
  expectation = 0
  for (j in 1:nostratum) {
    doses <- Z[which(index == match_index[j])]
    resp <- Q[which(index == match_index[j])]
    ns <- length(doses)
    permZ <- matrix(doses[permutations(ns,ns)],nrow=factorial(ns))
    zTr <- trans(permZ) %*% resp * weights[j]

    if (gamma == 0) {
      probs <- rep(1/nrow(permZ),nrow(permZ))
    } else {
      if (direct == "upper") {
        probs <- exp(gamma*(permZ%*% resp))
      } else if (direct == "lower") {
        probs <- exp(gamma*(permZ%*% (1-resp)))
      }
      sum_probs <- sum(probs)
      probs <- probs/sum_probs
    }
    match_mean = sum(probs * zTr)
    expectation = expectation + match_mean
    variance = variance + sum(probs * zTr^2) - (match_mean)^2

  }


  # Check normal deviate
  if (direct == "upper") {
    deviate = (obsT - expectation) / sqrt(variance)
    return (list(obsT = obsT, exp = expectation, var = variance,
                 deviate = deviate))
  } else if (direct == "lower") {
    deviate = (obsT - expectation) / sqrt(variance)
    return (list(obsT = obsT, exp = expectation, var = variance,
                 deviate = deviate))
  }


}

#' Computes deviation from uniform distribution in total variation
#' distance for a given amount of unmeasured confounding and a greater than
#' alternative with a binary outcome.
#'
#' @param Z A length N vector of (nonnegative) observed doses.
#' @param Q A length N vector of observed binary outcomes.
#' @param index A length N vector of indices indicating matched set membership.
#' @param gamma The nonnegative sensitivity parameter; gamma = 0 means no
#' unmeasured confounding.
#' @param direct The direction of the test - "upper" or "lower"; default is upper.
#'
#' @return A vector of length equaling the number of matched sets consisting
#' of the TV distance from the uniform for each matched set at gamma level of
#' unmeasured confounding for the worst-case.
#' @import gtools
#' @export
#'
#' @examples
#' # Load the data
#' data <- treat_out_match
#' # compute total variation distances.
#' total_variation <- dev_TV(data$treat, data$complain,
#' data$match_ind, gamma = log(1.5))
dev_TV <- function(Z, Q, index, gamma, direct = "upper") {
  # Matched set indices
  match_index = unique(index)

  # number of matched sets
  nostratum <- length(unique(index))

  # Vector to record deviations from uniform
  deviations <- rep(0,nostratum)

  if (gamma == 0) {
    return(deviations)
  }
  for (j in 1:nostratum) {

    doses <- Z[which(index == match_index[j])]
    resp <- Q[which(index == match_index[j])]
    ns <- length(doses)
    permZ <- matrix(doses[permutations(ns,ns)],nrow=factorial(ns))

    if (direct == "upper") {
      probs <- exp(gamma*(permZ%*% resp))
    } else if (direct == "lower") {
      probs <- exp(gamma*(permZ%*% (1-resp)))
    }
    sum_probs <- sum(probs)
    probs <- probs/sum_probs
    deviations[j] <- 0.5 * sum(abs(probs-1/factorial(ns)))

  }

  return(deviations)
}
