
# General version for general statistics
#' Title
#'
#' @param Z A length N vector of (nonnegative) observed doses.
#' @param Q A length N vector of observed binary outcomes.
#' @param index A length N vector of indices indicating matched set membership.
#' @param mc An integer for the total number of Monte-Carlo samples.
#' @param gamma The nonnegative sensitivity parameter; gamma = 0 means no
#' unmeasured confounding.
#' @param trans The transformation of the doses to use for the test statistic.
#' Default is the identity function.
#' @param direct The direction of the test - "upper" or "lower"; default is upper.
#'
#' @return A length mc vector consisting of the Monte-Carlo samples of the
#' test statistic.
#' @import gtools
#' @export
#'
#' @examples
dose_sensitivity_mc_gen <- function(Z, Q, index, mc, gamma,
                                    trans = identity, direct = "upper") {
  if(any(Z) < 0) {
    stop("Negative doses")
  }
  if(all(Q %in% c(0,1))) {
    stop("Non-binary outcomes")
  }
  if (gamma) {
    stop("Negative gamma")
  }
  if (!(direct %in% c("upper","lower"))) {
    stop("invalid direct")
  }

  # Matched set indices
  match_index = unique(index)
  # Greater than obsT indicator
  indicator <- rep(NA, mc)
  # number of matched sets
  nostratum <- length(unique(index))
  # worst case distribution in each matched set
  worst_case <- vector('list',nostratum)
  for (j in 1:nostratum) {
    doses <- Z[which(index == match_index[j])]
    resp <- Q[which(index == match_index[j])]
    ns <- length(doses)
    permZ <- matrix(doses[permutations(ns,ns)],nrow=factorial(ns))
    zTr <- (trans(permZ)) %*% resp
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
  return(mc_samps)

}


# Normal approximation, general transformation
#' Title
#'
#' @param Z A length N vector of (nonnegative) observed doses.
#' @param Q A length N vector of observed binary outcomes.
#' @param index A length N vector of indices indicating matched set membership.
#' @param gamma The nonnegative sensitivity parameter; gamma = 0 means no
#' unmeasured confounding.
#' @param trans The transformation of the doses to use for the test statistic.
#' Default is the identity function.
#' @param obsT The observed value of the test statistic; default is NULL
#' @param direct The direction of the test - "upper" or "lower"; default is upper.
#'
#' @return The normal approximation deviate.
#' @import gtools
#' @export
#'
#' @examples
normal_test_gen <- function(Z, Q, index, gamma, trans = identity,
                            obsT = NULL, direct = "upper") {
  if(any(Z) < 0) {
    stop("Negative doses")
  }
  if(all(Q %in% c(0,1))) {
    stop("Non-binary outcomes")
  }
  if (gamma) {
    stop("Negative gamma")
  }
  if (!(direct %in% c("upper","lower"))) {
    stop("invalid direct")
  }

  if(is.null(obsT)) {
    obsT <- sum(trans(Z) * Q)
  }

  # Matched set indices
  match_index = unique(index)

  # number of matched sets
  nostratum <- length(unique(index))

  variance = 0
  expectation = 0
  for (j in 1:nostratum) {
    doses <- Z[which(index == match_index[j])]
    resp <- Q[which(index == match_index[j])]
    ns <- length(doses)
    permZ <- matrix(doses[permutations(ns,ns)],nrow=factorial(ns))
    zTr <- trans(permZ) %*% resp

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
    return (deviate)
  } else if (direct == "lower") {
    deviate = (obsT - expectation) / sqrt(variance)
    return (deviate)
  }


}

# Deviation from uniform function in TV distance
#' Title
#'
#' @param Z A length N vector of (nonnegative) observed doses.
#' @param Q A length N vector of observed binary outcomes.
#' @param index A length N vector of indices indicating matched set membership.
#' @param gamma The nonnegative sensitivity parameter; gamma = 0 means no
#' unmeasured confounding.
#' @param direct The direction of the test - "upper" or "lower"; default is upper.
#'
#' @return A vector of length equaling the number of matched sets consisting
#' of the TV distance from the uniform for gamma level of unmeasured confounding.
#' @import gtools
#' @export
#'
#' @examples
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
