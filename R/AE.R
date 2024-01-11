#' Separable algorithm for threshold attributable effect in a sensitivity
#' analysis with at most one over-exposed unit in each matched set.
#' For a greater than alternative, finds the 'a' matched sets that
#' most decrease the mean and/or variance.
#'
#' @param Z A length N vector of (nonnegative) observed doses.
#' @param Q A length N vector of observed binary outcomes.
#' @param index A length N vector of indices indicating matched set membership.
#' @param gamma The nonnegative sensitivity parameter; gamma = 0 means no
#' unmeasured confounding.
#' @param thresh The dose threshold for the TAE.
#' @param a The number of attributable effects to test for.
#' @param trans The transformation of the doses to use for the test statistic.
#' Default is the identity function.
#' @param mc Number of monte-carlo samples if testing the sharp null, i.e. a = 0.
#'
#' @return Either "reject" if the value a is deemed
#' not plausible/compatible, "feasible" if the value a is deemed so,
#'  else a list containing a p-value and dataframe of
#' matched sets that have contribution to the test statistic
#' sorted in order of smallest mean reduction followed by
#' smallest variance reduction.
#' @import dplyr
#' @export
#'
#' @examples
#' # Load the data
#' data <- treat_out_match
#' # Solve by the separable algorithm
#' solution <- binary_thresh_attribute(data$treat, data$complain, data$match_ind,
#' gamma = 0, thresh = log(3.5), a = 5, trans = identity)
#'
binary_thresh_attribute <- function(Z, Q, index, gamma, thresh = 0, a = 1, trans = identity,
                                   mc = 50000) {
  if(!all(Q %in% c(0,1))) {
    stop("Non-binary outcomes")
  }
  if (gamma < 0) {
    stop("Negative gamma")
  }
  # The value of the test statistic under observed dose and response
  Zthresh = trans(Z) > thresh
  TT = sum(Zthresh*Q)
  # Reject if a too large.
  if (a >= TT) {
    return("Reject")
  }
  if (a == 0) {
    return(list(result=dose_sensitivity_mc_gen(Z, Q, index, mc = mc, gamma,
                                               trans = function(x){trans(x)>thresh},
                                               direct = "upper")))
  }
  # Find the indices where the contribution from the matched set is 1.
  data  = data.frame(Z = Z, Q = Q, index = index)
  contribute = data %>% group_by(index) %>%
    summarise(contrib = sum((trans(Z) > thresh)*Q))
  contribute_ind = contribute$index[which(contribute$contrib == 1)]

  # Matched set indices
  match_index = unique(index)

  # number of matched sets
  nostratum <- length(unique(index))

  # Means and variances if attribute or not.
  lambda_1 = rep(0,nostratum)
  lambda_2 = rep(0,nostratum)
  omega_1 = rep(0,nostratum)
  omega_2 = rep(0,nostratum)

  # Initialize for each matched set.
  for (j in 1:nostratum) {
    # doses in set j.
    doses <- Z[which(index == match_index[j])]
    # responses in set j.
    resp <- Q[which(index == match_index[j])]
    # find the index of subject that is potentially attributable.
    attr_ind <- which((trans(doses) > thresh)*resp == 1)
    # response under no dose if attributable
    attr_resp <- resp
    attr_resp[attr_ind] <- 0
    # number of subjects
    ns <- length(doses)

    permZ <- matrix(doses[permutations(ns,ns)],nrow=factorial(ns))

    # For each potential response vector under zero dose, compute obsT, exp, var

    zTr2 <- (trans(permZ) > thresh) %*% as.numeric(resp)
    zTr1 <- (trans(permZ) > thresh) %*% as.numeric(attr_resp)

    probs_upper2 <- exp(gamma*(permZ%*%as.numeric(resp)))
    probs_upper2 <- probs_upper2/sum(probs_upper2)
    probs_upper1 <- exp(gamma*(permZ%*%as.numeric(attr_resp)))
    probs_upper1 <- probs_upper1/sum(probs_upper1)


    lambda_1[j] = sum(probs_upper1*zTr1)
    lambda_2[j] = sum(probs_upper2*zTr2)
    omega_1[j] = lambda_1[j]*(1-lambda_1[j])
    omega_2[j] = lambda_2[j]*(1-lambda_2[j])

  }
  # data frame to pick the smallest mean and var differences
  vari = NULL
  pick = data.frame(mean = lambda_2-lambda_1, vari = omega_2 - omega_1, index = match_index)
  # contributing sets
  rel = pick[which( pick$index %in% contribute_ind),]
  # non-contributing sets
  rest = pick[-which( pick$index %in% contribute_ind),]
  # arrange the data by smallest mean change, if ties smallest var change
  rel = rel %>% arrange(mean, vari)
  # the indices to attribute to be least rejectable
  least_reject_ind = rel$index[1:a]
  # Make pi according to the a chosen sets in which to attribute
  pi =  lambda_1*(match_index %in% least_reject_ind) +
    lambda_2*(1-match_index %in% least_reject_ind)
  # If the worst-case mean is larger than observed value, return feasible.
  if (sum(pi) >= TT - a) {
    return("Feasible")
  }
  k = TT-a
  p_val = 1 - stats::pnorm((k-sum(pi))/sqrt(sum(pi*(1-pi))))
  return(list(p=p_val,cand = rel))
}

#' Inference for general attributable effects in sensitivity analysis with
#' continuous exposures and binary outcomes. Gurobi must be installed to use
#' this function.
#'
#' @param Z  A length N vector of (nonnegative) observed doses.
#' @param Q A length N vector of observed binary outcomes.
#' @param index A length N vector of indices indicating matched set membership.
#' @param gamma The nonnegative sensitivity parameter; gamma = 0 means no
#' unmeasured confounding.
#' @param alpha Level of the test.
#' @param monotone Whether to impose a monotonicity constraint on
#' the potential outcomes
#' @param Delta A numeric for the attributable effect to be tested for.
#' @param sign The sign of monotonicity - "positive" or "negative".
#' @param trans The transformation of the doses to use for the test statistic.
#' Default is the identity function.
#' @param direction A string indicating the direction of testing the attributable
#' effect with respect to Delta; "greater", "equal", or "less"
#' @param M The numeric penalty parameter to ensure correct sign relationship.
#' @param eps precision parameter for the objective function if solving feasible
#' = "Yes"
#' @param mv_threshold The number of allowed monotonicity violations.
#' @param baseline The baseline dose level for defining the threshold attributable
#' effect.
#' @param relax Whether to solve the continuous relaxation.
#' @param feasible Whether to look for a feasible solution or to find the optimal
#' solution.
#' @param alternative For constraining the rejection region using the test statistic
#' with baseline potential outcomes, should it be constructed with "lower" alternative,
#' "upper" alternative, or "both."
#'
#' @return A list containing the following:
#'
#' \item{sol}{A gurobi object containing the solution to the optimization. If
#' feasible is TRUE, then failing to find a solution indiciates rejection. If
#' feasible is FALSE, then optimal value greater than zero indiciates rejection.}
#' \item{attribut}{The attributable effect attained by the least rejectable
#'  allocation of potential outcomes and unmeasured confounders.}
#' \item{search}{A list of length the number of matched sets containing the
#' matrix of compatible baseline potential outcomes in each matched set.}
#' \item{null_exp}{The null expectation of the pivot attained by the least
#' rejectable allocation of potential outcomes and unmeasured confounders.}
#' \item{obsT}{The value of the pivot attained by the least rejectable
#'  allocation of potential outcomes and unmeasured confounders.}
#' \item{nps}{A vector with length the number of matched sets containing the
#' size of each matched set.}
#' \item{start.ind}{A vector with length the number of matched sets containing
#' the starting index of the decision variables pertaining to each matched set
#' in the optimization program.}
#' \item{exp_upper}{The expectation of the pivot attained by the unmeasured
#' confounders equalling the baseline potential outcome.}
#' \item{exp_lower}{The expectation of the pivot attained by the unmeasured
#' confounders equalling 1 minus the baseline potential outcome.}
#' \item{var_upper}{The variance of the pivot attained by the unmeasured
#' confounders equalling the baseline potential outcome.}
#' \item{var_lower}{The variance of the pivot attained by the unmeasured
#' confounders equalling 1 minus the baseline potential outcome.}
#' \item{dose}{A list of length the number of matched sets containing the
#' vector of doses observed in each matched set.}
#' \item{model}{The initialized gurobi model.}
#'
#' @export
#'
#' @examples
#' # To run the following example, Gurobi must be installed.
#' # Load the data
#' data <- treat_out_match
#' # Make a threshold at log(3.5) transformation function.
#' above = function(Z) { return(Z > log(3.5)) }
#' # Solve the mixed-integer program.
#' solution = dose_attributable_general(data$treat,
#' data$complain, data$match_ind, gamma=log(1),
#' alpha = 0.1, monotone = TRUE, trans = above,
#' Delta = 5, direction = "less", M = 10000, eps = 0.0001,
#' alternative = "upper", relax = FALSE, feasible = FALSE)
#'
dose_attributable_general <- function(Z, Q, index, gamma, alpha = 0.05, monotone = TRUE,
                                      Delta, sign = "positive", trans = identity,
                                      direction = "equal", M = 10000, eps = 0.01,
                                      alternative = "both",
                                      mv_threshold = NA, baseline = 0,
                                      relax = FALSE, feasible = TRUE) {
  # error-checking
  if (!requireNamespace("gurobi", quietly = TRUE)) {
    return("Please install gurobi")
  }
  if(!all(Q %in% c(0,1))) {
    stop("Non-binary outcomes")
  }
  if (gamma < 0) {
    stop("Negative gamma")
  }
  if (!(sign %in% c("positive","negative"))) {
    stop("invalid sign")
  }
  if (!(direction %in% c("greater","equal","less"))) {
    stop("invalid direction")
  }
  if (!(alternative %in% c("lower","upper","both"))) {
    stop("invalid alternative")
  }
  # The value of the test statistic under observed dose and response
  TT = sum(trans(Z)*Q)

  # Matched set indices
  match_index = unique(index)

  # number of matched sets
  nostratum <- length(unique(index))

  # observed test stat contribution
  obsT = c()
  # expectation contribution (one for u=rc, one for u=1-rc)
  exp_upper = c()
  exp_lower = c()
  # variance contribution
  var_upper = c()
  var_lower = c()
  # null expectation contribution
  null_exp = c()

  # number of monotonicity violations
  mv = c()
  if (is.na(mv_threshold)) {
    mv_threshold = length(Z)
  }

  # Variable to track number of decision variables in each matched set.
  nps = rep(0,nostratum)
  # Collect the potential response vectors to return
  search_space = vector(mode='list', length = nostratum)
  dose_space = vector(mode = 'list', length = nostratum)
  # Initialize for each matched set.
  for (j in 1:nostratum) {

    doses <- Z[which(index == match_index[j])]
    dose_space[[j]] = doses
    resp <- Q[which(index == match_index[j])]
    ns <- length(doses)

    permZ <- matrix(doses[permutations(ns,ns)],nrow=factorial(ns))

    # vector of responses that equal observed repsonse if dose 0, 0 otherwise
    check_comp = (doses <= baseline)*resp

    # responses under zero dose, all binary
    incomp = c()
    RZ = expand.grid(replicate(ns, 0:1, simplify = FALSE))
    # Get rid of incompatible ones based on observed data
    for (l in 1:nrow(RZ)) {
      pot_resp = RZ[l,]
      if (!all((pot_resp*(doses<=baseline))==check_comp)) {
        incomp = c(incomp,l)
      } else {
        if (monotone & sign == "positive" & !all(pot_resp<=resp)) {
          incomp = c(incomp,l)
        } else if(monotone & sign == "negative" & !all(pot_resp>=resp)) {
          incomp = c(incomp,l)
        }
      }
    }
    # Get rid of incompatible
    if (length(incomp) > 0) {
      RZ = RZ[-incomp,]
    }
    search_space[[j]] = RZ
    nps[j] = nrow(RZ)


    # For each potential response vector under zero dose, compute obsT, exp, var
    for (l in 1:nrow(RZ)) {
      pot_resp = RZ[l,]
      zTr <- trans(permZ) %*% as.numeric(pot_resp)
      if (gamma == 0) {
        probs_upper <- rep(1/nrow(permZ),nrow(permZ))
        probs_lower <- rep(1/nrow(permZ),nrow(permZ))
      } else {
        probs_upper <- exp(gamma*(permZ%*%as.numeric(pot_resp)))
        probs_lower <- exp(gamma*(permZ%*% (1-as.numeric(pot_resp))))
        probs_upper <- probs_upper/sum(probs_upper)
        probs_lower <- probs_lower/sum(probs_lower)
      }

      if (monotone & sign == "positive") {
        mv = c(mv, sum(pot_resp > resp))
      } else if (monotone & sign == "negative") {
        mv = c(mv, sum(pot_resp < resp))
      }


      obsT = c(obsT, (trans(doses)) %*% as.numeric(pot_resp))
      exp_upper = c(exp_upper, sum(probs_upper*zTr))
      exp_lower = c(exp_lower, sum(probs_lower*zTr))
      var_upper = c(var_upper, sum(probs_upper*(zTr^2))-(sum(probs_upper*zTr))^2)
      var_lower = c(var_lower, sum(probs_lower*(zTr^2))-(sum(probs_lower*zTr))^2)
      null_exp = c(null_exp, mean(zTr))
    }
  }

  # Make constraints
  # Each matched set only one contribution constraint
  nvars = sum(nps)+1+2
  Lin.eq <- matrix(0,nostratum,nvars)
  lin.eq.rhs <- rep(1,nostratum)
  start.ind = cumsum(c(1,nps))[-(nostratum+1)]

  for (j in 1:nostratum) {
    Lin.eq[j,start.ind[j]:(start.ind[j]+nps[j]-1)]<-1
  }
  # The attributable dose is greater than Delta or less than Delta
  Lin.ineq = matrix(0,nrow=6, nvars)

  lin.ineq.rhs <- rep(0,6)
  if (direction == "greater") {
    Lin.ineq[1,] = c(obsT,0,0,0)
    lin.ineq.rhs[1] = TT - Delta
  } else if (direction == "less") {
    Lin.ineq[1,] = -c(obsT,0,0,0)
    lin.ineq.rhs[1] = -TT + Delta
  } else if(direction == "equal") {
    Lin.ineq[1,] = c(obsT,0,0,0)
    lin.ineq.rhs[1] = TT - Delta
    Lin.ineq[6,] = -c(obsT,0,0,0)
    lin.ineq.rhs[6] = -TT + Delta
  }

  # Penalties for lower
  if (alternative %in% c("lower","both")) {
    Lin.ineq[2,] = c(obsT-exp_lower, 0,-M, 0)
    lin.ineq.rhs[2] = 0
    Lin.ineq[3,] = c(exp_lower-obsT, 0, M, 0)
    lin.ineq.rhs[3] = M
  }

  # Penalties for upper
  if (alternative %in% c("upper","both")) {
    Lin.ineq[4,] = c(obsT-exp_upper, 0, 0, M)
    lin.ineq.rhs[4] = M
    Lin.ineq[5,] = c(exp_upper-obsT, 0, 0, -M)
    lin.ineq.rhs[5] = 0
  }

  sense = c(rep("=", nostratum), rep("<=", 6))
  # Add monotonocity violation constraint
  mv.Lin.ineq <- c(mv,0,0,0)
  Lin.ineq <- rbind(Lin.ineq, mv.Lin.ineq)
  lin.ineq.rhs = c(lin.ineq.rhs, mv_threshold)
  sense = c(sense, "<=")

  # Upper and lower u quadratic constraints
  kap = stats::qchisq(1-2*alpha,1)

  # First quadratic constraint: for u = 1-r case
  qc1 <- list()
  qc1$Qc <- outer(c(obsT-exp_lower, 0, 0, 0),c(obsT-exp_lower, 0, 0, 0))
  qc1$q = c(-kap*var_lower,-1, -M, 0)
  qc1$rhs <- 0.0


  # Second quadratic constraint: for u = r case
  qc2 <- list()
  qc2$Qc <- outer(c(obsT-exp_upper, 0, 0, 0),c(obsT-exp_upper, 0, 0, 0))
  qc2$q = c(-kap*var_upper,-1, 0, -M)
  qc2$rhs <- 0.0

  print("initialized")
  # Solve the program, both max and min
  # Min
  model = list()
  model$A = rbind(Lin.eq, Lin.ineq)
  if (alternative == "upper") {
    model$quadcon <- list(qc2)
  } else if (alternative == "lower") {
    model$quadcon <- list(qc1)
  } else if (alternative == "both") {
    model$quadcon <- list(qc1, qc2)
  }

  model$lb = rep(-Inf,nvars)
  if (feasible) {
    model$ub = c(rep(Inf,nvars-3),-eps,Inf,Inf)
    model$obj = rep(0,nvars)
    model$vtypes = c(rep("B", nvars-3),"C","B","B")
  } else {
    model$ub = rep(Inf,nvars)
    model$obj = c(rep(0,nvars-3),1,0,0)
    if (relax) {
      model$vtypes = rep("C", nvars)
    } else{
      model$vtypes = c(rep("B", nvars-3),"C","B","B")
    }
  }

  model$sense = sense
  model$rhs = c(lin.eq.rhs,lin.ineq.rhs)

  model$modelsense = "min"
  #iis = gurobi_iis(model)
  solm = gurobi::gurobi(model,params = list(MIPgap = 0.01,
                                    WorkLimit = 1000))


  return(list(sol=solm, attribut = TT-sum(c(obsT,0,0,0)*solm$x),
              search = search_space, null_exp = sum(c(null_exp,0,0,0)*solm$x),
              obsT = sum(c(obsT,0,0,0)*solm$x),nps = nps, start.ind = start.ind,
              exp_upper = sum(c(exp_upper,0,0,0)*solm$x),
              exp_lower = sum(c(exp_lower,0,0,0)*solm$x),
              var_upper = sum(c(var_upper,0,0,0)*solm$x),
              var_lower = sum(c(var_lower,0,0,0)*solm$x),
              dose = dose_space, model = model))

}

#' Inference for threshold attributable effects in sensitivity analysis with
#' continuous exposures and binary outcomes. Gurobi must be installed to use
#' this function.
#'
#' @param Z  A length N vector of (nonnegative) observed doses.
#' @param Q A length N vector of observed binary outcomes.
#' @param index A length N vector of indices indicating matched set membership.
#' @param gamma The nonnegative sensitivity parameter; gamma = 0 means no
#' unmeasured confounding.
#' @param alpha Level of the test.
#' @param monotone Whether to impose a monotonicity constraint on
#' the potential outcomes
#' @param Delta A numeric for the attributable effect to be tested for.
#' @param sign The sign of monotonicity - "positive" or "negative".
#' @param direction A string indicating the direction of testing the attributable
#' effect with respect to Delta; "greater", "equal", or "less"
#' @param threshold The threshold for the TAE.
#' @param M The numeric penalty parameter to ensure correct sign relationship.
#' @param eps precision parameter for the objective function if solving feasible
#' = "Yes"
#' @param mv_threshold The number of allowed monotonicity violations.
#' @param baseline The baseline dose level for defining the threshold attributable
#' effect.
#' @param relax Whether to solve the continuous relaxation.
#' @param feasible Whether to look for a feasible solution or to find the optima.
#'
#' @return A list containing the following:
#'
#' \item{sol}{A gurobi object containing the solution to the optimization. If
#' feasible is TRUE, then failing to find a solution indiciates rejection. If
#' feasible is FALSE, then optimal value greater than zero indiciates rejection.}
#' \item{attribut}{The attributable effect attained by the least rejectable
#'  allocation of potential outcomes and unmeasured confounders.}
#' \item{search}{A list of length the number of matched sets containing the
#' matrix of compatible baseline potential outcomes in each matched set.}
#' \item{null_exp}{The null expectation of the pivot attained by the least
#' rejectable allocation of potential outcomes and unmeasured confounders.}
#' \item{obsT}{The value of the pivot attained by the least rejectable
#'  allocation of potential outcomes and unmeasured confounders.}
#' \item{nps}{A vector with length the number of matched sets containing the
#' size of each matched set.}
#' \item{start.ind}{A vector with length the number of matched sets containing
#' the starting index of the decision
#' variables pertaining to each matched set in the optimization program.}
#' \item{exp_upper}{The expectation of the pivot attained by the unmeasured
#' confounders equalling the baseline potential outcome.}
#' \item{var_upper}{The variance of the pivot attained by the unmeasured
#' confounders equalling the baseline potential outcome.}
#' \item{dose}{A list of length the number of matched sets containing the
#' vector of doses observed in each matched set.}
#' \item{model}{The initialized gurobi model.}
#' @export
#'
#' @examples
#' # To run the following example, Gurobi must be installed.
#' # Load the data
#' data <- treat_out_match
#' # Solve the mixed-integer program.
#' solution = dose_thresh_attributable_one_sided(data$treat,
#' data$complain, data$match_ind,
#' gamma=log(1), alpha = 0.1, monotone = TRUE, Delta = 5,
#'  direction = "less", threshold = log(3.5),M = 10000,
#'  eps = 0.0001,relax = FALSE, feasible = FALSE)
dose_thresh_attributable_one_sided <- function(Z, Q, index, gamma, alpha = 0.05, monotone = TRUE,
                                               Delta, sign = "positive", direction = "equal",
                                               threshold = 0,
                                               M = 10000, eps = 0.01, mv_threshold = NA,
                                               baseline = 0,
                                               relax = FALSE, feasible = TRUE) {
  # error-checking
  if (!requireNamespace("gurobi", quietly = TRUE)) {
    return("Please install gurobi")
  }
  if(!all(Q %in% c(0,1))) {
    stop("Non-binary outcomes")
  }
  if (gamma < 0) {
    stop("Negative gamma")
  }
  if (!(sign %in% c("positive","negative"))) {
    stop("invalid sign")
  }
  if (!(direction %in% c("greater","equal","less"))) {
    stop("invalid direction")
  }
  # The value of the test statistic under observed dose and response
  Zthresh = Z > threshold
  TT = sum(Zthresh*Q)

  # Matched set indices
  match_index = unique(index)

  # number of matched sets
  nostratum <- length(unique(index))

  # observed test stat contribution
  obsT = c()

  # expectation contribution (one for u=rc, one for u=1-rc)
  exp_upper = c()

  # variance contribution
  var_upper = c()

  # null expectation contribution
  null_exp = c()

  # number of monotonicity violations
  mv = c()
  if (is.na(mv_threshold)) {
    mv_threshold = length(Z)
  }

  # Variable to track number of decision variables in each matched set.
  nps = rep(0,nostratum)
  # Collect the potential response vectors to return
  search_space = vector(mode='list', length = nostratum)
  dose_space = vector(mode = 'list', length = nostratum)
  # Initialize for each matched set.
  for (j in 1:nostratum) {

    doses <- Z[which(index == match_index[j])]
    dose_space[[j]] = doses
    resp <- Q[which(index == match_index[j])]
    ns <- length(doses)

    permZ <- matrix(doses[permutations(ns,ns)],nrow=factorial(ns))

    # vector of responses that equal observed response if dose 0, 0 otherwise
    check_comp = (doses <= baseline)*resp

    # responses under zero dose, all binary
    incomp = c()
    RZ = expand.grid(replicate(ns, 0:1, simplify = FALSE))
    # Get rid of incompatible ones based on observed data
    for (l in 1:nrow(RZ)) {
      pot_resp = RZ[l,]
      if (!all((pot_resp*(doses<=baseline))==check_comp)) {
        incomp = c(incomp,l)
      } else {
        if (monotone & sign == "positive" & !all(pot_resp<=resp)) {
          incomp = c(incomp,l)
        } else if(monotone & sign == "negative" & !all(pot_resp>=resp)) {
          incomp = c(incomp,l)
        }
      }
    }
    # Get rid of incompatible
    if (length(incomp) > 0) {
      RZ = RZ[-incomp,]
    }
    search_space[[j]] = RZ
    nps[j] = nrow(RZ)


    # For each potential response vector under baseline dose, compute obsT, exp, var
    for (l in 1:nrow(RZ)) {
      pot_resp = RZ[l,]
      zTr <- (permZ > threshold) %*% as.numeric(pot_resp)
      if (gamma == 0) {
        probs_upper <- rep(1/nrow(permZ),nrow(permZ))
      } else {
        probs_upper <- exp(gamma*(permZ%*%as.numeric(pot_resp)))
        probs_upper <- probs_upper/sum(probs_upper)

      }

      if (monotone & sign == "positive") {
        mv = c(mv, sum(pot_resp > resp))
      } else if (monotone & sign == "negative") {
        mv = c(mv, sum(pot_resp < resp))
      }

      # Make the relevant vectors
      obsT = c(obsT, (doses>threshold) %*% as.numeric(pot_resp))
      exp_upper = c(exp_upper, sum(probs_upper*zTr))
      var_upper = c(var_upper, sum(probs_upper*(zTr^2))-(sum(probs_upper*zTr))^2)
      null_exp = c(null_exp, mean(zTr))
    }
  }

  # Make constraints
  # Each matched set only one contribution constraint
  nvars = sum(nps)+1+2
  Lin.eq <- matrix(0,nostratum,nvars)
  lin.eq.rhs <- rep(1,nostratum)
  start.ind = cumsum(c(1,nps))[-(nostratum+1)]

  for (j in 1:nostratum) {
    Lin.eq[j,start.ind[j]:(start.ind[j]+nps[j]-1)]<-1
  }
  # The attributable dose is greater than Delta or less than Delta, or equal to Delta

  if (direction != "equal") {
    Lin.ineq = matrix(0,nrow=3, nvars)
    lin.ineq.rhs <- rep(0,3)
    if (direction == "greater") {
      Lin.ineq[1,] = c(obsT,0,0,0)
      lin.ineq.rhs[1] = TT - Delta
    } else if (direction == "less") {
      Lin.ineq[1,] = -c(obsT,0,0,0)
      lin.ineq.rhs[1] = -TT + Delta
    }
    # Penalties for upper
    Lin.ineq[2,] = c(obsT-exp_upper, 0, 0, M)
    lin.ineq.rhs[2] = M
    Lin.ineq[3,] = c(exp_upper-obsT, 0, 0, -M)
    lin.ineq.rhs[3] = 0
    sense = c(rep("=",nostratum),rep("<=",3))
  } else {
    Lin.ineq = matrix(0,nrow=4, nvars)
    lin.ineq.rhs <- rep(0,4)

    Lin.ineq[1,] = c(obsT,0,0,0)
    lin.ineq.rhs[1] = TT - Delta

    Lin.ineq[2,] = -c(obsT,0,0,0)
    lin.ineq.rhs[2] = -TT + Delta

    # Penalties for upper
    Lin.ineq[3,] = c(obsT-exp_upper, 0, 0, M)
    lin.ineq.rhs[3] = M
    Lin.ineq[4,] = c(exp_upper-obsT, 0, 0, -M)
    lin.ineq.rhs[4] = 0
    sense = c(rep("=",nostratum),rep("<=",4))
  }
  # Add monotonocity violation constraint
  mv.Lin.ineq <- c(mv,0,0,0)
  Lin.ineq <- rbind(Lin.ineq, mv.Lin.ineq)
  lin.ineq.rhs = c(lin.ineq.rhs, mv_threshold)
  sense = c(sense, "<=")



  # Upper and lower u quadratic constraints
  kap = stats::qchisq(1-2*alpha,1)


  # Second quadratic constraint: for u = r case
  qc2 <- list()
  qc2$Qc <- outer(c(obsT-exp_upper, 0, 0, 0),c(obsT-exp_upper, 0, 0, 0))
  qc2$q = c(-kap*var_upper,-1, 0, -M)
  qc2$rhs <- 0.0

  print("initialized")
  # Solve the program, both max and min
  # Min
  model = list()
  model$A = rbind(Lin.eq, Lin.ineq)
  model$quadcon <- list(qc2)
  model$lb = rep(-Inf,nvars)
  if (feasible) {
    model$ub = c(rep(Inf,nvars-3),-eps,Inf,Inf)
    model$obj = rep(0,nvars)
    model$vtypes = c(rep("B", nvars-3),"C","B","B")
  } else {
    model$ub = rep(Inf,nvars)
    model$obj = c(rep(0,nvars-3),1,0,0)
    if (relax) {
      model$vtypes = rep("C", nvars)
    } else{
      model$vtypes = c(rep("B", nvars-3),"C","B","B")
    }
  }

  model$sense = sense
  model$rhs = c(lin.eq.rhs,lin.ineq.rhs)

  model$modelsense = "min"
  #iis = gurobi_iis(model)
  solm = gurobi::gurobi(model,params = list(MIPgap = 0.01,
                                    WorkLimit = 1000))


  return(list(sol=solm, attribut = TT-sum(c(obsT,0,0,0)*solm$x),
              search = search_space, null_exp = sum(c(null_exp,0,0,0)*solm$x),
              obsT = sum(c(obsT,0,0,0)*solm$x),nps = nps, start.ind = start.ind,
              exp_upper = sum(c(exp_upper,0,0,0)*solm$x),
              var_upper = sum(c(var_upper,0,0,0)*solm$x),
              dose = dose_space, model = model))

}
