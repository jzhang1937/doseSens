# Helper functions
# from a given vector, return all unique permutations of it
# Source:
# https://stackoverflow.com/questions/5671149/permute-all-unique-enumerations-of-a-vector-in-r
uniqueperm <- function(d) {
  dat <- factor(d)
  N <- length(dat)
  n <- tabulate(dat)
  ng <- length(n)
  if(ng==1) return(d)
  a <- N-c(0,cumsum(n))[-(ng+1)]
  foo <- lapply(1:ng, function(i) matrix(utils::combn(a[i],n[i]),nrow=n[i]))
  out <- matrix(NA, nrow=N, ncol=prod(sapply(foo, ncol)))
  xxx <- c(0,cumsum(sapply(foo, nrow)))
  xxx <- cbind(xxx[-length(xxx)]+1, xxx[-1])
  miss <- matrix(1:N,ncol=1)
  for(i in seq_len(length(foo)-1)) {
    l1 <- foo[[i]]
    nn <- ncol(miss)
    miss <- matrix(rep(miss, ncol(l1)), nrow=nrow(miss))
    k <- (rep(0:(ncol(miss)-1), each=nrow(l1)))*nrow(miss) +
      l1[,rep(1:ncol(l1), each=nn)]
    out[xxx[i,1]:xxx[i,2],] <- matrix(miss[k], ncol=ncol(miss))
    miss <- matrix(miss[-k], ncol=ncol(miss))
  }
  k <- length(foo)
  out[xxx[k,1]:xxx[k,2],] <- miss
  out <- out[rank(as.numeric(dat), ties.method="first"),]
  foo <- cbind(as.vector(out), as.vector(col(out)))
  out[foo] <- d
  t(out)
}


#' A helper that takes a gurobi model object outputted from
#' dose_attributable_general or dose_thresh_attributable_one_sided
#' and changes the Delta parameter. Saves computation time by directly editing
#' the constraint involving Delta without having to reinitialize the other
#' constraints that are kept the same. Outputs a list analogous to output from
#' dose_attributable_general or dose_thresh_attributable_one_sided.
#'
#' @param model A gurobi model object outputted from dose_attributable_general.
#' @param Delta The new Delta to test for.
#' @param direction The new direction to test
#' @param TT The observed test statistic.
#'
#' @return A gurobi model and solution.
#'
change_Delta <- function(model, Delta, direction = "equal", TT) {
  requireNamespace("gurobi", quietly = TRUE)
  new.rhs <- model$rhs
  last <- length(model$rhs)
  if (direction == "equal") {
    new.rhs[last-4] <- TT - Delta
    new.rhs[last-3] <- -TT + Delta
  } else if (direction == "less") {
    new.rhs[last-3] <- -TT + Delta
  } else if (direction == "greater") {
    new.rhs[last-3] <- TT - Delta
  }
  model$rhs <- new.rhs
  solm = gurobi::gurobi(model,params = list(MIPgap = 0.01,
                                    WorkLimit = 1000))


  return(list(sol=solm, model = model))
}
