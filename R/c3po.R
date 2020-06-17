r2d2 <- function(refdata,
                bc1d,
                icond = c(1),
                lag_search = 0,
                lag_keep = 0) {
  # Looks for blocks of size lag_search+1
  # Keeps only the lag_keep+1 last values of the block
  if (lag_search < lag_keep) {
    stop("lag_search has to be greater or equal to lag_keep")
  }

  P <- ncol(bc1d) # number of var (including icond !!)
  Ntimes_BC <- nrow(bc1d)
  Ntimes_REF <-  nrow(refdata)
  stopifnot(lag_search <= min(Ntimes_BC, Ntimes_REF))
  if (lag_search > min(Ntimes_BC, Ntimes_REF)) {
    stop("lag_search has to be less or equal to the number of observations in refdata or bc1d")
  }
  r2d2_bc = .Call(C_r2d2, refdata, bc1d, icond - 1, lag_search, lag_keep)
  return(
    list(
      r2d2_bc = r2d2_bc[[4]],
      visited_time = r2d2_bc[[7]],
      time_bestanalogue = r2d2_bc[[5]] + 1,
      dist_bestanalogue = sqrt(r2d2_bc[[6]])
    )
  )
}

