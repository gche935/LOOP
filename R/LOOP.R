
## ----- Sub-function gmc & gm - Calculate the group-mean centered and group mean ----- ##
gmc <- function(x) { x - mean(x, na.rm=TRUE) }
gm  <- function(x) { mean(x, na.rm=TRUE) }
## -------------------------------------------------------- ##


## ----- Sub-function mL2 - Create the group-mean centered and group mean variables ----- ##
mL2 <- function(data.source, id, mL2.variables) {
  data_centered <- data.source %>%
  group_by(id) %>%
  mutate(across(all_of(mL2.variables), ~gmc(.x), .names = "gmc_{.col}")) %>%
  mutate(across(all_of(mL2.variables), ~gm(.x), .names = "gm_{.col}")) %>%
  ungroup()
}

## -------------------------------------------------------- ##



## ----- Sub-Function Invariance for testing invariance of parameters ----- ##

Invariance <- function(parEst, pest2, pest3, no.path, no.waves, lag, Type1, Type1Adj, X, Y, Z, W) {

  mcmc <- MASS::mvrnorm(n=1000000, mu=pest2, Sigma=pest3, tol = 1e-6)  # Run 1,000,000 simulations
  names(pest2) <-colnames(pest3)  # Save Parameter Names to Estimated Parameters
  b.no <- nrow(mcmc)  # No. of successful simulated samples


  ## -- Differences in path coefficients pXX -- ##
  if (any(parEst[,4] == "pXX21")) {

    # -- Lag = 1 -- #
    for (j in 3:no.waves) {
      for (i in j:no.waves) {
        mcmcA <- (mcmc[, paste("pXX", i, i-1, sep="")] - mcmc[, paste("pXX", i-j+2, i-j+1, sep="")])
        mcmc <- cbind(mcmc, mcmcA)
        colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pXX", i, i-1, "-pXX", i-j+2, i-j+1, sep="")
        pest2A <- (pest2[paste("pXX", i, i-1, sep="")] - pest2[paste("pXX", i-j+2, i-j+1, sep="")])
        names(pest2A) <- paste("pXX", i, i-1, "-pXX", i-j+2, i-j+1, sep="")
        pest2 <- append(pest2, pest2A)

        mcmcA <- (mcmc[, paste("pYY", i, i-1, sep="")] - mcmc[, paste("pYY", i-j+2, i-j+1, sep="")])
        mcmc <- cbind(mcmc, mcmcA)
        colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pYY", i, i-1, "-pYY", i-j+2, i-j+1, sep="")
        pest2A <- (pest2[paste("pYY", i, i-1, sep="")] - pest2[paste("pYY", i-j+2, i-j+1, sep="")])
        names(pest2A) <- paste("pYY", i, i-1, "-pYY", i-j+2, i-j+1, sep="")
        pest2 <- append(pest2, pest2A)

        if (Z != "NULL") {
          mcmcA <- (mcmc[, paste("pZZ", i, i-1, sep="")] - mcmc[, paste("pZZ", i-j+2, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pZZ", i, i-1, "-pZZ", i-j+2, i-j+1, sep="")
          pest2A <- (pest2[paste("pZZ", i, i-1, sep="")] - pest2[paste("pZZ", i-j+2, i-j+1, sep="")])
          names(pest2A) <- paste("pZZ", i, i-1, "-pZZ", i-j+2, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)
        } # end (if Z)

        if (W != "NULL") {
          mcmcA <- (mcmc[, paste("pWW", i, i-1, sep="")] - mcmc[, paste("pWW", i-j+2, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pWW", i, i-1, "-pWW", i-j+2, i-j+1, sep="")
          pest2A <- (pest2[paste("pWW", i, i-1, sep="")] - pest2[paste("pWW", i-j+2, i-j+1, sep="")])
          names(pest2A) <- paste("pWW", i, i-1, "-pWW", i-j+2, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)
        } # end (if W != "NULL")
      } # end (for i)
    } # end (for j)

    # -- Lag = 2 -- #
    if (lag == 2 & no.waves > 3) {
      for (j in 4:no.waves) {
        for (i in j:no.waves) {
          mcmcA <- (mcmc[, paste("pXX", i, i-2, sep="")] - mcmc[, paste("pXX", i-j+3, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pXX", i, i-2, "-pXX", i-j+3, i-j+1, sep="")
          pest2A <- (pest2[paste("pXX", i, i-2, sep="")] - pest2[paste("pXX", i-j+3, i-j+1, sep="")])
          names(pest2A) <- paste("pXX", i, i-2, "-pXX", i-j+3, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("pYY", i, i-2, sep="")] - mcmc[, paste("pYY", i-j+3, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pYY", i, i-2, "-pYY", i-j+3, i-j+1, sep="")
          pest2A <- (pest2[paste("pYY", i, i-2, sep="")] - pest2[paste("pYY", i-j+3, i-j+1, sep="")])
          names(pest2A) <- paste("pYY", i, i-2, "-pYY", i-j+3, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          if (Z != "NULL") {
            mcmcA <- (mcmc[, paste("pZZ", i, i-2, sep="")] - mcmc[, paste("pZZ", i-j+3, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pZZ", i, i-2, "-pZZ", i-j+3, i-j+1, sep="")
            pest2A <- (pest2[paste("pZZ", i, i-2, sep="")] - pest2[paste("pZZ", i-j+3, i-j+1, sep="")])
            names(pest2A) <- paste("pZZ", i, i-2, "-pZZ", i-j+3, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)
          } # end (if Z != "NULL")

          if (W != "NULL") {
            mcmcA <- (mcmc[, paste("pWW", i, i-2, sep="")] - mcmc[, paste("pWW", i-j+3, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pWW", i, i-2, "-pWW", i-j+3, i-j+1, sep="")
            pest2A <- (pest2[paste("pWW", i, i-2, sep="")] - pest2[paste("pWW", i-j+3, i-j+1, sep="")])
            names(pest2A) <- paste("pWW", i, i-2, "-pWW", i-j+3, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)
          } # end (if W != "NULL")
        } # end (for i)
      } # end (for j)
    } # end (lag == 2)

    # -- Lag = 3 -- #
    if (lag == 3 & no.waves > 4) {
      for (j in 5:no.waves) {
        for (i in j:no.waves) {
          mcmcA <- (mcmc[, paste("pXX", i, i-3, sep="")] - mcmc[, paste("pXX", i-j+4, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pXX", i, i-3, "-pXX", i-j+4, i-j+1, sep="")
          pest2A <- (pest2[paste("pXX", i, i-3, sep="")] - pest2[paste("pXX", i-j+4, i-j+1, sep="")])
          names(pest2A) <- paste("pXX", i, i-3, "-pXX", i-j+4, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("pYY", i, i-3, sep="")] - mcmc[, paste("pYY", i-j+4, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pYY", i, i-3, "-pYY", i-j+4, i-j+1, sep="")
          pest2A <- (pest2[paste("pYY", i, i-3, sep="")] - pest2[paste("pYY", i-j+4, i-j+1, sep="")])
          names(pest2A) <- paste("pYY", i, i-3, "-pYY", i-j+4, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          if (Z != "NULL") {
            mcmcA <- (mcmc[, paste("pZZ", i, i-3, sep="")] - mcmc[, paste("pZZ", i-j+4, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pZZ", i, i-3, "-pZZ", i-j+4, i-j+1, sep="")
            pest2A <- (pest2[paste("pZZ", i, i-3, sep="")] - pest2[paste("pZZ", i-j+4, i-j+1, sep="")])
            names(pest2A) <- paste("pZZ", i, i-3, "-pZZ", i-j+4, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)
          } # end (if Z != "NULL")

          if (W != "NULL") {
            mcmcA <- (mcmc[, paste("pWW", i, i-3, sep="")] - mcmc[, paste("pWW", i-j+4, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pWW", i, i-3, "-pWW", i-j+4, i-j+1, sep="")
            pest2A <- (pest2[paste("pWW", i, i-3, sep="")] - pest2[paste("pWW", i-j+4, i-j+1, sep="")])
            names(pest2A) <- paste("pWW", i, i-3, "-pWW", i-j+4, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)
          } # end (if W != "NULL")
        } # end (for i)
      } # end (for j)
    } # end (lag == 3)

    # -- Lag = 4 -- #
    if (lag == 4 & no.waves > 5) {
      for (j in 6:no.waves) {
        for (i in j:no.waves) {
          mcmcA <- (mcmc[, paste("pXX", i, i-4, sep="")] - mcmc[, paste("pXX", i-j+5, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pXX", i, i-4, "-pXX", i-j+5, i-j+1, sep="")
          pest2A <- (pest2[paste("pXX", i, i-4, sep="")] - pest2[paste("pXX", i-j+5, i-j+1, sep="")])
          names(pest2A) <- paste("pXX", i, i-4, "-pXX", i-j+5, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("pYY", i, i-4, sep="")] - mcmc[, paste("pYY", i-j+5, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pYY", i, i-4, "-pYY", i-j+5, i-j+1, sep="")
          pest2A <- (pest2[paste("pYY", i, i-4, sep="")] - pest2[paste("pYY", i-j+5, i-j+1, sep="")])
          names(pest2A) <- paste("pYY", i, i-4, "-pYY", i-j+5, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          if (Z != "NULL") {
            mcmcA <- (mcmc[, paste("pZZ", i, i-4, sep="")] - mcmc[, paste("pZZ", i-j+5, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pZZ", i, i-4, "-pZZ", i-j+5, i-j+1, sep="")
            pest2A <- (pest2[paste("pZZ", i, i-4, sep="")] - pest2[paste("pZZ", i-j+5, i-j+1, sep="")])
            names(pest2A) <- paste("pZZ", i, i-4, "-pZZ", i-j+5, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)
          } # end (if Z != "NULL")

          if (W != "NULL") {
            mcmcA <- (mcmc[, paste("pWW", i, i-4, sep="")] - mcmc[, paste("pWW", i-j+5, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pWW", i, i-4, "-pWW", i-j+5, i-j+1, sep="")
            pest2A <- (pest2[paste("pWW", i, i-4, sep="")] - pest2[paste("pWW", i-j+5, i-j+1, sep="")])
            names(pest2A) <- paste("pWW", i, i-4, "-pWW", i-j+5, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)
          } # end (if W != "NULL")
        } # end (for i)
      } # end (for j)
    } # end (lag == 4)
  } # end (if pXX21)
  ## ----- end (Difference in path coefficients) ----- ##


  ## -- Differences in path coefficients -- ##
  if (any(parEst[,4] == "pXY21")) {

    # -- Lag = 1 -- #
    for (j in 3:no.waves) {
      for (i in j:no.waves) {
        mcmcA <- (mcmc[, paste("pXY", i, i-1, sep="")] - mcmc[, paste("pXY", i-j+2, i-j+1, sep="")])
        mcmc <- cbind(mcmc, mcmcA)
        colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pXY", i, i-1, "-pXY", i-j+2, i-j+1, sep="")
        pest2A <- (pest2[paste("pXY", i, i-1, sep="")] - pest2[paste("pXY", i-j+2, i-j+1, sep="")])
        names(pest2A) <- paste("pXY", i, i-1, "-pXY", i-j+2, i-j+1, sep="")
        pest2 <- append(pest2, pest2A)

        mcmcA <- (mcmc[, paste("pYX", i, i-1, sep="")] - mcmc[, paste("pYX", i-j+2, i-j+1, sep="")])
        mcmc <- cbind(mcmc, mcmcA)
        colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pYX", i, i-1, "-pYX", i-j+2, i-j+1, sep="")
        pest2A <- (pest2[paste("pYX", i, i-1, sep="")] - pest2[paste("pYX", i-j+2, i-j+1, sep="")])
        names(pest2A) <- paste("pYX", i, i-1, "-pYX", i-j+2, i-j+1, sep="")
        pest2 <- append(pest2, pest2A)

        if (Z != "NULL") {
          mcmcA <- (mcmc[, paste("pXZ", i, i-1, sep="")] - mcmc[, paste("pXZ", i-j+2, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pXZ", i, i-1, "-pXZ", i-j+2, i-j+1, sep="")
          pest2A <- (pest2[paste("pXZ", i, i-1, sep="")] - pest2[paste("pXZ", i-j+2, i-j+1, sep="")])
          names(pest2A) <- paste("pXZ", i, i-1, "-pXZ", i-j+2, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("pYZ", i, i-1, sep="")] - mcmc[, paste("pYZ", i-j+2, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pYZ", i, i-1, "-pYZ", i-j+2, i-j+1, sep="")
          pest2A <- (pest2[paste("pYZ", i, i-1, sep="")] - pest2[paste("pYZ", i-j+2, i-j+1, sep="")])
          names(pest2A) <- paste("pYZ", i, i-1, "-pYZ", i-j+2, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("pZX", i, i-1, sep="")] - mcmc[, paste("pZX", i-j+2, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pZX", i, i-1, "-pZX", i-j+2, i-j+1, sep="")
          pest2A <- (pest2[paste("pZX", i, i-1, sep="")] - pest2[paste("pZX", i-j+2, i-j+1, sep="")])
          names(pest2A) <- paste("pZX", i, i-1, "-pZX", i-j+2, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("pZY", i, i-1, sep="")] - mcmc[, paste("pZY", i-j+2, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pZY", i, i-1, "-pZY", i-j+2, i-j+1, sep="")
          pest2A <- (pest2[paste("pZY", i, i-1, sep="")] - pest2[paste("pZY", i-j+2, i-j+1, sep="")])
          names(pest2A) <- paste("pZY", i, i-1, "-pZY", i-j+2, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)
        } # end (if Z != "NULL")

        if (W != "NULL") {
          mcmcA <- (mcmc[, paste("pXW", i, i-1, sep="")] - mcmc[, paste("pXW", i-j+2, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pXW", i, i-1, "-pXW", i-j+2, i-j+1, sep="")
          pest2A <- (pest2[paste("pXW", i, i-1, sep="")] - pest2[paste("pXW", i-j+2, i-j+1, sep="")])
          names(pest2A) <- paste("pXW", i, i-1, "-pXW", i-j+2, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("pYW", i, i-1, sep="")] - mcmc[, paste("pYW", i-j+2, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pYW", i, i-1, "-pYW", i-j+2, i-j+1, sep="")
          pest2A <- (pest2[paste("pYW", i, i-1, sep="")] - pest2[paste("pYW", i-j+2, i-j+1, sep="")])
          names(pest2A) <- paste("pYW", i, i-1, "-pYW", i-j+2, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("pZW", i, i-1, sep="")] - mcmc[, paste("pZW", i-j+2, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pZW", i, i-1, "-pZW", i-j+2, i-j+1, sep="")
          pest2A <- (pest2[paste("pZW", i, i-1, sep="")] - pest2[paste("pZW", i-j+2, i-j+1, sep="")])
          names(pest2A) <- paste("pZW", i, i-1, "-pZW", i-j+2, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("pWX", i, i-1, sep="")] - mcmc[, paste("pWX", i-j+2, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pWX", i, i-1, "-pWX", i-j+2, i-j+1, sep="")
          pest2A <- (pest2[paste("pWX", i, i-1, sep="")] - pest2[paste("pWX", i-j+2, i-j+1, sep="")])
          names(pest2A) <- paste("pWX", i, i-1, "-pWX", i-j+2, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("pWY", i, i-1, sep="")] - mcmc[, paste("pWY", i-j+2, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pWY", i, i-1, "-pWY", i-j+2, i-j+1, sep="")
          pest2A <- (pest2[paste("pWY", i, i-1, sep="")] - pest2[paste("pWY", i-j+2, i-j+1, sep="")])
          names(pest2A) <- paste("pWY", i, i-1, "-pWY", i-j+2, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("pWZ", i, i-1, sep="")] - mcmc[, paste("pWZ", i-j+2, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pWZ", i, i-1, "-pWZ", i-j+2, i-j+1, sep="")
          pest2A <- (pest2[paste("pWZ", i, i-1, sep="")] - pest2[paste("pWZ", i-j+2, i-j+1, sep="")])
          names(pest2A) <- paste("pWZ", i, i-1, "-pWZ", i-j+2, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)
        } # end (if W != "NULL")
      } # end (for i)
    } # end (for j)

    # -- Lag = 2 -- #
    if (lag == 2 & no.waves > 3) {
      for (j in 4:no.waves) {
        for (i in j:no.waves) {
          mcmcA <- (mcmc[, paste("pXY", i, i-2, sep="")] - mcmc[, paste("pXY", i-j+3, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pXY", i, i-2, "-pXY", i-j+3, i-j+1, sep="")
          pest2A <- (pest2[paste("pXY", i, i-2, sep="")] - pest2[paste("pXY", i-j+3, i-j+1, sep="")])
          names(pest2A) <- paste("pXY", i, i-2, "-pXY", i-j+3, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("pYX", i, i-2, sep="")] - mcmc[, paste("pYX", i-j+3, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pYX", i, i-2, "-pYX", i-j+3, i-j+1, sep="")
          pest2A <- (pest2[paste("pYX", i, i-2, sep="")] - pest2[paste("pYX", i-j+3, i-j+1, sep="")])
          names(pest2A) <- paste("pYX", i, i-2, "-pYX", i-j+3, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          if (Z != "NULL") {
            mcmcA <- (mcmc[, paste("pXZ", i, i-2, sep="")] - mcmc[, paste("pXZ", i-j+3, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pXZ", i, i-2, "-pXZ", i-j+3, i-j+1, sep="")
            pest2A <- (pest2[paste("pXZ", i, i-2, sep="")] - pest2[paste("pXZ", i-j+3, i-j+1, sep="")])
            names(pest2A) <- paste("pXZ", i, i-2, "-pXZ", i-j+3, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)

            mcmcA <- (mcmc[, paste("pYZ", i, i-2, sep="")] - mcmc[, paste("pYZ", i-j+3, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pYZ", i, i-2, "-pYZ", i-j+3, i-j+1, sep="")
            pest2A <- (pest2[paste("pYZ", i, i-2, sep="")] - pest2[paste("pYZ", i-j+3, i-j+1, sep="")])
            names(pest2A) <- paste("pYZ", i, i-2, "-pYZ", i-j+3, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)

            mcmcA <- (mcmc[, paste("pZX", i, i-2, sep="")] - mcmc[, paste("pZX", i-j+3, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pZX", i, i-2, "-pZX", i-j+3, i-j+1, sep="")
            pest2A <- (pest2[paste("pZX", i, i-2, sep="")] - pest2[paste("pZX", i-j+3, i-j+1, sep="")])
            names(pest2A) <- paste("pZX", i, i-2, "-pZX", i-j+3, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)

            mcmcA <- (mcmc[, paste("pZY", i, i-2, sep="")] - mcmc[, paste("pZY", i-j+3, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pZY", i, i-2, "-pZY", i-j+3, i-j+1, sep="")
            pest2A <- (pest2[paste("pZY", i, i-2, sep="")] - pest2[paste("pZY", i-j+3, i-j+1, sep="")])
            names(pest2A) <- paste("pZY", i, i-2, "-pZY", i-j+3, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)
          } # end (if Z != "NULL")

          if (W != "NULL") {
            mcmcA <- (mcmc[, paste("pXW", i, i-2, sep="")] - mcmc[, paste("pXW", i-j+3, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pXW", i, i-2, "-pXW", i-j+3, i-j+1, sep="")
            pest2A <- (pest2[paste("pXW", i, i-2, sep="")] - pest2[paste("pXW", i-j+3, i-j+1, sep="")])
            names(pest2A) <- paste("pXW", i, i-2, "-pXW", i-j+3, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)

            mcmcA <- (mcmc[, paste("pYW", i, i-2, sep="")] - mcmc[, paste("pYW", i-j+3, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pYW", i, i-2, "-pYW", i-j+3, i-j+1, sep="")
            pest2A <- (pest2[paste("pYW", i, i-2, sep="")] - pest2[paste("pYW", i-j+3, i-j+1, sep="")])
            names(pest2A) <- paste("pYW", i, i-2, "-pYW", i-j+3, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)

            mcmcA <- (mcmc[, paste("pZW", i, i-2, sep="")] - mcmc[, paste("pZW", i-j+3, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pZW", i, i-2, "-pZW", i-j+3, i-j+1, sep="")
            pest2A <- (pest2[paste("pZW", i, i-2, sep="")] - pest2[paste("pZW", i-j+3, i-j+1, sep="")])
            names(pest2A) <- paste("pZW", i, i-2, "-pZW", i-j+3, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)

            mcmcA <- (mcmc[, paste("pWX", i, i-2, sep="")] - mcmc[, paste("pWX", i-j+3, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pWX", i, i-2, "-pWX", i-j+3, i-j+1, sep="")
            pest2A <- (pest2[paste("pWX", i, i-2, sep="")] - pest2[paste("pWX", i-j+3, i-j+1, sep="")])
            names(pest2A) <- paste("pWX", i, i-2, "-pWX", i-j+3, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)

            mcmcA <- (mcmc[, paste("pWY", i, i-2, sep="")] - mcmc[, paste("pWY", i-j+3, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pWY", i, i-2, "-pWY", i-j+3, i-j+1, sep="")
            pest2A <- (pest2[paste("pWY", i, i-2, sep="")] - pest2[paste("pWY", i-j+3, i-j+1, sep="")])
            names(pest2A) <- paste("pWY", i, i-2, "-pWY", i-j+3, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)

            mcmcA <- (mcmc[, paste("pWZ", i, i-2, sep="")] - mcmc[, paste("pWZ", i-j+3, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pWZ", i, i-2, "-pWZ", i-j+3, i-j+1, sep="")
            pest2A <- (pest2[paste("pWZ", i, i-2, sep="")] - pest2[paste("pWZ", i-j+3, i-j+1, sep="")])
            names(pest2A) <- paste("pWZ", i, i-2, "-pWZ", i-j+3, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)
          } # end (if W != "NULL")
        } # end (for i)
      } # end (for j)
    } # end (lag == 2)

    # -- Lag = 3 -- #
    if (lag == 3 & no.waves > 4) {
      for (j in 5:no.waves) {
        for (i in j:no.waves) {
          mcmcA <- (mcmc[, paste("pXY", i, i-3, sep="")] - mcmc[, paste("pXY", i-j+4, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pXY", i, i-3, "-pXY", i-j+4, i-j+1, sep="")
          pest2A <- (pest2[paste("pXY", i, i-3, sep="")] - pest2[paste("pXY", i-j+4, i-j+1, sep="")])
          names(pest2A) <- paste("pXY", i, i-3, "-pXY", i-j+4, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("pYX", i, i-3, sep="")] - mcmc[, paste("pYX", i-j+4, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pYX", i, i-3, "-pYX", i-j+4, i-j+1, sep="")
          pest2A <- (pest2[paste("pYX", i, i-3, sep="")] - pest2[paste("pYX", i-j+4, i-j+1, sep="")])
          names(pest2A) <- paste("pYX", i, i-3, "-pYX", i-j+4, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          if (Z != "NULL") {
            mcmcA <- (mcmc[, paste("pXZ", i, i-3, sep="")] - mcmc[, paste("pXZ", i-j+4, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pXZ", i, i-3, "-pXZ", i-j+4, i-j+1, sep="")
            pest2A <- (pest2[paste("pXZ", i, i-3, sep="")] - pest2[paste("pXZ", i-j+4, i-j+1, sep="")])
            names(pest2A) <- paste("pXZ", i, i-3, "-pXZ", i-j+4, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)

            mcmcA <- (mcmc[, paste("pYZ", i, i-3, sep="")] - mcmc[, paste("pYZ", i-j+4, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pYZ", i, i-3, "-pYZ", i-j+4, i-j+1, sep="")
            pest2A <- (pest2[paste("pYZ", i, i-3, sep="")] - pest2[paste("pYZ", i-j+4, i-j+1, sep="")])
            names(pest2A) <- paste("pYZ", i, i-3, "-pYZ", i-j+4, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)

            mcmcA <- (mcmc[, paste("pZX", i, i-3, sep="")] - mcmc[, paste("pZX", i-j+4, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pZX", i, i-3, "-pZX", i-j+4, i-j+1, sep="")
            pest2A <- (pest2[paste("pZX", i, i-3, sep="")] - pest2[paste("pZX", i-j+4, i-j+1, sep="")])
            names(pest2A) <- paste("pZX", i, i-3, "-pZX", i-j+4, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)

            mcmcA <- (mcmc[, paste("pZY", i, i-3, sep="")] - mcmc[, paste("pZY", i-j+4, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pZY", i, i-3, "-pZY", i-j+4, i-j+1, sep="")
            pest2A <- (pest2[paste("pZY", i, i-3, sep="")] - pest2[paste("pZY", i-j+4, i-j+1, sep="")])
            names(pest2A) <- paste("pZY", i, i-3, "-pZY", i-j+4, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)
          } # end (if Z != "NULL")

          if (W != "NULL") {
            mcmcA <- (mcmc[, paste("pXW", i, i-3, sep="")] - mcmc[, paste("pXW", i-j+4, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pXW", i, i-3, "-pXW", i-j+4, i-j+1, sep="")
            pest2A <- (pest2[paste("pXW", i, i-3, sep="")] - pest2[paste("pXW", i-j+4, i-j+1, sep="")])
            names(pest2A) <- paste("pXW", i, i-3, "-pXW", i-j+4, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)

            mcmcA <- (mcmc[, paste("pYW", i, i-3, sep="")] - mcmc[, paste("pYW", i-j+4, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pYW", i, i-3, "-pYW", i-j+4, i-j+1, sep="")
            pest2A <- (pest2[paste("pYW", i, i-3, sep="")] - pest2[paste("pYW", i-j+4, i-j+1, sep="")])
            names(pest2A) <- paste("pYW", i, i-3, "-pYW", i-j+4, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)

            mcmcA <- (mcmc[, paste("pZW", i, i-3, sep="")] - mcmc[, paste("pZW", i-j+4, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pZW", i, i-3, "-pZW", i-j+4, i-j+1, sep="")
            pest2A <- (pest2[paste("pZW", i, i-3, sep="")] - pest2[paste("pZW", i-j+4, i-j+1, sep="")])
            names(pest2A) <- paste("pZW", i, i-3, "-pZW", i-j+4, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)

            mcmcA <- (mcmc[, paste("pWX", i, i-3, sep="")] - mcmc[, paste("pWX", i-j+4, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pWX", i, i-3, "-pWX", i-j+4, i-j+1, sep="")
            pest2A <- (pest2[paste("pWX", i, i-3, sep="")] - pest2[paste("pWX", i-j+4, i-j+1, sep="")])
            names(pest2A) <- paste("pWX", i, i-3, "-pWX", i-j+4, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)

            mcmcA <- (mcmc[, paste("pWY", i, i-3, sep="")] - mcmc[, paste("pWY", i-j+4, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pWY", i, i-3, "-pWY", i-j+4, i-j+1, sep="")
            pest2A <- (pest2[paste("pWY", i, i-3, sep="")] - pest2[paste("pWY", i-j+4, i-j+1, sep="")])
            names(pest2A) <- paste("pWY", i, i-3, "-pWY", i-j+4, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)

            mcmcA <- (mcmc[, paste("pWZ", i, i-3, sep="")] - mcmc[, paste("pWZ", i-j+4, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pWZ", i, i-3, "-pWZ", i-j+4, i-j+1, sep="")
            pest2A <- (pest2[paste("pWZ", i, i-3, sep="")] - pest2[paste("pWZ", i-j+4, i-j+1, sep="")])
            names(pest2A) <- paste("pWZ", i, i-3, "-pWZ", i-j+4, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)
          } # end (if W != "NULL")
        } # end (for i)
      } # end (for j)
    } # end (lag == 3)

    # -- Lag = 4 -- #
    if (lag == 4 & no.waves > 5) {
      for (j in 6:no.waves) {
        for (i in j:no.waves) {
          mcmcA <- (mcmc[, paste("pXY", i, i-4, sep="")] - mcmc[, paste("pXY", i-j+5, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pXY", i, i-4, "-pXY", i-j+5, i-j+1, sep="")
          pest2A <- (pest2[paste("pXY", i, i-4, sep="")] - pest2[paste("pXY", i-j+5, i-j+1, sep="")])
          names(pest2A) <- paste("pXY", i, i-4, "-pXY", i-j+5, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("pYX", i, i-4, sep="")] - mcmc[, paste("pYX", i-j+5, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pYX", i, i-4, "-pYX", i-j+5, i-j+1, sep="")
          pest2A <- (pest2[paste("pYX", i, i-4, sep="")] - pest2[paste("pYX", i-j+5, i-j+1, sep="")])
          names(pest2A) <- paste("pYX", i, i-4, "-pYX", i-j+5, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          if (Z != "NULL") {
            mcmcA <- (mcmc[, paste("pXZ", i, i-4, sep="")] - mcmc[, paste("pXZ", i-j+5, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pXZ", i, i-4, "-pXZ", i-j+5, i-j+1, sep="")
            pest2A <- (pest2[paste("pXZ", i, i-4, sep="")] - pest2[paste("pXZ", i-j+5, i-j+1, sep="")])
            names(pest2A) <- paste("pXZ", i, i-4, "-pXZ", i-j+5, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)

            mcmcA <- (mcmc[, paste("pYZ", i, i-4, sep="")] - mcmc[, paste("pYZ", i-j+5, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pYZ", i, i-4, "-pYZ", i-j+5, i-j+1, sep="")
            pest2A <- (pest2[paste("pYZ", i, i-4, sep="")] - pest2[paste("pYZ", i-j+5, i-j+1, sep="")])
            names(pest2A) <- paste("pYZ", i, i-4, "-pYZ", i-j+5, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)

            mcmcA <- (mcmc[, paste("pZX", i, i-4, sep="")] - mcmc[, paste("pZX", i-j+5, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pZX", i, i-4, "-pZX", i-j+5, i-j+1, sep="")
            pest2A <- (pest2[paste("pZX", i, i-4, sep="")] - pest2[paste("pZX", i-j+5, i-j+1, sep="")])
            names(pest2A) <- paste("pZX", i, i-4, "-pZX", i-j+5, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)

            mcmcA <- (mcmc[, paste("pZY", i, i-4, sep="")] - mcmc[, paste("pZY", i-j+5, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pZY", i, i-4, "-pZY", i-j+5, i-j+1, sep="")
            pest2A <- (pest2[paste("pZY", i, i-4, sep="")] - pest2[paste("pZY", i-j+5, i-j+1, sep="")])
            names(pest2A) <- paste("pZY", i, i-4, "-pZY", i-j+5, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)
          } # end (if Z != "NULL")

          if (W != "NULL") {
            mcmcA <- (mcmc[, paste("pXW", i, i-4, sep="")] - mcmc[, paste("pXW", i-j+5, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pXW", i, i-4, "-pXW", i-j+5, i-j+1, sep="")
            pest2A <- (pest2[paste("pXW", i, i-4, sep="")] - pest2[paste("pXW", i-j+5, i-j+1, sep="")])
            names(pest2A) <- paste("pXW", i, i-4, "-pXW", i-j+5, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)

            mcmcA <- (mcmc[, paste("pYW", i, i-4, sep="")] - mcmc[, paste("pYW", i-j+5, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pYW", i, i-4, "-pYW", i-j+5, i-j+1, sep="")
            pest2A <- (pest2[paste("pYW", i, i-4, sep="")] - pest2[paste("pYW", i-j+5, i-j+1, sep="")])
            names(pest2A) <- paste("pYW", i, i-4, "-pYW", i-j+5, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)

            mcmcA <- (mcmc[, paste("pZW", i, i-4, sep="")] - mcmc[, paste("pZW", i-j+5, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pZW", i, i-4, "-pZW", i-j+5, i-j+1, sep="")
            pest2A <- (pest2[paste("pZW", i, i-4, sep="")] - pest2[paste("pZW", i-j+5, i-j+1, sep="")])
            names(pest2A) <- paste("pZW", i, i-4, "-pZW", i-j+5, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)

            mcmcA <- (mcmc[, paste("pWX", i, i-4, sep="")] - mcmc[, paste("pWX", i-j+5, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pWX", i, i-4, "-pWX", i-j+5, i-j+1, sep="")
            pest2A <- (pest2[paste("pWX", i, i-4, sep="")] - pest2[paste("pWX", i-j+5, i-j+1, sep="")])
            names(pest2A) <- paste("pWX", i, i-4, "-pWX", i-j+5, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)

            mcmcA <- (mcmc[, paste("pWY", i, i-4, sep="")] - mcmc[, paste("pWY", i-j+5, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pWY", i, i-4, "-pWY", i-j+5, i-j+1, sep="")
            pest2A <- (pest2[paste("pWY", i, i-4, sep="")] - pest2[paste("pWY", i-j+5, i-j+1, sep="")])
            names(pest2A) <- paste("pWY", i, i-4, "-pWY", i-j+5, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)

            mcmcA <- (mcmc[, paste("pWZ", i, i-4, sep="")] - mcmc[, paste("pWZ", i-j+5, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("pWZ", i, i-4, "-pWZ", i-j+5, i-j+1, sep="")
            pest2A <- (pest2[paste("pWZ", i, i-4, sep="")] - pest2[paste("pWZ", i-j+5, i-j+1, sep="")])
            names(pest2A) <- paste("pWZ", i, i-4, "-pWZ", i-j+5, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)
          } # end (if W != "NULL")
        } # end (for i)
      } # end (for j)
    } # end (lag == 4)
  } # end (if pXY21)
  ## ----- end (Difference in path coefficients) ----- ##


  ## -- Differences in proportional change dXX -- ##
  if (any(parEst[,4] == "dXX21")) {

    # -- Lag = 1 -- #
    for (j in 3:no.waves) {
      for (i in j:no.waves) {
        mcmcA <- (mcmc[, paste("dXX", i, i-1, sep="")] - mcmc[, paste("dXX", i-j+2, i-j+1, sep="")])
        mcmc <- cbind(mcmc, mcmcA)
        colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dXX", i, i-1, "-dXX", i-j+2, i-j+1, sep="")
        pest2A <- (pest2[paste("dXX", i, i-1, sep="")] - pest2[paste("dXX", i-j+2, i-j+1, sep="")])
        names(pest2A) <- paste("dXX", i, i-1, "-dXX", i-j+2, i-j+1, sep="")
        pest2 <- append(pest2, pest2A)

        mcmcA <- (mcmc[, paste("dYY", i, i-1, sep="")] - mcmc[, paste("dYY", i-j+2, i-j+1, sep="")])
        mcmc <- cbind(mcmc, mcmcA)
        colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dYY", i, i-1, "-dYY", i-j+2, i-j+1, sep="")
        pest2A <- (pest2[paste("dYY", i, i-1, sep="")] - pest2[paste("dYY", i-j+2, i-j+1, sep="")])
        names(pest2A) <- paste("dYY", i, i-1, "-dYY", i-j+2, i-j+1, sep="")
        pest2 <- append(pest2, pest2A)

        if (Z != "NULL") {
          mcmcA <- (mcmc[, paste("dZZ", i, i-1, sep="")] - mcmc[, paste("dZZ", i-j+2, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dZZ", i, i-1, "-dZZ", i-j+2, i-j+1, sep="")
          pest2A <- (pest2[paste("dZZ", i, i-1, sep="")] - pest2[paste("dZZ", i-j+2, i-j+1, sep="")])
          names(pest2A) <- paste("dZZ", i, i-1, "-dZZ", i-j+2, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)
        } # end (if Z)

        if (W != "NULL") {
          mcmcA <- (mcmc[, paste("dWW", i, i-1, sep="")] - mcmc[, paste("dWW", i-j+2, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dWW", i, i-1, "-dWW", i-j+2, i-j+1, sep="")
          pest2A <- (pest2[paste("dWW", i, i-1, sep="")] - pest2[paste("dWW", i-j+2, i-j+1, sep="")])
          names(pest2A) <- paste("dWW", i, i-1, "-dWW", i-j+2, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)
        } # end (if W != "NULL")
      } # end (for i)
    } # end (for j)

    # -- Lag = 2 -- #
    if (lag == 2 & no.waves > 3) {
      for (j in 4:no.waves) {
        for (i in j:no.waves) {
          mcmcA <- (mcmc[, paste("dXX", i, i-2, sep="")] - mcmc[, paste("dXX", i-j+3, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dXX", i, i-2, "-dXX", i-j+3, i-j+1, sep="")
          pest2A <- (pest2[paste("dXX", i, i-2, sep="")] - pest2[paste("dXX", i-j+3, i-j+1, sep="")])
          names(pest2A) <- paste("dXX", i, i-2, "-dXX", i-j+3, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("dYY", i, i-2, sep="")] - mcmc[, paste("dYY", i-j+3, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dYY", i, i-2, "-dYY", i-j+3, i-j+1, sep="")
          pest2A <- (pest2[paste("dYY", i, i-2, sep="")] - pest2[paste("dYY", i-j+3, i-j+1, sep="")])
          names(pest2A) <- paste("dYY", i, i-2, "-dYY", i-j+3, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          if (Z != "NULL") {
            mcmcA <- (mcmc[, paste("dZZ", i, i-2, sep="")] - mcmc[, paste("dZZ", i-j+3, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dZZ", i, i-2, "-dZZ", i-j+3, i-j+1, sep="")
            pest2A <- (pest2[paste("dZZ", i, i-2, sep="")] - pest2[paste("dZZ", i-j+3, i-j+1, sep="")])
            names(pest2A) <- paste("dZZ", i, i-2, "-dZZ", i-j+3, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)
          } # end (if Z != "NULL")

          if (W != "NULL") {
            mcmcA <- (mcmc[, paste("dWW", i, i-2, sep="")] - mcmc[, paste("dWW", i-j+3, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dWW", i, i-2, "-dWW", i-j+3, i-j+1, sep="")
            pest2A <- (pest2[paste("dWW", i, i-2, sep="")] - pest2[paste("dWW", i-j+3, i-j+1, sep="")])
            names(pest2A) <- paste("dWW", i, i-2, "-dWW", i-j+3, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)
          } # end (if W != "NULL")
        } # end (for i)
      } # end (for j)
    } # end (lag == 2)

    # -- Lag = 3 -- #
    if (lag == 3 & no.waves > 4) {
      for (j in 5:no.waves) {
        for (i in j:no.waves) {
          mcmcA <- (mcmc[, paste("dXX", i, i-3, sep="")] - mcmc[, paste("dXX", i-j+4, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dXX", i, i-3, "-dXX", i-j+4, i-j+1, sep="")
          pest2A <- (pest2[paste("dXX", i, i-3, sep="")] - pest2[paste("dXX", i-j+4, i-j+1, sep="")])
          names(pest2A) <- paste("dXX", i, i-3, "-dXX", i-j+4, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("dYY", i, i-3, sep="")] - mcmc[, paste("dYY", i-j+4, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dYY", i, i-3, "-dYY", i-j+4, i-j+1, sep="")
          pest2A <- (pest2[paste("dYY", i, i-3, sep="")] - pest2[paste("dYY", i-j+4, i-j+1, sep="")])
          names(pest2A) <- paste("dYY", i, i-3, "-dYY", i-j+4, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          if (Z != "NULL") {
            mcmcA <- (mcmc[, paste("dZZ", i, i-3, sep="")] - mcmc[, paste("dZZ", i-j+4, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dZZ", i, i-3, "-dZZ", i-j+4, i-j+1, sep="")
            pest2A <- (pest2[paste("dZZ", i, i-3, sep="")] - pest2[paste("dZZ", i-j+4, i-j+1, sep="")])
            names(pest2A) <- paste("dZZ", i, i-3, "-dZZ", i-j+4, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)
          } # end (if Z != "NULL")

          if (W != "NULL") {
            mcmcA <- (mcmc[, paste("dWW", i, i-3, sep="")] - mcmc[, paste("dWW", i-j+4, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dWW", i, i-3, "-dWW", i-j+4, i-j+1, sep="")
            pest2A <- (pest2[paste("dWW", i, i-3, sep="")] - pest2[paste("dWW", i-j+4, i-j+1, sep="")])
            names(pest2A) <- paste("dWW", i, i-3, "-dWW", i-j+4, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)
          } # end (if W != "NULL")
        } # end (for i)
      } # end (for j)
    } # end (lag == 3)

    # -- Lag = 4 -- #
    if (lag == 4 & no.waves > 5) {
      for (j in 6:no.waves) {
        for (i in j:no.waves) {
          mcmcA <- (mcmc[, paste("dXX", i, i-4, sep="")] - mcmc[, paste("dXX", i-j+5, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dXX", i, i-4, "-dXX", i-j+5, i-j+1, sep="")
          pest2A <- (pest2[paste("dXX", i, i-4, sep="")] - pest2[paste("dXX", i-j+5, i-j+1, sep="")])
          names(pest2A) <- paste("dXX", i, i-4, "-dXX", i-j+5, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("dYY", i, i-4, sep="")] - mcmc[, paste("dYY", i-j+5, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dYY", i, i-4, "-dYY", i-j+5, i-j+1, sep="")
          pest2A <- (pest2[paste("dYY", i, i-4, sep="")] - pest2[paste("dYY", i-j+5, i-j+1, sep="")])
          names(pest2A) <- paste("dYY", i, i-4, "-dYY", i-j+5, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          if (Z != "NULL") {
            mcmcA <- (mcmc[, paste("dZZ", i, i-4, sep="")] - mcmc[, paste("dZZ", i-j+5, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dZZ", i, i-4, "-dZZ", i-j+5, i-j+1, sep="")
            pest2A <- (pest2[paste("dZZ", i, i-4, sep="")] - pest2[paste("dZZ", i-j+5, i-j+1, sep="")])
            names(pest2A) <- paste("dZZ", i, i-4, "-dZZ", i-j+5, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)
          } # end (if Z != "NULL")

          if (W != "NULL") {
            mcmcA <- (mcmc[, paste("dWW", i, i-4, sep="")] - mcmc[, paste("dWW", i-j+5, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dWW", i, i-4, "-dWW", i-j+5, i-j+1, sep="")
            pest2A <- (pest2[paste("dWW", i, i-4, sep="")] - pest2[paste("dWW", i-j+5, i-j+1, sep="")])
            names(pest2A) <- paste("dWW", i, i-4, "-dWW", i-j+5, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)
          } # end (if W != "NULL")
        } # end (for i)
      } # end (for j)
    } # end (lag == 4)
  } # end (if dXX21)
  ## ----- end (Difference in proportional change) ----- ##


  ## -- Differences in proportional change dXY -- ##
  if (any(parEst[,4] == "dXY21")) {

    # -- Lag = 1 -- #
    for (j in 3:no.waves) {
      for (i in j:no.waves) {
        mcmcA <- (mcmc[, paste("dXY", i, i-1, sep="")] - mcmc[, paste("dXY", i-j+2, i-j+1, sep="")])
        mcmc <- cbind(mcmc, mcmcA)
        colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dXY", i, i-1, "-dXY", i-j+2, i-j+1, sep="")
        pest2A <- (pest2[paste("dXY", i, i-1, sep="")] - pest2[paste("dXY", i-j+2, i-j+1, sep="")])
        names(pest2A) <- paste("dXY", i, i-1, "-dXY", i-j+2, i-j+1, sep="")
        pest2 <- append(pest2, pest2A)

        mcmcA <- (mcmc[, paste("dYX", i, i-1, sep="")] - mcmc[, paste("dYX", i-j+2, i-j+1, sep="")])
        mcmc <- cbind(mcmc, mcmcA)
        colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dYX", i, i-1, "-dYX", i-j+2, i-j+1, sep="")
        pest2A <- (pest2[paste("dYX", i, i-1, sep="")] - pest2[paste("dYX", i-j+2, i-j+1, sep="")])
        names(pest2A) <- paste("dYX", i, i-1, "-dYX", i-j+2, i-j+1, sep="")
        pest2 <- append(pest2, pest2A)

        if (Z != "NULL") {
          mcmcA <- (mcmc[, paste("dXZ", i, i-1, sep="")] - mcmc[, paste("dXZ", i-j+2, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dXZ", i, i-1, "-dXZ", i-j+2, i-j+1, sep="")
          pest2A <- (pest2[paste("dXZ", i, i-1, sep="")] - pest2[paste("dXZ", i-j+2, i-j+1, sep="")])
          names(pest2A) <- paste("dXZ", i, i-1, "-dXZ", i-j+2, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("dYZ", i, i-1, sep="")] - mcmc[, paste("dYZ", i-j+2, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dYZ", i, i-1, "-dYZ", i-j+2, i-j+1, sep="")
          pest2A <- (pest2[paste("dYZ", i, i-1, sep="")] - pest2[paste("dYZ", i-j+2, i-j+1, sep="")])
          names(pest2A) <- paste("dYZ", i, i-1, "-dYZ", i-j+2, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("dZX", i, i-1, sep="")] - mcmc[, paste("dZX", i-j+2, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dZX", i, i-1, "-dZX", i-j+2, i-j+1, sep="")
          pest2A <- (pest2[paste("dZX", i, i-1, sep="")] - pest2[paste("dZX", i-j+2, i-j+1, sep="")])
          names(pest2A) <- paste("dZX", i, i-1, "-dZX", i-j+2, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("dZY", i, i-1, sep="")] - mcmc[, paste("dZY", i-j+2, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dZY", i, i-1, "-dZY", i-j+2, i-j+1, sep="")
          pest2A <- (pest2[paste("dZY", i, i-1, sep="")] - pest2[paste("dZY", i-j+2, i-j+1, sep="")])
          names(pest2A) <- paste("dZY", i, i-1, "-dZY", i-j+2, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)
        } # end (if Z != "NULL")

        if (W != "NULL") {
          mcmcA <- (mcmc[, paste("dXW", i, i-1, sep="")] - mcmc[, paste("dXW", i-j+2, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dXW", i, i-1, "-dXW", i-j+2, i-j+1, sep="")
          pest2A <- (pest2[paste("dXW", i, i-1, sep="")] - pest2[paste("dXW", i-j+2, i-j+1, sep="")])
          names(pest2A) <- paste("dXW", i, i-1, "-dXW", i-j+2, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("dYW", i, i-1, sep="")] - mcmc[, paste("dYW", i-j+2, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dYW", i, i-1, "-dYW", i-j+2, i-j+1, sep="")
          pest2A <- (pest2[paste("dYW", i, i-1, sep="")] - pest2[paste("dYW", i-j+2, i-j+1, sep="")])
          names(pest2A) <- paste("dYW", i, i-1, "-dYW", i-j+2, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("dZW", i, i-1, sep="")] - mcmc[, paste("dZW", i-j+2, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dZW", i, i-1, "-dZW", i-j+2, i-j+1, sep="")
          pest2A <- (pest2[paste("dZW", i, i-1, sep="")] - pest2[paste("dZW", i-j+2, i-j+1, sep="")])
          names(pest2A) <- paste("dZW", i, i-1, "-dZW", i-j+2, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("dWX", i, i-1, sep="")] - mcmc[, paste("dWX", i-j+2, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dWX", i, i-1, "-dWX", i-j+2, i-j+1, sep="")
          pest2A <- (pest2[paste("dWX", i, i-1, sep="")] - pest2[paste("dWX", i-j+2, i-j+1, sep="")])
          names(pest2A) <- paste("dWX", i, i-1, "-dWX", i-j+2, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("dWY", i, i-1, sep="")] - mcmc[, paste("dWY", i-j+2, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dWY", i, i-1, "-dWY", i-j+2, i-j+1, sep="")
          pest2A <- (pest2[paste("dWY", i, i-1, sep="")] - pest2[paste("dWY", i-j+2, i-j+1, sep="")])
          names(pest2A) <- paste("dWY", i, i-1, "-dWY", i-j+2, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("dWZ", i, i-1, sep="")] - mcmc[, paste("dWZ", i-j+2, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dWZ", i, i-1, "-dWZ", i-j+2, i-j+1, sep="")
          pest2A <- (pest2[paste("dWZ", i, i-1, sep="")] - pest2[paste("dWZ", i-j+2, i-j+1, sep="")])
          names(pest2A) <- paste("dWZ", i, i-1, "-dWZ", i-j+2, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)
        } # end (if W != "NULL")
      } # end (for i)
    } # end (for j)

    # -- Lag = 2 -- #
    if (lag == 2 & no.waves > 3) {
      for (j in 4:no.waves) {
        for (i in j:no.waves) {
          mcmcA <- (mcmc[, paste("dXY", i, i-2, sep="")] - mcmc[, paste("dXY", i-j+3, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dXY", i, i-2, "-dXY", i-j+3, i-j+1, sep="")
          pest2A <- (pest2[paste("dXY", i, i-2, sep="")] - pest2[paste("dXY", i-j+3, i-j+1, sep="")])
          names(pest2A) <- paste("dXY", i, i-2, "-dXY", i-j+3, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("dYX", i, i-2, sep="")] - mcmc[, paste("dYX", i-j+3, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dYX", i, i-2, "-dYX", i-j+3, i-j+1, sep="")
          pest2A <- (pest2[paste("dYX", i, i-2, sep="")] - pest2[paste("dYX", i-j+3, i-j+1, sep="")])
          names(pest2A) <- paste("dYX", i, i-2, "-dYX", i-j+3, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          if (Z != "NULL") {
            mcmcA <- (mcmc[, paste("dXZ", i, i-2, sep="")] - mcmc[, paste("dXZ", i-j+3, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dXZ", i, i-2, "-dXZ", i-j+3, i-j+1, sep="")
            pest2A <- (pest2[paste("dXZ", i, i-2, sep="")] - pest2[paste("dXZ", i-j+3, i-j+1, sep="")])
            names(pest2A) <- paste("dXZ", i, i-2, "-dXZ", i-j+3, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)

            mcmcA <- (mcmc[, paste("dYZ", i, i-2, sep="")] - mcmc[, paste("dYZ", i-j+3, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dYZ", i, i-2, "-dYZ", i-j+3, i-j+1, sep="")
            pest2A <- (pest2[paste("dYZ", i, i-2, sep="")] - pest2[paste("dYZ", i-j+3, i-j+1, sep="")])
            names(pest2A) <- paste("dYZ", i, i-2, "-dYZ", i-j+3, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)

            mcmcA <- (mcmc[, paste("dZX", i, i-2, sep="")] - mcmc[, paste("dZX", i-j+3, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dZX", i, i-2, "-dZX", i-j+3, i-j+1, sep="")
            pest2A <- (pest2[paste("dZX", i, i-2, sep="")] - pest2[paste("dZX", i-j+3, i-j+1, sep="")])
            names(pest2A) <- paste("dZX", i, i-2, "-dZX", i-j+3, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)

            mcmcA <- (mcmc[, paste("dZY", i, i-2, sep="")] - mcmc[, paste("dZY", i-j+3, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dZY", i, i-2, "-dZY", i-j+3, i-j+1, sep="")
            pest2A <- (pest2[paste("dZY", i, i-2, sep="")] - pest2[paste("dZY", i-j+3, i-j+1, sep="")])
            names(pest2A) <- paste("dZY", i, i-2, "-dZY", i-j+3, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)
          } # end (if Z != "NULL")

          if (W != "NULL") {
            mcmcA <- (mcmc[, paste("dXW", i, i-2, sep="")] - mcmc[, paste("dXW", i-j+3, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dXW", i, i-2, "-dXW", i-j+3, i-j+1, sep="")
            pest2A <- (pest2[paste("dXW", i, i-2, sep="")] - pest2[paste("dXW", i-j+3, i-j+1, sep="")])
            names(pest2A) <- paste("dXW", i, i-2, "-dXW", i-j+3, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)

            mcmcA <- (mcmc[, paste("dYW", i, i-2, sep="")] - mcmc[, paste("dYW", i-j+3, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dYW", i, i-2, "-dYW", i-j+3, i-j+1, sep="")
            pest2A <- (pest2[paste("dYW", i, i-2, sep="")] - pest2[paste("dYW", i-j+3, i-j+1, sep="")])
            names(pest2A) <- paste("dYW", i, i-2, "-dYW", i-j+3, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)

            mcmcA <- (mcmc[, paste("dZW", i, i-2, sep="")] - mcmc[, paste("dZW", i-j+3, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dZW", i, i-2, "-dZW", i-j+3, i-j+1, sep="")
            pest2A <- (pest2[paste("dZW", i, i-2, sep="")] - pest2[paste("dZW", i-j+3, i-j+1, sep="")])
            names(pest2A) <- paste("dZW", i, i-2, "-dZW", i-j+3, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)

            mcmcA <- (mcmc[, paste("dWX", i, i-2, sep="")] - mcmc[, paste("dWX", i-j+3, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dWX", i, i-2, "-dWX", i-j+3, i-j+1, sep="")
            pest2A <- (pest2[paste("dWX", i, i-2, sep="")] - pest2[paste("dWX", i-j+3, i-j+1, sep="")])
            names(pest2A) <- paste("dWX", i, i-2, "-dWX", i-j+3, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)

            mcmcA <- (mcmc[, paste("dWY", i, i-2, sep="")] - mcmc[, paste("dWY", i-j+3, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dWY", i, i-2, "-dWY", i-j+3, i-j+1, sep="")
            pest2A <- (pest2[paste("dWY", i, i-2, sep="")] - pest2[paste("dWY", i-j+3, i-j+1, sep="")])
            names(pest2A) <- paste("dWY", i, i-2, "-dWY", i-j+3, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)

            mcmcA <- (mcmc[, paste("dWZ", i, i-2, sep="")] - mcmc[, paste("dWZ", i-j+3, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dWZ", i, i-2, "-dWZ", i-j+3, i-j+1, sep="")
            pest2A <- (pest2[paste("dWZ", i, i-2, sep="")] - pest2[paste("dWZ", i-j+3, i-j+1, sep="")])
            names(pest2A) <- paste("dWZ", i, i-2, "-dWZ", i-j+3, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)
          } # end (if W != "NULL")
        } # end (for i)
      } # end (for j)
    } # end (lag == 2)

    # -- Lag = 3 -- #
    if (lag == 3 & no.waves > 4) {
      for (j in 5:no.waves) {
        for (i in j:no.waves) {
          mcmcA <- (mcmc[, paste("dXY", i, i-3, sep="")] - mcmc[, paste("dXY", i-j+4, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dXY", i, i-3, "-dXY", i-j+4, i-j+1, sep="")
          pest2A <- (pest2[paste("dXY", i, i-3, sep="")] - pest2[paste("dXY", i-j+4, i-j+1, sep="")])
          names(pest2A) <- paste("dXY", i, i-3, "-dXY", i-j+4, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("dYX", i, i-3, sep="")] - mcmc[, paste("dYX", i-j+4, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dYX", i, i-3, "-dYX", i-j+4, i-j+1, sep="")
          pest2A <- (pest2[paste("dYX", i, i-3, sep="")] - pest2[paste("dYX", i-j+4, i-j+1, sep="")])
          names(pest2A) <- paste("dYX", i, i-3, "-dYX", i-j+4, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          if (Z != "NULL") {
            mcmcA <- (mcmc[, paste("dXZ", i, i-3, sep="")] - mcmc[, paste("dXZ", i-j+4, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dXZ", i, i-3, "-dXZ", i-j+4, i-j+1, sep="")
            pest2A <- (pest2[paste("dXZ", i, i-3, sep="")] - pest2[paste("dXZ", i-j+4, i-j+1, sep="")])
            names(pest2A) <- paste("dXZ", i, i-3, "-dXZ", i-j+4, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)

            mcmcA <- (mcmc[, paste("dYZ", i, i-3, sep="")] - mcmc[, paste("dYZ", i-j+4, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dYZ", i, i-3, "-dYZ", i-j+4, i-j+1, sep="")
            pest2A <- (pest2[paste("dYZ", i, i-3, sep="")] - pest2[paste("dYZ", i-j+4, i-j+1, sep="")])
            names(pest2A) <- paste("dYZ", i, i-3, "-dYZ", i-j+4, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)

            mcmcA <- (mcmc[, paste("dZX", i, i-3, sep="")] - mcmc[, paste("dZX", i-j+4, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dZX", i, i-3, "-dZX", i-j+4, i-j+1, sep="")
            pest2A <- (pest2[paste("dZX", i, i-3, sep="")] - pest2[paste("dZX", i-j+4, i-j+1, sep="")])
            names(pest2A) <- paste("dZX", i, i-3, "-dZX", i-j+4, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)

            mcmcA <- (mcmc[, paste("dZY", i, i-3, sep="")] - mcmc[, paste("dZY", i-j+4, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dZY", i, i-3, "-dZY", i-j+4, i-j+1, sep="")
            pest2A <- (pest2[paste("dZY", i, i-3, sep="")] - pest2[paste("dZY", i-j+4, i-j+1, sep="")])
            names(pest2A) <- paste("dZY", i, i-3, "-dZY", i-j+4, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)
          } # end (if Z != "NULL")

          if (W != "NULL") {
            mcmcA <- (mcmc[, paste("dXW", i, i-3, sep="")] - mcmc[, paste("dXW", i-j+4, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dXW", i, i-3, "-dXW", i-j+4, i-j+1, sep="")
            pest2A <- (pest2[paste("dXW", i, i-3, sep="")] - pest2[paste("dXW", i-j+4, i-j+1, sep="")])
            names(pest2A) <- paste("dXW", i, i-3, "-dXW", i-j+4, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)

            mcmcA <- (mcmc[, paste("dYW", i, i-3, sep="")] - mcmc[, paste("dYW", i-j+4, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dYW", i, i-3, "-dYW", i-j+4, i-j+1, sep="")
            pest2A <- (pest2[paste("dYW", i, i-3, sep="")] - pest2[paste("dYW", i-j+4, i-j+1, sep="")])
            names(pest2A) <- paste("dYW", i, i-3, "-dYW", i-j+4, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)

            mcmcA <- (mcmc[, paste("dZW", i, i-3, sep="")] - mcmc[, paste("dZW", i-j+4, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dZW", i, i-3, "-dZW", i-j+4, i-j+1, sep="")
            pest2A <- (pest2[paste("dZW", i, i-3, sep="")] - pest2[paste("dZW", i-j+4, i-j+1, sep="")])
            names(pest2A) <- paste("dZW", i, i-3, "-dZW", i-j+4, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)

            mcmcA <- (mcmc[, paste("dWX", i, i-3, sep="")] - mcmc[, paste("dWX", i-j+4, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dWX", i, i-3, "-dWX", i-j+4, i-j+1, sep="")
            pest2A <- (pest2[paste("dWX", i, i-3, sep="")] - pest2[paste("dWX", i-j+4, i-j+1, sep="")])
            names(pest2A) <- paste("dWX", i, i-3, "-dWX", i-j+4, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)

            mcmcA <- (mcmc[, paste("dWY", i, i-3, sep="")] - mcmc[, paste("dWY", i-j+4, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dWY", i, i-3, "-dWY", i-j+4, i-j+1, sep="")
            pest2A <- (pest2[paste("dWY", i, i-3, sep="")] - pest2[paste("dWY", i-j+4, i-j+1, sep="")])
            names(pest2A) <- paste("dWY", i, i-3, "-dWY", i-j+4, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)

            mcmcA <- (mcmc[, paste("dWZ", i, i-3, sep="")] - mcmc[, paste("dWZ", i-j+4, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dWZ", i, i-3, "-dWZ", i-j+4, i-j+1, sep="")
            pest2A <- (pest2[paste("dWZ", i, i-3, sep="")] - pest2[paste("dWZ", i-j+4, i-j+1, sep="")])
            names(pest2A) <- paste("dWZ", i, i-3, "-dWZ", i-j+4, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)
          } # end (if W != "NULL")
        } # end (for i)
      } # end (for j)
    } # end (lag == 3)

    # -- Lag = 4 -- #
    if (lag == 4 & no.waves > 5) {
      for (j in 6:no.waves) {
        for (i in j:no.waves) {
          mcmcA <- (mcmc[, paste("dXY", i, i-4, sep="")] - mcmc[, paste("dXY", i-j+5, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dXY", i, i-4, "-dXY", i-j+5, i-j+1, sep="")
          pest2A <- (pest2[paste("dXY", i, i-4, sep="")] - pest2[paste("dXY", i-j+5, i-j+1, sep="")])
          names(pest2A) <- paste("dXY", i, i-4, "-dXY", i-j+5, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("dYX", i, i-4, sep="")] - mcmc[, paste("dYX", i-j+5, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dYX", i, i-4, "-dYX", i-j+5, i-j+1, sep="")
          pest2A <- (pest2[paste("dYX", i, i-4, sep="")] - pest2[paste("dYX", i-j+5, i-j+1, sep="")])
          names(pest2A) <- paste("dYX", i, i-4, "-dYX", i-j+5, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          if (Z != "NULL") {
            mcmcA <- (mcmc[, paste("dXZ", i, i-4, sep="")] - mcmc[, paste("dXZ", i-j+5, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dXZ", i, i-4, "-dXZ", i-j+5, i-j+1, sep="")
            pest2A <- (pest2[paste("dXZ", i, i-4, sep="")] - pest2[paste("dXZ", i-j+5, i-j+1, sep="")])
            names(pest2A) <- paste("dXZ", i, i-4, "-dXZ", i-j+5, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)

            mcmcA <- (mcmc[, paste("dYZ", i, i-4, sep="")] - mcmc[, paste("dYZ", i-j+5, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dYZ", i, i-4, "-dYZ", i-j+5, i-j+1, sep="")
            pest2A <- (pest2[paste("dYZ", i, i-4, sep="")] - pest2[paste("dYZ", i-j+5, i-j+1, sep="")])
            names(pest2A) <- paste("dYZ", i, i-4, "-dYZ", i-j+5, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)

            mcmcA <- (mcmc[, paste("dZX", i, i-4, sep="")] - mcmc[, paste("dZX", i-j+5, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dZX", i, i-4, "-dZX", i-j+5, i-j+1, sep="")
            pest2A <- (pest2[paste("dZX", i, i-4, sep="")] - pest2[paste("dZX", i-j+5, i-j+1, sep="")])
            names(pest2A) <- paste("dZX", i, i-4, "-dZX", i-j+5, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)

            mcmcA <- (mcmc[, paste("dZY", i, i-4, sep="")] - mcmc[, paste("dZY", i-j+5, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dZY", i, i-4, "-dZY", i-j+5, i-j+1, sep="")
            pest2A <- (pest2[paste("dZY", i, i-4, sep="")] - pest2[paste("dZY", i-j+5, i-j+1, sep="")])
            names(pest2A) <- paste("dZY", i, i-4, "-dZY", i-j+5, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)
          } # end (if Z != "NULL")

          if (W != "NULL") {
            mcmcA <- (mcmc[, paste("dXW", i, i-4, sep="")] - mcmc[, paste("dXW", i-j+5, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dXW", i, i-4, "-dXW", i-j+5, i-j+1, sep="")
            pest2A <- (pest2[paste("dXW", i, i-4, sep="")] - pest2[paste("dXW", i-j+5, i-j+1, sep="")])
            names(pest2A) <- paste("dXW", i, i-4, "-dXW", i-j+5, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)

            mcmcA <- (mcmc[, paste("dYW", i, i-4, sep="")] - mcmc[, paste("dYW", i-j+5, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dYW", i, i-4, "-dYW", i-j+5, i-j+1, sep="")
            pest2A <- (pest2[paste("dYW", i, i-4, sep="")] - pest2[paste("dYW", i-j+5, i-j+1, sep="")])
            names(pest2A) <- paste("dYW", i, i-4, "-dYW", i-j+5, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)

            mcmcA <- (mcmc[, paste("dZW", i, i-4, sep="")] - mcmc[, paste("dZW", i-j+5, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dZW", i, i-4, "-dZW", i-j+5, i-j+1, sep="")
            pest2A <- (pest2[paste("dZW", i, i-4, sep="")] - pest2[paste("dZW", i-j+5, i-j+1, sep="")])
            names(pest2A) <- paste("dZW", i, i-4, "-dZW", i-j+5, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)

            mcmcA <- (mcmc[, paste("dWX", i, i-4, sep="")] - mcmc[, paste("dWX", i-j+5, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dWX", i, i-4, "-dWX", i-j+5, i-j+1, sep="")
            pest2A <- (pest2[paste("dWX", i, i-4, sep="")] - pest2[paste("dWX", i-j+5, i-j+1, sep="")])
            names(pest2A) <- paste("dWX", i, i-4, "-dWX", i-j+5, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)

            mcmcA <- (mcmc[, paste("dWY", i, i-4, sep="")] - mcmc[, paste("dWY", i-j+5, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dWY", i, i-4, "-dWY", i-j+5, i-j+1, sep="")
            pest2A <- (pest2[paste("dWY", i, i-4, sep="")] - pest2[paste("dWY", i-j+5, i-j+1, sep="")])
            names(pest2A) <- paste("dWY", i, i-4, "-dWY", i-j+5, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)

            mcmcA <- (mcmc[, paste("dWZ", i, i-4, sep="")] - mcmc[, paste("dWZ", i-j+5, i-j+1, sep="")])
            mcmc <- cbind(mcmc, mcmcA)
            colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("dWZ", i, i-4, "-dWZ", i-j+5, i-j+1, sep="")
            pest2A <- (pest2[paste("dWZ", i, i-4, sep="")] - pest2[paste("dWZ", i-j+5, i-j+1, sep="")])
            names(pest2A) <- paste("dWZ", i, i-4, "-dWZ", i-j+5, i-j+1, sep="")
            pest2 <- append(pest2, pest2A)
          } # end (if W != "NULL")
        } # end (for i)
      } # end (for j)
    } # end (lag == 4)
  } # end (if dXY21)
  ## ----- end (Difference in proportional change) ----- ##


  ## -- Differences in covariance of latent variable residuals -- ##
  if (any(parEst[,4] == "eXY2")) {
    for (j in 3:no.waves) {
      for (i in j:no.waves) {
        mcmcA <- (mcmc[, paste("eXY", i, sep="")] - mcmc[, paste("eXY", i-j+2, sep="")])
        mcmc <- cbind(mcmc, mcmcA)
        colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("eXY", i, "-eXY", i-j+2, sep="")
        pest2A <- (pest2[paste("eXY", i, sep="")] - pest2[paste("eXY", i-j+2, sep="")])
        names(pest2A) <- paste("eXY", i, "-eXY", i-j+2, sep="")
        pest2 <- append(pest2, pest2A)

        if (Z != "NULL") {
          mcmcA <- (mcmc[, paste("eXZ", i, sep="")] - mcmc[, paste("eXZ", i-j+2, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("eXZ", i, "-eXZ", i-j+2, sep="")
          pest2A <- (pest2[paste("eXZ", i, sep="")] - pest2[paste("eXZ", i-j+2, sep="")])
          names(pest2A) <- paste("eXZ", i, "-eXZ", i-j+2, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("eYZ", i, sep="")] - mcmc[, paste("eYZ", i-j+2, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("eYZ", i, "-eYZ", i-j+2, sep="")
          pest2A <- (pest2[paste("eYZ", i, sep="")] - pest2[paste("eYZ", i-j+2, sep="")])
          names(pest2A) <- paste("eYZ", i, "-eYZ", i-j+2, sep="")
          pest2 <- append(pest2, pest2A)
        } # end (if Z != "NULL")

        if (W != "NULL") {
          mcmcA <- (mcmc[, paste("eXW", i, sep="")] - mcmc[, paste("eXW", i-j+2, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("eXW", i, "-eXW", i-j+2, sep="")
          pest2A <- (pest2[paste("eXW", i, sep="")] - pest2[paste("eXW", i-j+2, sep="")])
          names(pest2A) <- paste("eXW", i, "-eXW", i-j+2, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("eYW", i, sep="")] - mcmc[, paste("eYW", i-j+2, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("eYW", i, "-eYW", i-j+2, sep="")
          pest2A <- (pest2[paste("eYW", i, sep="")] - pest2[paste("eYW", i-j+2, sep="")])
          names(pest2A) <- paste("eYW", i, "-eYW", i-j+2, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("eZW", i, sep="")] - mcmc[, paste("eZW", i-j+2, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("eZW", i, "-eZW", i-j+2, sep="")
          pest2A <- (pest2[paste("eZW", i, sep="")] - pest2[paste("eZW", i-j+2, sep="")])
          names(pest2A) <- paste("eZW", i, "-eZW", i-j+2, sep="")
          pest2 <- append(pest2, pest2A)
        } # end (if W != "NULL")
      } # end (for i)
    } # end (for j)
  } # end (if eXY2)
  ## ----- end (Difference in covariance of latent variable residuals) ----- ##


  ## -- Differences in variance of latent variable residuals -- ##
  if (any(parEst[,4] == "eXX2")) {
    for (j in 3:no.waves) {
      for (i in j:no.waves) {
        mcmcA <- (mcmc[, paste("eXX", i, sep="")] - mcmc[, paste("eXX", i-j+2, sep="")])
        mcmc <- cbind(mcmc, mcmcA)
        colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("eXX", i, "-eXX", i-j+2, sep="")
        pest2A <- (pest2[paste("eXX", i, sep="")] - pest2[paste("eXX", i-j+2, sep="")])
        names(pest2A) <- paste("eXX", i, "-eXX", i-j+2, sep="")
        pest2 <- append(pest2, pest2A)

        mcmcA <- (mcmc[, paste("eYY", i, sep="")] - mcmc[, paste("eYY", i-j+2, sep="")])
        mcmc <- cbind(mcmc, mcmcA)
        colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("eYY", i, "-eYY", i-j+2, sep="")
        pest2A <- (pest2[paste("eYY", i, sep="")] - pest2[paste("eYY", i-j+2, sep="")])
        names(pest2A) <- paste("eYY", i, "-eYY", i-j+2, sep="")
        pest2 <- append(pest2, pest2A)

        if (Z != "NULL") {
          mcmcA <- (mcmc[, paste("eZZ", i, sep="")] - mcmc[, paste("eZZ", i-j+2, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("eZZ", i, "-eZZ", i-j+2, sep="")
          pest2A <- (pest2[paste("eZZ", i, sep="")] - pest2[paste("eZZ", i-j+2, sep="")])
          names(pest2A) <- paste("eZZ", i, "-eZZ", i-j+2, sep="")
          pest2 <- append(pest2, pest2A)
        } # end (if Z != "NULL")

        if (W != "NULL") {
          mcmcA <- (mcmc[, paste("eWW", i, sep="")] - mcmc[, paste("eWW", i-j+2, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("eWW", i, "-eWW", i-j+2, sep="")
          pest2A <- (pest2[paste("eWW", i, sep="")] - pest2[paste("eWW", i-j+2, sep="")])
          names(pest2A) <- paste("eWW", i, "-eWW", i-j+2, sep="")
          pest2 <- append(pest2, pest2A)
        } # end (if W != "NULL")
      } # end (for i)
    } # end (for j)
  } # end (if eXX2)
  ## ----- end (Difference in variance of latent variable residuals) ----- ##


  ## -- Differences in variance of impulses -- ##
  if (any(parEst[,4] == "iXX1")) {
    for (j in 2:no.waves) {
      for (i in j:no.waves) {
        mcmcA <- (mcmc[, paste("iXX", i, sep="")] - mcmc[, paste("iXX", i-j+1, sep="")])
        mcmc <- cbind(mcmc, mcmcA)
        colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("iXX", i, "-iXX", i-j+1, sep="")
        pest2A <- (pest2[paste("iXX", i, sep="")] - pest2[paste("iXX", i-j+1, sep="")])
        names(pest2A) <- paste("iXX", i, "-iXX", i-j+1, sep="")
        pest2 <- append(pest2, pest2A)

        mcmcA <- (mcmc[, paste("iYY", i, sep="")] - mcmc[, paste("iYY", i-j+1, sep="")])
        mcmc <- cbind(mcmc, mcmcA)
        colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("iYY", i, "-iYY", i-j+1, sep="")
        pest2A <- (pest2[paste("iYY", i, sep="")] - pest2[paste("iYY", i-j+1, sep="")])
        names(pest2A) <- paste("iYY", i, "-iYY", i-j+1, sep="")
        pest2 <- append(pest2, pest2A)

        if (Z != "NULL") {
          mcmcA <- (mcmc[, paste("iZZ", i, sep="")] - mcmc[, paste("iZZ", i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("iZZ", i, "-iZZ", i-j+1, sep="")
          pest2A <- (pest2[paste("iZZ", i, sep="")] - pest2[paste("iZZ", i-j+1, sep="")])
          names(pest2A) <- paste("iZZ", i, "-iZZ", i-j+1, sep="")
          pest2 <- append(pest2, pest2A)
        } # end (if Z != "NULL")

        if (W != "NULL") {
          mcmcA <- (mcmc[, paste("iWW", i, sep="")] - mcmc[, paste("iWW", i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("iWW", i, "-iWW", i-j+1, sep="")
          pest2A <- (pest2[paste("iWW", i, sep="")] - pest2[paste("iWW", i-j+1, sep="")])
          names(pest2A) <- paste("iWW", i, "-iWW", i-j+1, sep="")
          pest2 <- append(pest2, pest2A)
        } # end (if W != "NULL")
      } # end (for i)
    } # end (for j)
  } # end (if iXX1)
  ## ----- end (Difference in variance of latent variable residuals) ----- ##


  ## -- Differences in co-movements -- ##
  if (any(parEst[,4] == "iXY1")) {
    for (j in 2:no.waves) {
      for (i in j:no.waves) {
        mcmcA <- (mcmc[, paste("iXY", i, sep="")] - mcmc[, paste("iXY", i-j+1, sep="")])
        mcmc <- cbind(mcmc, mcmcA)
        colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("iXY", i, "-iXY", i-j+1, sep="")
        pest2A <- (pest2[paste("iXY", i, sep="")] - pest2[paste("iXY", i-j+1, sep="")])
        names(pest2A) <- paste("iXY", i, "-iXY", i-j+1, sep="")
        pest2 <- append(pest2, pest2A)

        if (Z != "NULL") {
          mcmcA <- (mcmc[, paste("iXZ", i, sep="")] - mcmc[, paste("iXZ", i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("iXZ", i, "-iXZ", i-j+1, sep="")
          pest2A <- (pest2[paste("iXZ", i, sep="")] - pest2[paste("iXZ", i-j+1, sep="")])
          names(pest2A) <- paste("iXZ", i, "-iXZ", i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("iYZ", i, sep="")] - mcmc[, paste("iYZ", i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("iYZ", i, "-iYZ", i-j+1, sep="")
          pest2A <- (pest2[paste("iYZ", i, sep="")] - pest2[paste("iYZ", i-j+1, sep="")])
          names(pest2A) <- paste("iYZ", i, "-iYZ", i-j+1, sep="")
          pest2 <- append(pest2, pest2A)
        } # end (if Z != "NULL")

        if (W != "NULL") {
          mcmcA <- (mcmc[, paste("iXW", i, sep="")] - mcmc[, paste("iXW", i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("iXW", i, "-iXW", i-j+1, sep="")
          pest2A <- (pest2[paste("iXW", i, sep="")] - pest2[paste("iXW", i-j+1, sep="")])
          names(pest2A) <- paste("iXW", i, "-iXW", i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("iYW", i, sep="")] - mcmc[, paste("iYW", i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("iYW", i, "-iYW", i-j+1, sep="")
          pest2A <- (pest2[paste("iYW", i, sep="")] - pest2[paste("iYW", i-j+1, sep="")])
          names(pest2A) <- paste("iYW", i, "-iYW", i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("iZW", i, sep="")] - mcmc[, paste("iZW", i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("iZW", i, "-iZW", i-j+1, sep="")
          pest2A <- (pest2[paste("iZW", i, sep="")] - pest2[paste("iZW", i-j+1, sep="")])
          names(pest2A) <- paste("iZW", i, "-iZW", i-j+1, sep="")
          pest2 <- append(pest2, pest2A)
        } # end (if W != "NULL")
      } # end (for i)
    } # end (for j)
  } # end (if iXY1)
  ## ----- end (Difference in covariance of latent variable residuals) ----- ##


  ## -- Differences in Indicator Residual Variance -- ##
  if (any(parEst[,4] == "eIXX1")) {
    for (j in 2:no.waves) {
      for (i in j:no.waves) {
        mcmcA <- (mcmc[, paste("eIXX", i, sep="")] - mcmc[, paste("eIXX", i-j+1, sep="")])
        mcmc <- cbind(mcmc, mcmcA)
        colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("eIXX", i, "-eIXX", i-j+1, sep="")
        pest2A <- (pest2[paste("eIXX", i, sep="")] - pest2[paste("eIXX", i-j+2, sep="")])
        names(pest2A) <- paste("eIXX", i, "-eIXX", i-j+1, sep="")
        pest2 <- append(pest2, pest2A)

        mcmcA <- (mcmc[, paste("eIYY", i, sep="")] - mcmc[, paste("eIYY", i-j+1, sep="")])
        mcmc <- cbind(mcmc, mcmcA)
        colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("eIYY", i, "-eIYY", i-j+1, sep="")
        pest2A <- (pest2[paste("eIYY", i, sep="")] - pest2[paste("eIYY", i-j+1, sep="")])
        names(pest2A) <- paste("eIYY", i, "-eIYY", i-j+1, sep="")
        pest2 <- append(pest2, pest2A)

        if (Z != "NULL") {
          mcmcA <- (mcmc[, paste("eIZZ", i, sep="")] - mcmc[, paste("eIZZ", i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("eIZZ", i, "-eIZZ", i-j+1, sep="")
          pest2A <- (pest2[paste("eIZZ", i, sep="")] - pest2[paste("eIZZ", i-j+1, sep="")])
          names(pest2A) <- paste("eIZZ", i, "-eIZZ", i-j+1, sep="")
          pest2 <- append(pest2, pest2A)
        } # end (if Z != "NULL")

        if (W != "NULL") {
          mcmcA <- (mcmc[, paste("eIWW", i, sep="")] - mcmc[, paste("eIWW", i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("eIWW", i, "-eIWW", i-j+1, sep="")
          pest2A <- (pest2[paste("eIWW", i, sep="")] - pest2[paste("eIWW", i-j+1, sep="")])
          names(pest2A) <- paste("eIWW", i, "-eIWW", i-j+1, sep="")
          pest2 <- append(pest2, pest2A)
        } # End (if W != "NULL")
      } # end (for i)
    } # end (for j)
  } else if (any(parEst[,4] == "eIXX2")) {
    for (j in 3:no.waves) {
      for (i in j:no.waves) {
        mcmcA <- (mcmc[, paste("eIXX", i, sep="")] - mcmc[, paste("eIXX", i-j+2, sep="")])
        mcmc <- cbind(mcmc, mcmcA)
        colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("eIXX", i, "-eIXX", i-j+2, sep="")
        pest2A <- (pest2[paste("eIXX", i, sep="")] - pest2[paste("eIXX", i-j+2, sep="")])
        names(pest2A) <- paste("eIXX", i, "-eIXX", i-j+2, sep="")
        pest2 <- append(pest2, pest2A)

        mcmcA <- (mcmc[, paste("eIYY", i, sep="")] - mcmc[, paste("eIYY", i-j+2, sep="")])
        mcmc <- cbind(mcmc, mcmcA)
        colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("eIYY", i, "-eIYY", i-j+2, sep="")
        pest2A <- (pest2[paste("eIYY", i, sep="")] - pest2[paste("eIYY", i-j+2, sep="")])
        names(pest2A) <- paste("eIYY", i, "-eIYY", i-j+2, sep="")
        pest2 <- append(pest2, pest2A)

        if (Z != "NULL") {
          mcmcA <- (mcmc[, paste("eIZZ", i, sep="")] - mcmc[, paste("eIZZ", i-j+2, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("eIZZ", i, "-eIZZ", i-j+2, sep="")
          pest2A <- (pest2[paste("eIZZ", i, sep="")] - pest2[paste("eIZZ", i-j+2, sep="")])
          names(pest2A) <- paste("eIZZ", i, "-eIZZ", i-j+2, sep="")
          pest2 <- append(pest2, pest2A)
        } # end (if Z != "NULL")

        if (W != "NULL") {
          mcmcA <- (mcmc[, paste("eIWW", i, sep="")] - mcmc[, paste("eIWW", i-j+2, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("eIWW", i, "-eIWW", i-j+2, sep="")
          pest2A <- (pest2[paste("eIWW", i, sep="")] - pest2[paste("eIWW", i-j+2, sep="")])
          names(pest2A) <- paste("eIWW", i, "-eIWW", i-j+2, sep="")
          pest2 <- append(pest2, pest2A)
        } # End (if W != "NULL")
      } # end (for i)
    } # end (for j)
  } # end (if eIXX)

  ## -------------------------------------------------------- ##


  ## -- Differences in indicator residual covariance -- ##
  if (any(parEst[,4] == "eIXY1")) {
    for (j in 2:no.waves) {
      for (i in j:no.waves) {
        mcmcA <- (mcmc[, paste0("eIXY", i, sep="")] - mcmc[, paste("eIXY", i-j+1, sep="")])
        mcmc <- cbind(mcmc, mcmcA)
        colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("eIXY", i, "-eIXY", i-j+1, sep="")
        pest2A <- (pest2[paste("eIXY", i, sep="")] - pest2[paste("eIXY", i-j+1, sep="")])
        names(pest2A) <- paste("eIXY", i, "-eIXY", i-j+1, sep="")
        pest2 <- append(pest2, pest2A)

        if (Z != "NULL") {
          mcmcA <- (mcmc[, paste("eIXZ", i, sep="")] - mcmc[, paste("eIXZ", i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("eIXZ", i, "-eIXZ", i-j+1, sep="")
          pest2A <- (pest2[paste("eIXZ", i, sep="")] - pest2[paste("eIXZ", i-j+1, sep="")])
          names(pest2A) <- paste("eIXZ", i, "-eIXZ", i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("eIYZ", i, sep="")] - mcmc[, paste("eIYZ", i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("eIYZ", i, "-eIYZ", i-j+1, sep="")
          pest2A <- (pest2[paste("eIYZ", i, sep="")] - pest2[paste("eIYZ", i-j+1, sep="")])
          names(pest2A) <- paste("eIYZ", i, "-eIYZ", i-j+1, sep="")
          pest2 <- append(pest2, pest2A)
        } # end (if Z != "NULL")

        if (W != "NULL") {
          mcmcA <- (mcmc[, paste("eIXW", i, sep="")] - mcmc[, paste("eIXW", i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("eIXW", i, "-eIXW", i-j+1, sep="")
          pest2A <- (pest2[paste("eIXW", i, sep="")] - pest2[paste("eIXW", i-j+1, sep="")])
          names(pest2A) <- paste("eIXW", i, "-eIXW", i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("eIYW", i, sep="")] - mcmc[, paste("eIYW", i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("eIYW", i, "-eIYW", i-j+1, sep="")
          pest2A <- (pest2[paste("eIYW", i, sep="")] - pest2[paste("eIYW", i-j+1, sep="")])
          names(pest2A) <- paste("eIYW", i, "-eIYW", i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("eIZW", i, sep="")] - mcmc[, paste("eIZW", i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("eIZW", i, "-eIZW", i-j+1, sep="")
          pest2A <- (pest2[paste("eIZW", i, sep="")] - pest2[paste("eIZW", i-j+1, sep="")])
          names(pest2A) <- paste("eIZW", i, "-eIZW", i-j+1, sep="")
          pest2 <- append(pest2, pest2A)
        } # end (if W != "NULL")
      } # end (for i)
    } # end (for j)
  } else if (any(parEst[,4] == "eIXY2")) {
    for (j in 3:no.waves) {
      for (i in j:no.waves) {
        mcmcA <- (mcmc[, paste0("eIXY", i, sep="")] - mcmc[, paste("eIXY", i-j+2, sep="")])
        mcmc <- cbind(mcmc, mcmcA)
        colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("eIXY", i, "-eIXY", i-j+2, sep="")
        pest2A <- (pest2[paste("eIXY", i, sep="")] - pest2[paste("eIXY", i-j+2, sep="")])
        names(pest2A) <- paste("eIXY", i, "-eIXY", i-j+2, sep="")
        pest2 <- append(pest2, pest2A)
        if (Z != "NULL") {
          mcmcA <- (mcmc[, paste("eIXZ", i, sep="")] - mcmc[, paste("eIXZ", i-j+2, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("eIXZ", i, "-eIXZ", i-j+2, sep="")
          pest2A <- (pest2[paste("eIXZ", i, sep="")] - pest2[paste("eIXZ", i-j+2, sep="")])
          names(pest2A) <- paste("eIXZ", i, "-eIXZ", i-j+2, sep="")
          pest2 <- append(pest2, pest2A)
          mcmcA <- (mcmc[, paste("eIYZ", i, sep="")] - mcmc[, paste("eIYZ", i-j+2, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("eIYZ", i, "-eIYZ", i-j+2, sep="")
          pest2A <- (pest2[paste("eIYZ", i, sep="")] - pest2[paste("eIYZ", i-j+2, sep="")])
          names(pest2A) <- paste("eIYZ", i, "-eIYZ", i-j+2, sep="")
          pest2 <- append(pest2, pest2A)
        } # end (if Z != "NULL")
        if (W != "NULL") {
          mcmcA <- (mcmc[, paste("eIXW", i, sep="")] - mcmc[, paste("eIXW", i-j+2, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("eIXW", i, "-eIXW", i-j+2, sep="")
          pest2A <- (pest2[paste("eIXW", i, sep="")] - pest2[paste("eIXW", i-j+2, sep="")])
          names(pest2A) <- paste("eIXW", i, "-eIXW", i-j+2, sep="")
          pest2 <- append(pest2, pest2A)
          mcmcA <- (mcmc[, paste("eIYW", i, sep="")] - mcmc[, paste("eIYW", i-j+2, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("eIYW", i, "-eIYW", i-j+2, sep="")
          pest2A <- (pest2[paste("eIYW", i, sep="")] - pest2[paste("eIYW", i-j+2, sep="")])
          names(pest2A) <- paste("eIYW", i, "-eIYW", i-j+2, sep="")
          pest2 <- append(pest2, pest2A)
          mcmcA <- (mcmc[, paste("eIZW", i, sep="")] - mcmc[, paste("eIZW", i-j+2, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("eIZW", i, "-eIZW", i-j+2, sep="")
          pest2A <- (pest2[paste("eIZW", i, sep="")] - pest2[paste("eIZW", i-j+2, sep="")])
          names(pest2A) <- paste("eIZW", i, "-eIZW", i-j+2, sep="")
          pest2 <- append(pest2, pest2A)
        } # end (if W != "NULL")
      } # end (for i)
    } # end (for j)
  } # end (if eIXY1)

  ## ----- end (Difference in indicator residual covariance) ----- ##


  ## -- Differences in Intercept of Latent Variables -- ##
  if (any(parEst[,4] == "MwX2")) {
    for (j in 3:no.waves) {
      for (i in j:no.waves) {
        mcmcA <- (mcmc[, paste("MwX", i, sep="")] - mcmc[, paste("MwX", i-j+2, sep="")])
        mcmc <- cbind(mcmc, mcmcA)
        colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("MwX", i, "-MwX", i-j+2, sep="")
        pest2A <- (pest2[paste("MwX", i, sep="")] - pest2[paste("MwX", i-j+2, sep="")])
        names(pest2A) <- paste("MwX", i, "-MwX", i-j+2, sep="")
        pest2 <- append(pest2, pest2A)

        mcmcA <- (mcmc[, paste("MwY", i, sep="")] - mcmc[, paste("MwY", i-j+2, sep="")])
        mcmc <- cbind(mcmc, mcmcA)
        colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("MwY", i, "-MwY", i-j+2, sep="")
        pest2A <- (pest2[paste("MwY", i, sep="")] - pest2[paste("MwY", i-j+2, sep="")])
        names(pest2A) <- paste("MwY", i, "-MwY", i-j+2, sep="")
        pest2 <- append(pest2, pest2A)

        if (Z != "NULL") {
          mcmcA <- (mcmc[, paste("MwZ", i, sep="")] - mcmc[, paste("MwZ", i-j+2, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("MwZ", i, "-MwZ", i-j+2, sep="")
          pest2A <- (pest2[paste("MwZ", i, sep="")] - pest2[paste("MwZ", i-j+2, sep="")])
          names(pest2A) <- paste("MwZ", i, "-MwZ", i-j+2, sep="")
          pest2 <- append(pest2, pest2A)
        } # end (if Z != "NULL")

        if (W != "NULL") {
          mcmcA <- (mcmc[, paste("MwW", i, sep="")] - mcmc[, paste("MwW", i-j+2, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("MwW", i, "-MwW", i-j+2, sep="")
          pest2A <- (pest2[paste("MwW", i, sep="")] - pest2[paste("MwW", i-j+2, sep="")])
          names(pest2A) <- paste("MwW", i, "-MwW", i-j+2, sep="")
          pest2 <- append(pest2, pest2A)
        } # End (if W != "NULL")
      } # end (for i)
    } # end (for j)
  } # end (if MwX2)
  ## ----- end (Differences in Intercepts) ----- ##


  ## -- Prepare pest2 file with BCCI and PCI p-values -- ##
  pest2 <- cbind(pest2, 0)
  colnames(pest2)[2] <- "pBCCI"
  pest2 <- cbind(pest2, 0)
  colnames(pest2)[3] <- "pPCI"

  no.pest2 <- nrow(pest2)

  for (i in 1: no.pest2) {
    estM <- pest2[i,1] # estimated parameter
    abM <- mcmc[, i] # simulated parameter
    zM = qnorm(sum(abM<estM)/b.no) # Bias-Corrected Factor

    # Calculate Bias-Corrected Probability
    if ((estM>0 & min(abM)>0) | (estM<0 & max(abM)<0)) {
      pest2[i, 2] = 0
    } else if (qnorm(sum(abM>0)/b.no)+2*zM<0) {
      pest2[i, 2] = 2*pnorm((qnorm(sum(abM>0)/b.no)+2*zM))
    } else {
      pest2[i, 2] = 2*pnorm(-1*(qnorm(sum(abM>0)/b.no)+2*zM))
    }

    # Calculate Percentile Probability
    if (quantile(abM,probs=0.5)>0) {
      pest2[i, 3] = 2*(sum(abM<0)/b.no)
    } else {
      pest2[i,3] = 2*(sum(abM>0)/b.no)
    }
  } # end (for i)
  ## ----- end (Prepare pest2 file) ----- ##


  ## ---- List and Delete - Path Coefficients pXX ---- ##
  if (any(parEst[,4] == "pXX21")) {
    cat(rep("\n",3), "## ===== Identification of invariant paths (Lag = 1 wave) ===== ##")
    no.path = no.waves - 1
    MIset <- no.waves - 3
    no.compare = (no.path - 1)*(no.path)/2
    if (Type1Adj == "BON") {
      p <- Type1/no.compare
    } else {
      p = Type1
    }
    cat("\n", sprintf("Critical p-value is: %.4f\n", p), "\n")

#    no.compare.M = (no.waves - 1)*(no.waves)/2
    MIset.M <- no.waves - 2

    LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X", b="X", LDlag=1)  ## List and Delete - Path XX ##
    LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Y", b="Y", LDlag=1)  ## List and Delete - Path YY ##
    if (Z != "NULL") {
      LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Z", b="Z", LDlag=1)  ## List and Delete - Path ZZ ##
    } # end (if Z)
    if (W != "NULL") {
      LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="W", b="W", LDlag=1)  ## List and Delete - Path WW ##
    } # end (if W)

    ## -- Lag = 2 -- ##
    if (lag == 2 & no.waves > 3) {
      cat(rep("\n",3), "## ===== Identification of invariant paths (Lag = 2 waves) ===== ##")
      no.compare = (no.path - 2)*(no.path - 1)/2
      MIset <- no.path - 3
      if (Type1Adj == "BON") {
        p <- Type1/no.compare
      } else {
        p = Type1
      }
      cat("\n", sprintf("Critical p-value is: %.4f\n", p), "\n")

      LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X", b="X", LDlag=2)  ## List and Delete - Path XX ##
      LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Y", b="Y", LDlag=2)  ## List and Delete - Path YY ##
      if (Z != "NULL") {
        LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Z", b="Z", LDlag=2)  ## List and Delete - Path ZZ ##
      } # end (if Z != "NULL")
      if (W != "NULL") {
        LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="W", b="W", LDlag=2)  ## List and Delete - Path WW ##
      } # end (if W != "NULL")
    } # end (Lag == 2)


    ## -- Lag = 3 -- ##
    if (lag == 3 & no.waves > 4) {
      cat(rep("\n",3), "## ===== Identification of invariant paths (Lag = 3 waves) ===== ##")
      no.compare = (no.path - 3)*(no.path - 2)/2
      MIset <- no.path - 4

      if (Type1Adj == "BON") {
        p <- Type1/no.compare
      } else {
        p = Type1
      }
      cat("\n", sprintf("Critical p-value is: %.4f\n", p), "\n")

      LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X", b="X", LDlag=3)  ## List and Delete - Path XX ##
      LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Y", b="Y", LDlag=3)  ## List and Delete - Path YY ##
      if (Z != "NULL") {
        LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Z", b="Z", LDlag=3)  ## List and Delete - Path ZZ ##
      } # end (if Z != "NULL")
      if (W != "NULL") {
        LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="W", b="W", LDlag=3)  ## List and Delete - Path WW ##
      } # end (if W != "NULL")
    } # end (lag == 3)
    ## ----- ##


    ## -- Lag = 4 -- ##
    if (lag == 4 & no.waves > 5) {
      cat(rep("\n",3), "## ===== Identification of invariant paths (Lag = 4 waves) ===== ##")
      no.compare = (no.path - 4)*(no.path - 3)/2
      MIset <- no.path - 5

      if (Type1Adj == "BON") {
        p <- Type1/no.compare
      } else {
        p = Type1
      } # end (if Type1Adj)
      cat("\n", sprintf("Critical p-value is: %.4f\n", p), "\n")

      LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X", b="X", LDlag=4)  ## List and Delete - Path XX ##
      LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Y", b="Y", LDlag=4)  ## List and Delete - Path YY ##
      if (Z != "NULL") {
        LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Z", b="Z", LDlag=4)  ## List and Delete - Path ZZ ##
      } # end (if Z != "NULL")
      if (W != "NULL") {
        LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="W", b="W", LDlag=4)  ## List and Delete - Path WW ##
      } # end (if W != "NULL")
    } # end (lag == 4)
  } # end (if pXX21)
  ## -------------------------------------------------------- ##


  ## ---- List and Delete - Path Coefficients pXY ---- ##
  if (any(parEst[,4] == "pXY21")) {
    cat(rep("\n",3), "## ===== Identification of invariant paths (Lag = 1 wave) ===== ##")
    no.path = no.waves - 1
    MIset <- no.waves - 3
    no.compare = (no.path - 1)*(no.path)/2
    if (Type1Adj == "BON") {
      p <- Type1/no.compare
    } else {
      p = Type1
    } # end (if Type1Adj)
    cat("\n", sprintf("Critical p-value is: %.4f\n", p), "\n")

#    no.compare.M = (no.waves - 1)*(no.waves)/2
    MIset.M <- no.waves - 2

    LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X", b="Y", LDlag=1)  ## List and Delete - Path XY ##
    LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Y", b="X", LDlag=1)  ## List and Delete - Path YX ##
    if (Z != "NULL") {
      LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X", b="Z", LDlag=1)  ## List and Delete - Path XZ ##
      LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Y", b="Z", LDlag=1)  ## List and Delete - Path YZ ##
      LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Z", b="X", LDlag=1)  ## List and Delete - Path ZX ##
      LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Z", b="Y", LDlag=1)  ## List and Delete - Path ZY ##
    } # end (if Z != "NULL")
    if (W != "NULL") {
      LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X", b="W", LDlag=1)  ## List and Delete - Path XW ##
      LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Y", b="W", LDlag=1)  ## List and Delete - Path YW ##
      LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Z", b="W", LDlag=1)  ## List and Delete - Path ZW ##
      LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="W", b="X", LDlag=1)  ## List and Delete - Path WX ##
      LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="W", b="Y", LDlag=1)  ## List and Delete - Path WY ##
      LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="W", b="Z", LDlag=1)  ## List and Delete - Path WZ ##
    } # end (if W != "NULL")


    ## -- Lag = 2 -- ##
    if (lag == 2 & no.waves > 3) {
      cat(rep("\n",3), "## ===== Identification of invariant paths (Lag = 2 waves) ===== ##")
      no.compare = (no.path - 2)*(no.path - 1)/2
      MIset <- no.path - 3
      if (Type1Adj == "BON") {
        p <- Type1/no.compare
      } else {
        p = Type1
      } # end (if Type1Adj)
      cat("\n", sprintf("Critical p-value is: %.4f\n", p), "\n")

      LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X", b="Y", LDlag=2)  ## List and Delete - Path XY ##
      LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Y", b="X", LDlag=2)  ## List and Delete - Path YX ##
      if (Z != "NULL") {
        LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X", b="Z", LDlag=2)  ## List and Delete - Path XZ ##
        LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Y", b="Z", LDlag=2)  ## List and Delete - Path YZ ##
        LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Z", b="X", LDlag=2)  ## List and Delete - Path ZX ##
        LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Z", b="Y", LDlag=2)  ## List and Delete - Path ZY ##
      } # end (if Z != "NULL")
      if (W != "NULL") {
        LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X", b="W", LDlag=2)  ## List and Delete - Path XW ##
        LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Y", b="W", LDlag=2)  ## List and Delete - Path YW ##
        LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Z", b="W", LDlag=2)  ## List and Delete - Path ZW ##
        LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="W", b="X", LDlag=2)  ## List and Delete - Path WX ##
        LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="W", b="Y", LDlag=2)  ## List and Delete - Path WY ##
        LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="W", b="Z", LDlag=2)  ## List and Delete - Path WZ ##
      } # end (if W != "NULL")
    } # end (Lag == 2)


    ## -- Lag = 3 -- ##
    if (lag == 3 & no.waves > 4) {
      cat(rep("\n",3), "## ===== Identification of invariant paths (Lag = 3 waves) ===== ##")
      no.compare = (no.path - 3)*(no.path - 2)/2
      MIset <- no.path - 4
      if (Type1Adj == "BON") {
        p <- Type1/no.compare
      } else {
        p = Type1
      } # end (if Type1Adj)
      cat("\n", sprintf("Critical p-value is: %.4f\n", p), "\n")

      LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X", b="Y", LDlag=3)  ## List and Delete - Path XY ##
      LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Y", b="X", LDlag=3)  ## List and Delete - Path YX ##
      if (Z != "NULL") {
        LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X", b="Z", LDlag=3)  ## List and Delete - Path XZ ##
        LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Y", b="Z", LDlag=3)  ## List and Delete - Path YZ ##
        LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Z", b="X", LDlag=3)  ## List and Delete - Path ZX ##
        LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Z", b="Y", LDlag=3)  ## List and Delete - Path ZY ##
      } # end (if Z != "NULL")
      if (W != "NULL") {
        LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X", b="W", LDlag=3)  ## List and Delete - Path XW ##
        LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Y", b="W", LDlag=3)  ## List and Delete - Path YW ##
        LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Z", b="W", LDlag=3)  ## List and Delete - Path ZW ##
        LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="W", b="X", LDlag=3)  ## List and Delete - Path WX ##
        LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="W", b="Y", LDlag=3)  ## List and Delete - Path WY ##
        LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="W", b="Z", LDlag=3)  ## List and Delete - Path WZ ##
      } # end (if W != "NULL")
    } # end (lag == 3)
    ## ----- ##


    ## -- Lag = 4 -- ##
    if (lag == 4 & no.waves > 5) {
      cat(rep("\n",3), "## ===== Identification of invariant paths (Lag = 4 waves) ===== ##")
      no.compare = (no.path - 4)*(no.path - 3)/2
      MIset <- no.path - 5
      if (Type1Adj == "BON") {
        p <- Type1/no.compare
      } else {
        p = Type1
      } # end (if Type1Adj)
      cat("\n", sprintf("Critical p-value is: %.4f\n", p), "\n")

      LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X", b="Y", LDlag=4)  ## List and Delete - Path XY ##
      LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Y", b="X", LDlag=4)  ## List and Delete - Path YX ##
      if (Z != "NULL") {
        LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X", b="Z", LDlag=4)  ## List and Delete - Path XZ ##
        LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Y", b="Z", LDlag=4)  ## List and Delete - Path YZ ##
        LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Z", b="X", LDlag=4)  ## List and Delete - Path ZX ##
        LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Z", b="Y", LDlag=4)  ## List and Delete - Path ZY ##
      } # end (if Z != "NULL")

      if (W != "NULL") {
        LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X", b="W", LDlag=4)  ## List and Delete - Path XW ##
        LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Y", b="W", LDlag=4)  ## List and Delete - Path YW ##
        LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Z", b="W", LDlag=4)  ## List and Delete - Path ZW ##
        LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="W", b="X", LDlag=4)  ## List and Delete - Path WX ##
        LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="W", b="Y", LDlag=4)  ## List and Delete - Path WY ##
        LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="W", b="Z", LDlag=4)  ## List and Delete - Path WZ ##
      } # end (if W != "NULL")
    } # end (lag == 4)
  } # end (if pXY21)
  ## -------------------------------------------------------- ##


  ## ---- List & Delete - Proportion Change dXX ---- ##
  if (any(parEst[,4] == "dXX21")) {
    cat(rep("\n",3), "## ===== Identification of invariant proportional change (Lag = 1 wave) ===== ##")
    no.path = no.waves - 1
    MIset <- no.waves - 3
    no.compare = (no.path - 1)*(no.path)/2
    if (Type1Adj == "BON") {
      p <- Type1/no.compare
    } else {
      p = Type1
    } # end (if Type1Adj)
    cat("\n", sprintf("Critical p-value is: %.4f\n", p), "\n")

#    no.compare.M = (no.waves - 1)*(no.waves)/2
    MIset.M <- no.waves - 2

    LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X", b="X", LDlag=1)  ## List & Delete - Proportion Change XX ##
    LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Y", b="Y", LDlag=1)  ## List & Delete - Proportion Change YY ##
    if (Z != "NULL") {
      LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Z", b="Z", LDlag=1)  ## List & Delete - Proportion Change ZZ ##
    } # end (if Z)
    if (W != "NULL") {
      LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="W", b="W", LDlag=1)  ## List & Delete - Proportion Change WW ##
    } # end (if W)

    ## -- Lag = 2 -- ##
    if (lag == 2 & no.waves > 3) {
      cat(rep("\n",3), "## ===== Identification of invariant  proportional change (Lag = 2 waves) ===== ##")
      no.compare = (no.path - 2)*(no.path - 1)/2
      MIset <- no.path - 3
      if (Type1Adj == "BON") {
        p <- Type1/no.compare
      } else {
        p = Type1
      } # end (if Type1Adj)
      cat("\n", sprintf("Critical p-value is: %.4f\n", p), "\n")

      LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X", b="X", LDlag=2)  ## List & Delete - Proportion Change XX ##
      LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Y", b="Y", LDlag=2)  ## List & Delete - Proportion Change YY ##
      if (Z != "NULL") {
        LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Z", b="Z", LDlag=2)  ## List & Delete - Proportion Change ZZ ##
      } # end (if Z != "NULL")
      if (W != "NULL") {
        LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="W", b="W", LDlag=2)  ## List & Delete - Proportion Change WW ##
      } # end (if W != "NULL")
    } # end (Lag == 2)


    ## -- Lag = 3 -- ##
    if (lag == 3 & no.waves > 4) {
      cat(rep("\n",3), "## ===== Identification of invariant  proportional change (Lag = 3 waves) ===== ##")
      no.compare = (no.path - 3)*(no.path - 2)/2
      MIset <- no.path - 4
      if (Type1Adj == "BON") {
        p <- Type1/no.compare
      } else {
        p = Type1
      } # end (if Type1Adj)
      cat("\n", sprintf("Critical p-value is: %.4f\n", p), "\n")

      LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X", b="X", LDlag=3)  ## List & Delete - Proportion Change XX ##
      LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Y", b="Y", LDlag=3)  ## List & Delete - Proportion Change YY ##
      if (Z != "NULL") {
        LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Z", b="Z", LDlag=3)  ## List & Delete - Proportion Change ZZ ##
      } # end (if Z != "NULL")
      if (W != "NULL") {
        LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="W", b="W", LDlag=3)  ## List & Delete - Proportion Change WW ##
      } # end (if W != "NULL")
    } # end (lag == 3)
    ## ----- ##


    ## -- Lag = 4 -- ##
    if (lag == 4 & no.waves > 5) {
      cat(rep("\n",3), "## ===== Identification of invariant  proportional change (Lag = 4 waves) ===== ##")
      no.compare = (no.path - 4)*(no.path - 3)/2
      MIset <- no.path - 5
      if (Type1Adj == "BON") {
        p <- Type1/no.compare
      } else {
        p = Type1
      } # end (if Type1Adj)
      cat("\n", sprintf("Critical p-value is: %.4f\n", p), "\n")

      LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X", b="X", LDlag=4)  ## List & Delete - Proportion Change XX ##
      LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Y", b="Y", LDlag=4)  ## List & Delete - Proportion Change YY ##
      if (Z != "NULL") {
        LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Z", b="Z", LDlag=4)  ## List & Delete - Proportion Change ZZ ##
      } # end (if Z != "NULL")
      if (W != "NULL") {
        LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="W", b="W", LDlag=4)  ## List & Delete - Proportion Change WW ##
      } # end (if W != "NULL")
    } # end (lag == 4)
  } # end (if dXX21)
  ## -------------------------------------------------------- ##


  ## ---- List & Delete - Proportion Change Coefficients dXY ---- ##
  if (any(parEst[,4] == "dXY21")) {
    cat(rep("\n",3), "## ===== Identification of invariant  proportional change (Lag = 1 wave) ===== ##")
    no.path = no.waves - 1
    MIset <- no.waves - 3
    no.compare = (no.path - 1)*(no.path)/2
    if (Type1Adj == "BON") {
      p <- Type1/no.compare
    } else {
      p = Type1
    } # end (if Type1Adj)
    cat("\n", sprintf("Critical p-value is: %.4f\n", p), "\n")

#    no.compare.M = (no.waves - 1)*(no.waves)/2
    MIset.M <- no.waves - 2

    LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X", b="Y", LDlag=1)  ## List & Delete - Proportion Change XY ##
    LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Y", b="X", LDlag=1)  ## List & Delete - Proportion Change YX ##
    if (Z != "NULL") {
      LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X", b="Z", LDlag=1)  ## List & Delete - Proportion Change XZ ##
      LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Y", b="Z", LDlag=1)  ## List & Delete - Proportion Change YZ ##
      LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Z", b="X", LDlag=1)  ## List & Delete - Proportion Change ZX ##
      LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Z", b="Y", LDlag=1)  ## List & Delete - Proportion Change ZY ##
    } # end (if Z != "NULL")
    if (W != "NULL") {
      LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X", b="W", LDlag=1)  ## List & Delete - Proportion Change XW ##
      LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Y", b="W", LDlag=1)  ## List & Delete - Proportion Change YW ##
      LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Z", b="W", LDlag=1)  ## List & Delete - Proportion Change ZW ##
      LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="W", b="X", LDlag=1)  ## List & Delete - Proportion Change WX ##
      LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="W", b="Y", LDlag=1)  ## List & Delete - Proportion Change WY ##
      LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="W", b="Z", LDlag=1)  ## List & Delete - Proportion Change WZ ##
    } # end (if W != "NULL")


    ## -- Lag = 2 -- ##
    if (lag == 2 & no.waves > 3) {
      cat(rep("\n",3), "## ===== Identification of invariant  proportional change (Lag = 2 waves) ===== ##")
      no.compare = (no.path - 2)*(no.path - 1)/2
      MIset <- no.path - 3
      if (Type1Adj == "BON") {
        p <- Type1/no.compare
      } else {
        p = Type1
      } # end (if Type1Adj)
      cat("\n", sprintf("Critical p-value is: %.4f\n", p), "\n")

      LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X", b="Y", LDlag=2)  ## List & Delete - Proportion Change XY ##
      LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Y", b="X", LDlag=2)  ## List & Delete - Proportion Change YX ##
      if (Z != "NULL") {
        LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X", b="Z", LDlag=2)  ## List & Delete - Proportion Change XZ ##
        LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Y", b="Z", LDlag=2)  ## List & Delete - Proportion Change YZ ##
        LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Z", b="X", LDlag=2)  ## List & Delete - Proportion Change ZX ##
        LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Z", b="Y", LDlag=2)  ## List & Delete - Proportion Change ZY ##
      } # end (if Z != "NULL")
      if (W != "NULL") {
        LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X", b="W", LDlag=2)  ## List & Delete - Proportion Change XW ##
        LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Y", b="W", LDlag=2)  ## List & Delete - Proportion Change YW ##
        LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Z", b="W", LDlag=2)  ## List & Delete - Proportion Change ZW ##
        LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="W", b="X", LDlag=2)  ## List & Delete - Proportion Change WX ##
        LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="W", b="Y", LDlag=2)  ## List & Delete - Proportion Change WY ##
        LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="W", b="Z", LDlag=2)  ## List & Delete - Proportion Change WZ ##
      } # end (if W != "NULL")
    } # end (Lag == 2)


    ## -- Lag = 3 -- ##
    if (lag == 3 & no.waves > 4) {
      cat(rep("\n",3), "## ===== Identification of invariant  proportional change (Lag = 3 waves) ===== ##")
      no.compare = (no.path - 3)*(no.path - 2)/2
      MIset <- no.path - 4
      if (Type1Adj == "BON") {
        p <- Type1/no.compare
      } else {
        p = Type1
      } # end (if Type1Adj)
      cat("\n", sprintf("Critical p-value is: %.4f\n", p), "\n")

      LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X", b="Y", LDlag=3)  ## List & Delete - Proportion Change XY ##
      LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Y", b="X", LDlag=3)  ## List & Delete - Proportion Change YX ##
      if (Z != "NULL") {
        LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X", b="Z", LDlag=3)  ## List & Delete - Proportion Change XZ ##
        LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Y", b="Z", LDlag=3)  ## List & Delete - Proportion Change YZ ##
        LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Z", b="X", LDlag=3)  ## List & Delete - Proportion Change ZX ##
        LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Z", b="Y", LDlag=3)  ## List & Delete - Proportion Change ZY ##
      } # end (if Z != "NULL")
      if (W != "NULL") {
        LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X", b="W", LDlag=3)  ## List & Delete - Proportion Change XW ##
        LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Y", b="W", LDlag=3)  ## List & Delete - Proportion Change YW ##
        LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Z", b="W", LDlag=3)  ## List & Delete - Proportion Change ZW ##
        LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="W", b="X", LDlag=3)  ## List & Delete - Proportion Change WX ##
        LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="W", b="Y", LDlag=3)  ## List & Delete - Proportion Change WY ##
        LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="W", b="Z", LDlag=3)  ## List & Delete - Proportion Change WZ ##
      } # end (if W != "NULL")
    } # end (lag == 3)
    ## ----- ##


    ## -- Lag = 4 -- ##
    if (lag == 4 & no.waves > 5) {
      cat(rep("\n",3), "## ===== Identification of invariant  proportional change (Lag = 4 waves) ===== ##")
      no.compare = (no.path - 4)*(no.path - 3)/2
      MIset <- no.path - 5
      if (Type1Adj == "BON") {
        p <- Type1/no.compare
      } else {
        p = Type1
      } # end (if Type1Adj)
      cat("\n", sprintf("Critical p-value is: %.4f\n", p), "\n")

      LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X", b="Y", LDlag=4)  ## List & Delete - Proportion Change XY ##
      LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Y", b="X", LDlag=4)  ## List & Delete - Proportion Change YX ##
      if (Z != "NULL") {
        LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X", b="Z", LDlag=4)  ## List & Delete - Proportion Change XZ ##
        LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Y", b="Z", LDlag=4)  ## List & Delete - Proportion Change YZ ##
        LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Z", b="X", LDlag=4)  ## List & Delete - Proportion Change ZX ##
        LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Z", b="Y", LDlag=4)  ## List & Delete - Proportion Change ZY ##
      } # end (if Z != "NULL")

      if (W != "NULL") {
        LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X", b="W", LDlag=4)  ## List & Delete - Proportion Change XW ##
        LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Y", b="W", LDlag=4)  ## List & Delete - Proportion Change YW ##
        LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Z", b="W", LDlag=4)  ## List & Delete - Proportion Change ZW ##
        LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="W", b="X", LDlag=4)  ## List & Delete - Proportion Change WX ##
        LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="W", b="Y", LDlag=4)  ## List & Delete - Proportion Change WY ##
        LandD_PC(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="W", b="Z", LDlag=4)  ## List & Delete - Proportion Change WZ ##
      } # end (if W != "NULL")
    } # end (lag == 4)
  } # end (if dXY21)
  ## -------------------------------------------------------- ##


  ## ---- List and Delete - Residual variance eXX ---- ##
  if (any(parEst[,4] == "eXX2")) {
    cat(rep("\n",5), "## ===== Identification of invariant residual variance ===== ##")

    # -- Reset MISet and no.compare for residual variance --#
    no.path = no.waves - 1
    MIset <- no.waves - 3
    no.compare = (no.path - 1)*(no.path)/2
    if (Type1Adj == "BON") {
      p <- Type1/no.compare
    } else {
      p = Type1
    } # end (if Type1Adj)
    cat("\n", sprintf("Critical p-value is: %.4f\n", p), "\n")

    LandD_eXX(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X")  ## List and Delete - residual variance eXX ##
    LandD_eXX(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Y")  ## List and Delete - residual variance eYY ##

    if (Z != "NULL") {
      LandD_eXX(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Z")  ## List and Delete - residual variance eZZ ##
    } # end (if Z != "NULL")

    if (W != "NULL") {
      LandD_eXX(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="W")  ## List and Delete - residual variance eWW ##
    } # end (if W != "NULL")
  } # end (if eXX2)

  ## -------------------------------------------------------- ##


  ## ---- List and Delete - Residual covariance eXY ---- ##
  if (any(parEst[,4] == "eXY2")) {
    cat(rep("\n",5), "## ===== Identification of invariant residual covariance ===== ##")

    # -- Reset MISet and no.compare for residual covariance -- #
    no.path = no.waves - 1
    MIset <- no.waves - 3
    no.compare = (no.path - 1)*(no.path)/2
    if (Type1Adj == "BON") {
      p <- Type1/no.compare
    } else {
      p = Type1
    } # end (if Type1Adj)
    cat("\n", sprintf("Critical p-value is: %.4f\n", p), "\n")

    LandD_eXY(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X", b="Y")  ## List and Delete - eXY ##

    if (Z != "NULL") {
      LandD_eXY(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X", b="Z")  ## List and Delete - eXZ ##
      LandD_eXY(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Y", b="Z")  ## List and Delete - eYZ ##
    } # end (if Z != "NULL")

    if (W != "NULL") {
      LandD_eXY(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X", b="W")  ## List and Delete - eXW ##
      LandD_eXY(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Y", b="W")  ## List and Delete - eYW ##
      LandD_eXY(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Z", b="W")  ## List and Delete - eZW ##
    } # end (if W != "NULL")
  } # end (if eXY2)
  ## -------------------------------------------------------- ##


  ## ---- List and Delete - Impulse variance iXX ---- ##
  if (any(parEst[,4] == "iXX2")) {
    cat(rep("\n",5), "## ===== Identification of invariant impulse variance ===== ##")

    # -- Reset MISet and no.compare for impulse variance --#
    no.path = no.waves - 1
    MIset <- no.waves - 2
    no.compare = (no.waves - 1)*(no.waves)/2
    if (Type1Adj == "BON") {
      p <- Type1/no.compare
    } else {
      p = Type1
    } # end (if Type1Adj)
    cat("\n", sprintf("Critical p-value is: %.4f\n", p), "\n")

    LandD_iXX(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X")  ## List and Delete - impulse variance iX ##
    LandD_iXX(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Y")  ## List and Delete - impulse variance iY ##

    if (Z != "NULL") {
      LandD_iXX(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Z")  ## List and Delete - impulse variance iZ ##
    } # end (if Z != "NULL")

    if (W != "NULL") {
      LandD_iXX(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="W")  ## List and Delete - impulse variance iW ##
    } # end (if W != "NULL")
  } # end (if iXX2)

  ## -------------------------------------------------------- ##


  ## ---- List and Delete - Co-movements iXY ---- ##
  if (any(parEst[,4] == "iXY2")) {
    cat(rep("\n",5), "## ===== Identification of invariant co-movements ===== ##")

    # -- Reset MISet and no.compare for co-movements -- #
    no.path = no.waves - 1
    MIset <- no.waves - 2
    no.compare = (no.waves - 1)*(no.waves)/2
    if (Type1Adj == "BON") {
      p <- Type1/no.compare
    } else {
      p = Type1
    } # end (if Type1Adj)
    cat("\n", sprintf("Critical p-value is: %.4f\n", p), "\n")

    LandD_iXY(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X", b="Y")  ## List and Delete - iXY ##

    if (Z != "NULL") {
      LandD_iXY(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X", b="Z")  ## List and Delete - iXZ ##
      LandD_iXY(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Y", b="Z")  ## List and Delete - iYZ ##
    } # end (if Z != "NULL")

    if (W != "NULL") {
      LandD_iXY(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X", b="W")  ## List and Delete - iXW ##
      LandD_iXY(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Y", b="W")  ## List and Delete - iYW ##
      LandD_iXY(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Z", b="W")  ## List and Delete - iZW ##
    } # end (if W != "NULL")
  } # end (if iXY2)
  ## -------------------------------------------------------- ##


  ## ---- List and Delete - Indicator residual variance eIXX ---- ##
  if (any(parEst[,4] == "eIXX2")) {
    cat(rep("\n",5), "## ===== Identification of invariant indicator residual variance ===== ##")

    # -- Reset MISet and no.compare for residual variance --#
    no.path = no.waves - 1
    MIset <- no.waves - 3
    no.compare = (no.path - 1)*(no.path)/2
    if (any(parEst[,4] == "eIXX1")) {
      no.path = no.waves - 1
      MIset <- no.waves - 2
      no.compare = (no.waves - 1)*(no.waves)/2
    } # end (if eIXX1)
    if (Type1Adj == "BON") {
      p <- Type1/no.compare
    } else {
      p = Type1
    } # end (if Type1Adj)
    cat("\n", sprintf("Critical p-value is: %.4f\n", p), "\n")

    LandD_eIXX(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X")  ## List and Delete - indicator residual variance eIXX ##
    LandD_eIXX(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Y")  ## List and Delete - indicator residual variance eIYY ##

    if (Z != "NULL") {
      LandD_eIXX(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Z")  ## List and Delete - indicator residual variance eIZZ ##
    } # end (if Z != "NULL")

    if (W != "NULL") {
      LandD_eIXX(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="W")  ## List and Delete - indicator residual variance eIWW ##
    } # end (if W != "NULL")
  } # end (if eIXX2)

  ## -------------------------------------------------------- ##


  ## ---- List and Delete - Indicator residual covariance eIXY ---- ##
  if (any(parEst[,4] == "eIXY2")) {
    cat(rep("\n",5), "## ===== Identification of invariant residual covariance of indicators ===== ##")

    # -- Reset MISet and no.compare for residual variance --#
    no.path = no.waves - 1
    MIset <- no.waves - 3
    no.compare = (no.path - 1)*(no.path)/2
    if (any(parEst[,4] == "eIXY1")) {
      no.path = no.waves - 1
      MIset <- no.waves - 2
      no.compare = (no.waves - 1)*(no.waves)/2
    } # end (if eIXY1)
    if (Type1Adj == "BON") {
      p <- Type1/no.compare
    } else {
      p = Type1
    } # end (if Type1Adj)
    cat("\n", sprintf("Critical p-value is: %.4f\n", p), "\n")


    LandD_eIXY(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X", b="Y")  ## List and Delete - eIXY ##

    if (Z != "NULL") {
      LandD_eIXY(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X", b="Z")  ## List and Delete - eIXZ ##
      LandD_eIXY(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Y", b="Z")  ## List and Delete - eIYZ ##
    } # end (if Z != "NULL")

    if (W != "NULL") {
      LandD_eIXY(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X", b="W")  ## List and Delete - eIXW ##
      LandD_eIXY(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Y", b="W")  ## List and Delete - eIYW ##
      LandD_eIXY(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Z", b="W")  ## List and Delete - eIZW ##
    } # end (if W != "NULL")
  } # end (if eIXY2)

  ## -------------------------------------------------------- ##


  ## ---- List and Delete - Intercepts ---- ##
  if (any(parEst[,4] == "MwX2")) {
    cat(rep("\n",5), "## ===== Identification of invariant intercepts ===== ##")

    # -- Reset MISet and no.compare for intercepts --#
    no.path = no.waves - 1
    MIset <- no.waves - 3
    no.compare = (no.path - 1)*(no.path)/2
    if (Type1Adj == "BON") {
      p <- Type1/no.compare
    } else {
      p = Type1
    } # end (if Type1Adj)
    cat("\n", sprintf("Critical p-value is: %.4f\n", p), "\n")


    LandD_MEAN(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X")  ## List and Delete - Intercepts of X ##
    LandD_MEAN(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Y")  ## List and Delete - Intercepts of Y ##

    if (Z != "NULL") {
      LandD_MEAN(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Z")  ## List and Delete - Intercepts of Z ##
    } # end (if Z != "NULL")

    if (W != "NULL") {
      LandD_MEAN(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="W")  ## List and Delete - Intercepts of W ##
    } # end (if W != "NULL")
  } # end (if MwX1)

  ## -------------------------------------------------------- ##

} # end (Function Invariance)
## -------------------------------------------------------- ##





## ----- Sub-Function List and Delete for Path ----- ##
LandD_Path <- function(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X", b="Y", LDlag=1) {

  SumEst <- 0
  cat(rep("\n",3), paste("# -- Path ", a, b, " coefficients -- #", sep=""))
  aa <- get(a)
  bb <- get(b)
  for (i in 1: (no.path-LDlag+1)) {
    Clhs <- paste("w", bb, i+LDlag, sep="")
    Crhs <- paste("w", aa, i, sep="")
    TparEst <- parEst[parEst["lhs"] == Clhs & parEst["rhs"] == Crhs & parEst["op"] == "~",]
    p.TparEst <- paste("  p", a, b, i+LDlag, i, ":  ", Clhs, " ~ ", Crhs, " = ", format(round(TparEst["est"], digits=4), nsmall=4, scientific=FALSE),
                       ", p-value = ", format(round(TparEst["pvalue"], digits=4), nsmall=4, scientific=FALSE), sep="")
    cat("\n", p.TparEst)
    SumEst <- SumEst + TparEst["est"]
  }
  MeanEst <- SumEst/(no.path-LDlag+1)
  p.MeanEst <- paste("  Average across paths = ", format(round(MeanEst, digits=4), nsmall=4, scientific=FALSE), sep="")
  cat("\n", p.MeanEst)

  # -- Save pairs of non-invariant paths -- #
  NI.path <- t(matrix(0,ncol=no.compare, nrow=2))
  k = 1
  for (j in 3:(no.waves-LDlag+1)) {
    for (i in j:(no.waves-LDlag+1)) {
      Clhs <- paste0("p", a, b, (i+LDlag-1), i-1, "-p", a, b, i-j+LDlag+1, i-j+1, sep="")
      if (pest2[Clhs,3] < p) {
          NI.path[k, 1] <- i-1
          NI.path[k, 2] <- i-j+1
          k <- k + 1
      } # end (if pest2)
    } # end (for i)
  } # end (for j)
  Count.NI.path = k - 1

  # -- Select sets of invariant paths and print -- #
  cat(rep("\n",2), paste("# -- Sets of invariant Path ", a, b, " coefficients -- #", sep=""))

  for (k in 0:MIset) {
    Noset <- no.path - k - LDlag + 1
    NIset <- factorial(no.path-LDlag+1)/(factorial(no.path-k-LDlag+1)*factorial(k))
    mm <- t(combn((no.path-LDlag+1), Noset))
    for (j in 1:Count.NI.path) {
      p.j <- NI.path[j,1]
      q.j <- NI.path[j,2]
      for (i in 1:NIset) {
        count4 <- length(which(mm[i,] == p.j | mm[i,] == q.j))
        if (count4 > 1) {
          mm[i,] <- 0
        } # end (if count4)
      } # end (for i)
    } # end (for j)
    for (ii in 1:NIset) {
      count4 <- sum(mm[ii,])
      if (count4 > 0) {
        mm.p <- mm
        for (i in 1:(no.path-LDlag+1)) {
          mm.p[mm.p == i] <- paste0("p", a, b, i+LDlag, i, sep="")
        } # end (for i)
        cat("\n", "    ", mm.p[ii,])
      } # end (if count4)
    } # end (for ii)
  } # end (for k)
}  # end (Function LandD_Path)
## --------------------------------------------------------- ##




## ----- Sub-Function List and Delete for Proportional Change ----- ##
LandD_PC <- function(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X", b="Y", LDlag=1) {

  SumEst <- 0
  cat(rep("\n",3), paste("# -- Proportional change ", a, b, " coefficients -- #", sep=""))
  aa <- get(a)
  bb <- get(b)
  for (i in 1: (no.path-LDlag+1)) {
    Clhs <- paste("cs", bb, i+LDlag, sep="")
    Crhs <- paste("w", aa, i, sep="")
    TparEst <- parEst[parEst["lhs"] == Clhs & parEst["rhs"] == Crhs & parEst["op"] == "~",]
    p.TparEst <- paste("  d", a, b, i+LDlag, i, ":  ", Clhs, " ~ ", Crhs, " = ", format(round(TparEst["est"], digits=4), nsmall=4, scientific=FALSE),
                       ", p-value = ", format(round(TparEst["pvalue"], digits=4), nsmall=4, scientific=FALSE), sep="")
    cat("\n", p.TparEst)
    SumEst <- SumEst + TparEst["est"]
  }
  MeanEst <- SumEst/(no.path-LDlag+1)
  p.MeanEst <- paste("  Average across proportional change = ", format(round(MeanEst, digits=4), nsmall=4, scientific=FALSE), sep="")
  cat("\n", p.MeanEst)

  # -- Save pairs of non-invariant proportional change -- #
  NI.path <- t(matrix(0,ncol=no.compare, nrow=2))
  k = 1
  for (j in 3:(no.waves-LDlag+1)) {
    for (i in j:(no.waves-LDlag+1)) {
      Clhs <- paste0("d", a, b, (i+LDlag-1), i-1, "-d", a, b, i-j+LDlag+1, i-j+1, sep="")
      if (pest2[Clhs,3] < p) {
          NI.path[k, 1] <- i-1
          NI.path[k, 2] <- i-j+1
          k <- k + 1
      } # end (if pest2)
    } # end (for i)
  } # end (for j)
  Count.NI.path = k - 1

  # -- Select sets of invariant proportional change and print -- #
  cat(rep("\n",2), paste("# -- Sets of invariant Proportional Change ", a, b, " coefficients -- #", sep=""))

  for (k in 0:MIset) {
    Noset <- no.path - k - LDlag + 1
    NIset <- factorial(no.path-LDlag+1)/(factorial(no.path-k-LDlag+1)*factorial(k))
    mm <- t(combn((no.path-LDlag+1), Noset))
    for (j in 1:Count.NI.path) {
      p.j <- NI.path[j,1]
      q.j <- NI.path[j,2]
      for (i in 1:NIset) {
        count4 <- length(which(mm[i,] == p.j | mm[i,] == q.j))
        if (count4 > 1) {
          mm[i,] <- 0
        } # end (if count4)
      } # end (for i)
    } # end (for j)
    for (ii in 1:NIset) {
      count4 <- sum(mm[ii,])
      if (count4 > 0) {
        mm.p <- mm
        for (i in 1:(no.path-LDlag+1)) {
          mm.p[mm.p == i] <- paste0("d", a, b, i+LDlag, i, sep="")
        } # end (for i)
        cat("\n", "    ", mm.p[ii,])
      } # end (if count4)
    } # end (for ii)
  } # end (for k)
}  # end (Function LandD_PC)
## --------------------------------------------------------- ##




## ----- Sub-Function List and Delete for Residual Variance of latent variables ----- ##
LandD_eXX <- function(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X") {

  SumEst <- 0
  cat(rep("\n",3), paste("# -- Residual variance of e", a, a, " coefficients -- #", sep=""))
  aa <- get(a)
  for (i in 1: no.path) {
    Clhs <- paste("w", aa, i+1, sep="")
    Crhs <- paste("w", aa, i+1, sep="")
    TparEst <- parEst[parEst["lhs"] == Clhs & parEst["rhs"] == Crhs & parEst["op"] == "~~",]
    p.TparEst <- paste("  e", a, a, i+1, ":  ", Clhs, " ~~ ", Crhs, " = ", format(round(TparEst["est"], digits=4), nsmall=4, scientific=FALSE),
                       ", p-value = ", format(round(TparEst["pvalue"], digits=4), nsmall=4, scientific=FALSE), sep="")
    cat("\n", p.TparEst)
    SumEst <- SumEst + TparEst["est"]
  } # end (for i)
  MeanEst <- SumEst/(no.path)
  p.MeanEst <- paste("  Average across waves from T2 = ", format(round(MeanEst, digits=4), nsmall=4, scientific=FALSE), sep="")
  cat("\n", p.MeanEst)

  # -- Save pairs of non-invariant residual variances -- #
  NI.path <- t(matrix(0,ncol=no.compare, nrow=2))
  k = 1
  for (j in 3:no.waves) {
    for (i in j:no.waves) {
      Clhs <- paste("e", a, a, i, "-e", a, a, i-j+2, sep="")
      if (pest2[Clhs,3] < p) {
        NI.path[k, 1] <- i
        NI.path[k, 2] <- i-j+2
        k <- k + 1
      } # end (if pest2)
    } # end (for i)
  } # end (for j)
  Count.NI.path = k - 1

  ## -- Select sets of invariant residual variances and print -- ##
  cat(rep("\n",2), paste("# -- Sets of invariant residual variance e", a, a, " coefficients -- #", sep=""))

  for (k in 0:MIset) {
    Noset <- no.path - k
    NIset <- factorial(no.path)/(factorial(no.path-k)*factorial(k))
    mm <- t(combn((no.path), Noset))
    for (j in 1:Count.NI.path) {
      p.j <- NI.path[j,1]
      q.j <- NI.path[j,2]
      for (i in 1:NIset) {
        count4 <- length(which(mm[i,] == p.j | mm[i,] == q.j))
        if (count4 > 1) {
          mm[i,] <- 0
        } # end (if count4)
      } # end (for i)
    } # end (for j)
    for (ii in 1:NIset) {
      count4 <- sum(mm[ii,])
      if (count4 > 0) {
        mm.p <- mm
        for (i in 1:no.path) {
          mm.p[mm.p == i] <- paste("e", a, a, i+1, sep="")
        } # end (for i)
        cat("\n", "    ", mm.p[ii,])
      } # end (if count4)
    } # end (for ii)
  } # end (for k)
}  # end (Function LandD_eXX)
## ---------------------------------------------------------------- ##





## ----- Sub-Function List and Delete for Residual Covariance of latent variables ----- ##
LandD_eXY <- function(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X", b="Y") {

  SumEst <- 0
  cat(rep("\n",3), paste("# -- Residual covariance of e", a, b, " coefficients -- #", sep=""))
  aa <- get(a)
  bb <- get(b)
  for (i in 1: no.path) {
    Clhs <- paste("w", aa, i+1, sep="")
    Crhs <- paste("w", bb, i+1, sep="")
    TparEst <- parEst[parEst["lhs"] == Clhs & parEst["rhs"] == Crhs & parEst["op"] == "~~",]
    p.TparEst <- paste("  e", a, b, i+1, ":  ", Clhs, " ~~ ", Crhs, " = ", format(round(TparEst["est"], digits=4), nsmall=4, scientific=FALSE),
                       ", p-value = ", format(round(TparEst["pvalue"], digits=4), nsmall=4, scientific=FALSE), sep="")
    cat("\n", p.TparEst)
    SumEst <- SumEst + TparEst["est"]
  } # end (for i)
  MeanEst <- SumEst/(no.path)
  p.MeanEst <- paste("  Average across waves from T2 = ", format(round(MeanEst, digits=4), nsmall=4, scientific=FALSE), sep="")
  cat("\n", p.MeanEst)

  # -- Save pairs of non-invariant residual covariances -- #
  NI.path <- t(matrix(0,ncol=no.compare, nrow=2))
  k = 1
  for (j in 3:no.waves) {
    for (i in j:no.waves) {
      Clhs <- paste("e", a, b, i, "-e", a, b, i-j+2, sep="")
      if (pest2[Clhs,3] < p) {
        NI.path[k, 1] <- i
        NI.path[k, 2] <- i-j+2
        k <- k + 1
      } # end (if pest2)
    } # end (for i)
  } # end (for j)
  Count.NI.path = k - 1

  ## -- Select sets of invariant residual covariances and print -- ##
  cat(rep("\n",2), paste("# -- Sets of invariant residual covariance e", a, b, " coefficients -- #", sep=""))

  for (k in 0:MIset) {
    Noset <- no.path - k
    NIset <- factorial(no.path)/(factorial(no.path-k)*factorial(k))
    mm <- t(combn((no.path), Noset))
    for (j in 1:Count.NI.path) {
      p.j <- NI.path[j,1]
      q.j <- NI.path[j,2]
      for (i in 1:NIset) {
        count4 <- length(which(mm[i,] == p.j | mm[i,] == q.j))
        if (count4 > 1) {
          mm[i,] <- 0
        } # end (if count4)
      } # end (for i)
    } # end (for j)
    for (ii in 1:NIset) {
      count4 <- sum(mm[ii,])
      if (count4 > 0) {
        mm.p <- mm
        for (i in 1:no.path) {
          mm.p[mm.p == i] <- paste("e", a, b, i+1, sep="")
        } # end (for i)
        cat("\n", "    ", mm.p[ii,])
      } # end (if count4)
    } # end (for ii)
  } # end (for k)
}  # end (Function LandD_eXY)
## ---------------------------------------------------------------- ##





## ----- Sub-Function List and Delete for Impulse Variance ----- ##
LandD_iXX <- function(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X") {

  SumEst <- 0
  cat(rep("\n",3), paste("# -- Impulse variance of i", a, a, " coefficients -- #", sep=""))
  aa <- get(a)
  for (i in 1: no.waves) {
    Clhs <- paste("d", aa, i, sep="")
    Crhs <- paste("d", aa, i, sep="")
    TparEst <- parEst[parEst["lhs"] == Clhs & parEst["rhs"] == Crhs & parEst["op"] == "~~",]
    p.TparEst <- paste("  i", a, a, i, ":  ", Clhs, " ~~ ", Crhs, " = ", format(round(TparEst["est"], digits=4), nsmall=4, scientific=FALSE),
                       ", p-value = ", format(round(TparEst["pvalue"], digits=4), nsmall=4, scientific=FALSE), sep="")
    cat("\n", p.TparEst)
    SumEst <- SumEst + TparEst["est"]
  } # end (for i)
  MeanEst <- SumEst/(no.waves)
  p.MeanEst <- paste("  Average across waves from T1 = ", format(round(MeanEst, digits=4), nsmall=4, scientific=FALSE), sep="")
  cat("\n", p.MeanEst)

  # -- Save pairs of non-invariant impulse variances -- #
  NI.path <- t(matrix(0,ncol=no.compare, nrow=2))
  k = 1
  for (j in 2:no.waves) {
    for (i in j:no.waves) {
      Clhs <- paste("i", a, a, i, "-i", a, a, i-j+1, sep="")
      if (pest2[Clhs,3] < p) {
        NI.path[k, 1] <- i
        NI.path[k, 2] <- i-j+1
        k <- k + 1
      } # end (if pest2)
    } # end (for i)
  } # end (for j)
  Count.NI.path = k - 1

  ## -- Select sets of invariant impulse variances and print -- ##
  cat(rep("\n",2), paste("# -- Sets of invariant impulse variance i", a, a, " coefficients -- #", sep=""))

  for (k in 0:MIset) {
    Noset <- no.waves - k
    NIset <- factorial(no.waves)/(factorial(no.waves-k)*factorial(k))
    mm <- t(combn((no.waves), Noset))
    for (j in 1:Count.NI.path) {
      p.j <- NI.path[j,1]
      q.j <- NI.path[j,2]
      for (i in 1:NIset) {
        count4 <- length(which(mm[i,] == p.j | mm[i,] == q.j))
        if (count4 > 1) {
          mm[i,] <- 0
        } # end (if count4)
      } # end (for i)
    } # end (for j)
    for (ii in 1:NIset) {
      count4 <- sum(mm[ii,])
      if (count4 > 0) {
        mm.p <- mm
        for (i in 1:no.waves) {
          mm.p[mm.p == i] <- paste("i", a, a, i, sep="")
        } # end (for i)
        cat("\n", "    ", mm.p[ii,])
      } # end (if count4)
    } # end (for ii)
  } # end (for k)
}  # end (Function LandD_iXX)
## ---------------------------------------------------------------- ##





## ----- Sub-Function List and Delete for Co-Movements ----- ##
LandD_iXY <- function(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X", b="Y") {

  SumEst <- 0
  cat(rep("\n",3), paste("# -- Co-Movements of i", a, b, " coefficients -- #", sep=""))
  aa <- get(a)
  bb <- get(b)
  for (i in 1: no.waves) {
    Clhs <- paste("d", aa, i, sep="")
    Crhs <- paste("d", bb, i, sep="")
    TparEst <- parEst[parEst["lhs"] == Clhs & parEst["rhs"] == Crhs & parEst["op"] == "~~",]
    p.TparEst <- paste("  i", a, b, i, ":  ", Clhs, " ~~ ", Crhs, " = ", format(round(TparEst["est"], digits=4), nsmall=4, scientific=FALSE),
                       ", p-value = ", format(round(TparEst["pvalue"], digits=4), nsmall=4, scientific=FALSE), sep="")
    cat("\n", p.TparEst)
    SumEst <- SumEst + TparEst["est"]
  } # end (for i)
  MeanEst <- SumEst/(no.waves)
  p.MeanEst <- paste("  Average across waves from T1 = ", format(round(MeanEst, digits=4), nsmall=4, scientific=FALSE), sep="")
  cat("\n", p.MeanEst)

  # -- Save pairs of non-invariant co-movements -- #
  NI.path <- t(matrix(0,ncol=no.compare, nrow=2))
  k = 1
  for (j in 2:no.waves) {
    for (i in j:no.waves) {
      Clhs <- paste("i", a, b, i, "-i", a, b, i-j+1, sep="")
      if (pest2[Clhs,3] < p) {
        NI.path[k, 1] <- i
        NI.path[k, 2] <- i-j+1
        k <- k + 1
      } # end (if pest2)
    } # end (for i)
  } # end (for j)
  Count.NI.path = k - 1

  ## -- Select sets of invariant residual covariances and print -- ##
  cat(rep("\n",2), paste("# -- Sets of invariant co-movements i", a, b, " coefficients -- #", sep=""))

  for (k in 0:MIset) {
    Noset <- no.waves - k
    NIset <- factorial(no.waves)/(factorial(no.waves-k)*factorial(k))
    mm <- t(combn((no.waves), Noset))
    for (j in 1:Count.NI.path) {
      p.j <- NI.path[j,1]
      q.j <- NI.path[j,2]
      for (i in 1:NIset) {
        count4 <- length(which(mm[i,] == p.j | mm[i,] == q.j))
        if (count4 > 1) {
          mm[i,] <- 0
        } # end (if count4)
      } # end (for i)
    } # end (for j)
    for (ii in 1:NIset) {
      count4 <- sum(mm[ii,])
      if (count4 > 0) {
        mm.p <- mm
        for (i in 1:no.waves) {
          mm.p[mm.p == i] <- paste("i", a, b, i, sep="")
        } # end (for i)
        cat("\n", "    ", mm.p[ii,])
      } # end (if count4)
    } # end (for ii)
  } # end (for k)
}  # end (Function LandD_iXY)
## ---------------------------------------------------------------- ##





## ----- Sub-Function List and Delete for Residual Covariance of indicators ----- ##
LandD_eIXY <- function(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X", b="Y") {
if (any(parEst[,4] == "eIXY2")) {
  SumEst <- 0
  cat(rep("\n",3), paste("# -- Indicator residual covariance of cov", a, b, " coefficients -- #", sep=""))
  aa <- get(a)
  bb <- get(b)
  for (i in 1: no.path) {
    Clhs <- paste(aa, i+1, sep="")
    Crhs <- paste(bb, i+1, sep="")
    TparEst <- parEst[parEst["lhs"] == Clhs & parEst["rhs"] == Crhs & parEst["op"] == "~~",]
    p.TparEst <- paste("  eI", a, b, i+1, ":  ", Clhs, " ~~ ", Crhs, " = ", format(round(TparEst["est"], digits=4), nsmall=4, scientific=FALSE),
                       ", p-value = ", format(round(TparEst["pvalue"], digits=4), nsmall=4, scientific=FALSE), sep="")
    cat("\n", p.TparEst)
    SumEst <- SumEst + TparEst["est"]
  } # end (for i)
  MeanEst <- SumEst/(no.path)
  p.MeanEst <- paste("  Average across waves from T2 = ", format(round(MeanEst, digits=4), nsmall=4, scientific=FALSE), sep="")
  cat("\n", p.MeanEst)

  # -- Save pairs of non-invariant residual covariances -- #
  NI.path <- t(matrix(0,ncol=no.compare, nrow=2))
  k = 1
  for (j in 3:no.waves) {
    for (i in j:no.waves) {
      Clhs <- paste("eI", a, b, i, "-eI", a, b, i-j+2, sep="")
      if (pest2[Clhs,3] < p) {
        NI.path[k, 1] <- i
        NI.path[k, 2] <- i-j+2
        k <- k + 1
      } # end (if pest2)
    } # end (for i)
  } # end (for j)
  Count.NI.path = k - 1

  ## -- Select sets of invariant residual covariances and print -- ##
  cat(rep("\n",2), paste("# -- Sets of invariant indicator residual covariance eI", a, b, " coefficients -- #", sep=""))

  for (k in 0:MIset) {
    Noset <- no.path - k
    NIset <- factorial(no.path)/(factorial(no.path-k)*factorial(k))
    mm <- t(combn((no.path), Noset))
    for (j in 1:Count.NI.path) {
      p.j <- NI.path[j,1]
      q.j <- NI.path[j,2]
      for (i in 1:NIset) {
        count4 <- length(which(mm[i,] == p.j | mm[i,] == q.j))
        if (count4 > 1) {
          mm[i,] <- 0
        } # end (for count4)
      } # end (for i)
    } # end (for j)
    for (ii in 1:NIset) {
      count4 <- sum(mm[ii,])
      if (count4 > 0) {
        mm.p <- mm
        for (i in 1:no.waves) {
          mm.p[mm.p == i] <- paste("eI", a, b, i+1, sep="")
        } # end (for i)
        cat("\n", "    ", mm.p[ii,])
      } # end (if count4)
    } # end (for ii)
  } # end (for k)

} else if (any(parEst[,4] == "eIXY1")) {

  SumEst <- 0
  cat(rep("\n",3), paste("# -- Indicator residual covariance of cov", a, b, " coefficients -- #", sep=""))
  aa <- get(a)
  bb <- get(b)
  for (i in 1: no.waves) {
    Clhs <- paste(aa, i, sep="")
    Crhs <- paste(bb, i, sep="")
    TparEst <- parEst[parEst["lhs"] == Clhs & parEst["rhs"] == Crhs & parEst["op"] == "~~",]
    p.TparEst <- paste("  eI", a, b, i, ":  ", Clhs, " ~~ ", Crhs, " = ", format(round(TparEst["est"], digits=4), nsmall=4, scientific=FALSE),
                       ", p-value = ", format(round(TparEst["pvalue"], digits=4), nsmall=4, scientific=FALSE), sep="")
    cat("\n", p.TparEst)
    SumEst <- SumEst + TparEst["est"]
  } # end (for i)
  MeanEst <- SumEst/(no.waves)
  p.MeanEst <- paste("  Average across waves from T1 = ", format(round(MeanEst, digits=4), nsmall=4, scientific=FALSE), sep="")
  cat("\n", p.MeanEst)

  # -- Save pairs of non-invariant residual covariances -- #
  NI.path <- t(matrix(0,ncol=no.compare, nrow=2))
  k = 1
  for (j in 2:no.waves) {
    for (i in j:no.waves) {
      Clhs <- paste("eI", a, b, i, "-eI", a, b, i-j+1, sep="")
      if (pest2[Clhs,3] < p) {
        NI.path[k, 1] <- i
        NI.path[k, 2] <- i-j+1
        k <- k + 1
      } # end (if pest2)
    } # end (for i)
  } # end (for j)
  Count.NI.path = k - 1

  ## -- Select sets of invariant residual covariances and print -- ##
  cat(rep("\n",2), paste("# -- Sets of invariant indicator residual covariance eI", a, b, " coefficients -- #", sep=""))

  for (k in 0:MIset) {
    Noset <- no.waves - k
    NIset <- factorial(no.waves)/(factorial(no.waves-k)*factorial(k))
    mm <- t(combn((no.waves), Noset))
    for (j in 1:Count.NI.path) {
      p.j <- NI.path[j,1]
      q.j <- NI.path[j,2]
      for (i in 1:NIset) {
        count4 <- length(which(mm[i,] == p.j | mm[i,] == q.j))
        if (count4 > 1) {
          mm[i,] <- 0
        } # end (for count4)
      } # end (for i)
    } # end (for j)
    for (ii in 1:NIset) {
      count4 <- sum(mm[ii,])
      if (count4 > 0) {
        mm.p <- mm
        for (i in 1:no.waves) {
          mm.p[mm.p == i] <- paste("eI", a, b, i, sep="")
        } # end (for i)
        cat("\n", "    ", mm.p[ii,])
      } # end (if count4)
    } # end (for ii)
  } # end (for k)
} # end (if eIXY)
} # end (Function LandD_eIXY)
## ---------------------------------------------------------------- ##





## ----- Sub-Function List and Delete for Intercepts ----- ##
LandD_MEAN <- function(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X") {

  SumEst <- 0
  cat(rep("\n",3), paste("# -- Intercept of w", a, " -- #", sep=""))
  aa <- get(a)
  for (i in 1: no.path) {
    Clhs <- paste0("w", aa, i+1, sep="")
    TparEst <- parEst[parEst["lhs"] == Clhs & parEst["op"] == "~1",]
    p.TparEst <- paste0("  Mw", a, i+1, ":  Intercept of ", Clhs, " = ", format(round(TparEst["est"], digits=4), nsmall=4, scientific=FALSE),
                       ", p-value = ", format(round(TparEst["pvalue"], digits=4), nsmall=4, scientific=FALSE), sep="")
    cat("\n", p.TparEst)
    SumEst <- SumEst + TparEst["est"]
  } # end (for i)
  MeanEst <- SumEst/no.path
  p.MeanEst <- paste("  Average intercepts across waves from T2 = ", format(round(MeanEst, digits=4), nsmall=4, scientific=FALSE), sep="")
  cat("\n", p.MeanEst)

  # -- Save pairs of non-invariant Intercepts -- #
  NI.path <- t(matrix(0,ncol=no.compare, nrow=2))
  k = 1
  for (j in 3:no.waves) {
    for (i in j:no.waves) {
      Clhs <- paste("Mw", a, i, "-Mw", a, i-j+2, sep="")
      if (pest2[Clhs,3] < p) {
        NI.path[k, 1] <- i-1
        NI.path[k, 2] <- i-j+1
        k <- k + 1
      } # end (if pest2)
    } # end (for i)
  } # end (for j)
  Count.NI.path = k - 1

  # -- Select sets of invariant intercepts and print -- #
  cat(rep("\n",2), paste("# -- Sets of invariant intercepts ", aa, " -- #", sep=""))

  for (k in 0:MIset) {
    Noset <- no.path - k
    NIset <- factorial(no.path)/(factorial(no.path-k)*factorial(k))
    mm <- t(combn((no.path), Noset))
    for (j in 1:Count.NI.path) {
      p.j <- NI.path[j,1]
      q.j <- NI.path[j,2]
      for (i in 1:NIset) {
        count4 <- length(which(mm[i,] == p.j | mm[i,] == q.j))
        if (count4 > 1) {
          mm[i,] <- 0
        } # end (if count4)
      } # end (for i)
    } # end (for j)
    for (ii in 1:NIset) {
      count4 <- sum(mm[ii,])
      if (count4 > 0) {
        mm.p <- mm
        for (i in 1:no.path) {
          mm.p[mm.p == i] <- paste("Mw", a, i+1, sep="")
        } # end (for i)
        cat("\n", "    ", mm.p[ii,])
      } # end (if count4)
    } # end (for ii)
  } # end (for k)
}  # end (Function LandD_MEAN)
## ---------------------------------------------------------------- ##




## ----- Sub-Function List and Delete for Indicator Variance ----- ##
LandD_eIXX <- function(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X") {

if (any(parEst[,4] == "eIXX2")) {
  SumEst <- 0
  cat(rep("\n",3), paste("# -- Indicator Variance of ", a, " -- #", sep=""))
  aa <- get(a)
  for (i in 1: no.path) {
    Clhs <- paste(aa, i+1, sep="")
    TparEst <- parEst[parEst["lhs"] == Clhs & parEst["rhs"] == Clhs & parEst["op"] == "~~",]
    p.TparEst <- paste("  eI", a, a, i+1, ":  Indicator Variance of ", Clhs, " = ", format(round(TparEst["est"], digits=4), nsmall=4, scientific=FALSE),
                       ", p-value = ", format(round(TparEst["pvalue"], digits=4), nsmall=4, scientific=FALSE), sep="")
    cat("\n", p.TparEst)
    SumEst <- SumEst + TparEst["est"]
  } # end (for i)
  MeanEst <- SumEst/no.path
  p.MeanEst <- paste("  Average across waves from T2 = ", format(round(MeanEst, digits=4), nsmall=4, scientific=FALSE), sep="")
  cat("\n", p.MeanEst)

  # -- Save pairs of non-invariant Indicator Variance -- #
  NI.path <- t(matrix(0,ncol=no.compare, nrow=2))
  k = 1
  for (j in 3:no.waves) {
    for (i in j:no.waves) {
      Clhs <- paste("eI", a, a, i, "-eI", a, a, i-j+2, sep="")
      if (pest2[Clhs,3] < p) {
        NI.path[k, 1] <- i-1
        NI.path[k, 2] <- i-j+2
        k <- k + 1
      } # end (if pest2)
    } # end (for i)
  } # end (for j)
  Count.NI.path = k - 1

  # -- Select sets of invariant indicator variance and print -- #
  cat(rep("\n",2), paste("# -- Sets of invariant Indicator Variance ", aa, " -- #", sep=""))

  for (k in 0:MIset) {
    Noset <- no.path - k
    NIset <- factorial(no.path)/(factorial(no.path-k)*factorial(k))
    mm <- t(combn((no.path), Noset))
    for (j in 1:Count.NI.path) {
      p.j <- NI.path[j,1]
      q.j <- NI.path[j,2]
      for (i in 1:NIset) {
        count4 <- length(which(mm[i,] == p.j | mm[i,] == q.j))
        if (count4 > 1) {
          mm[i,] <- 0
        } # end (if count4)
      } # end (for i)
    } # end (for j)
    for (ii in 1:NIset) {
      count4 <- sum(mm[ii,])
      if (count4 > 0) {
        mm.p <- mm
        for (i in 1:no.waves) {
          mm.p[mm.p == i] <- paste("eI", a, a, i+1, sep="")
        } # end (for i)
        cat("\n", "    ", mm.p[ii,])
      } # end (if count4)
    } # end (for ii)
  } # end (for k)

} else if (any(parEst[,4] == "eIXX1")) {

  SumEst <- 0
  cat(rep("\n",3), paste("# -- Indicator Variance of ", a, " -- #", sep=""))
  aa <- get(a)
  for (i in 1: no.waves) {
    Clhs <- paste(aa, i, sep="")
    TparEst <- parEst[parEst["lhs"] == Clhs & parEst["rhs"] == Clhs & parEst["op"] == "~~",]
    p.TparEst <- paste("  eI", a, a, i, ":  Indicator Variance of ", Clhs, " = ", format(round(TparEst["est"], digits=4), nsmall=4, scientific=FALSE),
                       ", p-value = ", format(round(TparEst["pvalue"], digits=4), nsmall=4, scientific=FALSE), sep="")
    cat("\n", p.TparEst)
    SumEst <- SumEst + TparEst["est"]
  } # end (for i)
  MeanEst <- SumEst/no.waves
  p.MeanEst <- paste("  Average across waves = ", format(round(MeanEst, digits=4), nsmall=4, scientific=FALSE), sep="")
  cat("\n", p.MeanEst)

  # -- Save pairs of non-invariant Indicator Variance -- #
  NI.path <- t(matrix(0,ncol=no.compare, nrow=2))
  k = 1
  for (j in 2:no.waves) {
    for (i in j:no.waves) {
      Clhs <- paste("eI", a, a, i, "-eI", a, a, i-j+1, sep="")
      if (pest2[Clhs,3] < p) {
        NI.path[k, 1] <- i
        NI.path[k, 2] <- i-j+1
        k <- k + 1
      } # end (if pest2)
    } # end (for i)
  } # end (for j)
  Count.NI.path = k - 1

  # -- Select sets of invariant indicator variance and print -- #
  cat(rep("\n",2), paste("# -- Sets of invariant Indicator Variance ", aa, " -- #", sep=""))

  for (k in 0:MIset) {
    Noset <- no.waves - k
    NIset <- factorial(no.waves)/(factorial(no.waves-k)*factorial(k))
    mm <- t(combn((no.waves), Noset))
    for (j in 1:Count.NI.path) {
      p.j <- NI.path[j,1]
      q.j <- NI.path[j,2]
      for (i in 1:NIset) {
        count4 <- length(which(mm[i,] == p.j | mm[i,] == q.j))
        if (count4 > 1) {
          mm[i,] <- 0
        } # end (if count4)
      } # end (for i)
    } # end (for j)
    for (ii in 1:NIset) {
      count4 <- sum(mm[ii,])
      if (count4 > 0) {
        mm.p <- mm
        for (i in 1:no.waves) {
          mm.p[mm.p == i] <- paste("eI", a, a, i, sep="")
        } # end (for i)
        cat("\n", "    ", mm.p[ii,])
      } # end (if count4)
    } # end (for ii)
  } # end (for k)
} # end (if eIXX)
}  # end (Sub-Function LandD_eIXX)
## ---------------------------------------------------------------- ##





# ==================== Creating Function "CLPM" ==================== #
#' Function CLPM (Cross-Lagged Panel Model)
#'
#' Cross-Lagged Panel Model (CLPM)
#'
#' @param data.source name of data.frame
#' @param no.waves number of waves (minimum = 3, must be grater than lag)
#' @param lag number of waves between two lags (minimum = 1, maximum = 4)
#' @param Type1 Overall Type I error rate (default is 0.05) for comparing estimated parameters in the List and Delete method.
#' @param Type1Adj Adjustment of Type I error rate for multiple tests for each estimated parameter. Default is "BON" (Bonferroni adjustment -- Type I error rate/no. of pairwise comparisons), can also be "NULL" (without adjustment).
#' @param X name of variable X.
#' @param Y name of variable Y.
#' @param Z name of variable Z.
#' @param W name of variable W (Z must come before W).
#'
#' @return Cross-Lagged Panel Model (CLPM) outputs.
#' @export
#' @examples
#'
#' ## -- Example -- ##
#'
#' CLPM(data.source="Data_A", 7, 2, Type1=0.05, Type1Adj="BON", X="EXPOSE.", Y="INTENS.")
#'

CLPM <- function(data.source, no.waves, lag=1, Type1=0.05, Type1Adj="BON", X, Y, Z="NULL", W = "NULL") {

  ## -- Check inputs -- ##

  if (no.waves < 3) stop("Minimum number of waves is 3")

  if (lag < 1) stop("Minimum number of lag is 1")
  if (lag > 4) stop("Maximum number of lag is 4")

  if ((no.waves - lag) < 1) stop("Number of waves must be greater than lag plus 1")

  if (Type1 > 0.05) stop("Type I error rate > 0.05 is not recommended")
  if (Type1 < 0.0001) stop("Type I error rate < 0.0001 is not recommended")

  if (Z == "NULL" & W != "NULL") stop("Z must be defined before W")

  Type1Adj <- toupper(Type1Adj)
  match.arg(Type1Adj, c("BON", "NULL"))


  ## ----- Creating Model CLPM ----- ##
  sink('CLPM.txt')
    cat("\n", "# == Specify the model (CLPM) == #", "\n")
    cat("\n", "CLPM <- '")

    cat(rep("\n",2), "  # -- Create latent variables -- #")
    for (i in 1:no.waves) {
      cat("\n", paste("  w", X, i, " =~ 1*", X, i, sep=""))
      cat("\n", paste("  w", Y, i, " =~ 1*", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  w", Z, i, " =~ 1*", Z, i, sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  w", W, i, " =~ 1*", W, i, sep=""))
      } # end (if W)
    } # end (for i)

    cat(rep("\n",2), "  # -- Constrain residual variance of indicators to zero -- #")
    for (i in 1:no.waves) {
      cat("\n", paste("  ", X, i, " ~~ 0*", X, i, sep=""))
      cat("\n", paste("  ", Y, i, " ~~ 0*", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  ", Z, i, " ~~ 0*", Y, i, sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  ", W, i, " ~~ 0*", W, i, sep=""))
      } # end (if W)
    } # end (for i)

    cat(rep("\n",2), "  # -- Estimate variance of latent variables at first wave -- #")
    cat("\n", paste("  w", X, "1 ~~ eXX1*w", X, "1", sep=""))
    cat("\n", paste("  w", Y, "1 ~~ eYY1*w", Y, "1", sep=""))
    if (Z != "NULL") {
      cat("\n", paste("  w", Z, "1 ~~ eZZ1*w", Z, "1", sep=""))
    } # end (if Z)
    if (W != "NULL") {
      cat("\n", paste("  w", W, "1 ~~ eWW1*w", W, "1", sep=""))
    } # end (if W)

    cat(rep("\n",2), "  # -- Estimate residual variance of latent variables -- #")
    cat("\n", "  ###################################################################")
    cat("\n", "  # Remove the subscripts for eXX for invariant residual covariance #")
    cat("\n", "  ###################################################################")
    for (i in 2:no.waves) {
      cat("\n", paste("  w", X, i, " ~~ eXX", i, "*w", X, i, sep=""))
      cat("\n", paste("  w", Y, i, " ~~ eYY", i, "*w", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  w", Z, i, " ~~ eZZ", i, "*w", Z, i, sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  w", W, i, " ~~ eWW", i, "*w", W, i, sep=""))
      } # end (if W)
    } # end (for i)

    cat(rep("\n",2), "  # -- Estimate covariances between residuals of latent variables -- #")
    cat("\n", "  ###################################################################")
    cat("\n", "  # Remove the subscripts for eXY for invariant residual covariance #")
    cat("\n", "  ###################################################################")
    for (i in 2:no.waves) {
      cat("\n", paste("  w", X, i, " ~~ eXY", i, "*w", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  w", X, i, " ~~ eXZ", i, "*w", Z, i, sep=""))
        cat("\n", paste("  w", Y, i, " ~~ eYZ", i, "*w", Z, i, sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  w", X, i, " ~~ eXW", i, "*w", W, i, sep=""))
        cat("\n", paste("  w", Y, i, " ~~ eYW", i, "*w", W, i, sep=""))
        cat("\n", paste("  w", Z, i, " ~~ eZW", i, "*w", W, i, sep=""))
      } # end (if W)
    } # end ((for i)

    cat(rep("\n",2), "  # -- Estimate covariance between latent variables at first wave -- #")
    cat("\n", "   w", X, "1 ~~ w", Y, "1", sep="")
    if (Z != "NULL") {
      cat("\n", "    w", X, "1 ~~ w", Z, "1", sep="")
      cat("\n", "    w", Y, "1 ~~ w", Z, "1", sep="")
    } # end (if Z)
    if (W != "NULL") {
      cat("\n", "    w", X, "1 ~~ w", W, "1", sep="")
      cat("\n", "    w", Y, "1 ~~ w", W, "1", sep="")
      cat("\n", "    w", Z, "1 ~~ w", W, "1", sep="")
    } # end (if W)

    for (j in 1:lag) {
      cat(rep("\n",2), paste0("  # -- Estimate lagged effects between latent variables (Lag = ", j, " wave) -- #", sep =""))
      cat("\n", "  #############################################")
      cat("\n", "  # Remove the subscripts for invariant paths #")
      cat("\n", "  #############################################")
      for (i in (j+1):no.waves) {
        if (Z == "NULL") {
          cat("\n", paste("  w", X,i, " ~ pXX", i,i-j, "*w", X,i-j, " + pYX", i,i-j, "*w", Y,i-j, sep=""))
          cat("\n", paste("  w", Y,i, " ~ pXY", i,i-j, "*w", X,i-j, " + pYY", i,i-j, "*w", Y,i-j, sep=""))
        } else if (W != "NULL") {
          cat("\n", paste("  w", X,i, " ~ pXX", i,i-j, "*w", X,i-j, " + pYX", i,i-j, "*w", Y,i-j, " + pZX", i,i-j, "*w", Z,i-j, " + pWX", i,i-j, "*w", W,i-j, sep=""))
          cat("\n", paste("  w", Y,i, " ~ pXY", i,i-j, "*w", X,i-j, " + pYY", i,i-j, "*w", Y,i-j, " + pZY", i,i-j, "*w", Z,i-j, " + pWY", i,i-j, "*w", W,i-j, sep=""))
          cat("\n", paste("  w", Z,i, " ~ pXZ", i,i-j, "*w", X,i-j, " + pYZ", i,i-j, "*w", Y,i-j, " + pZZ", i,i-j, "*w", Z,i-j, " + pWZ", i,i-j, "*w", W,i-j, sep=""))
          cat("\n", paste("  w", W,i, " ~ pXW", i,i-j, "*w", X,i-j, " + pYW", i,i-j, "*w", Y,i-j, " + pZW", i,i-j, "*w", Z,i-j, " + pWW", i,i-j, "*w", W,i-j, sep=""))
        } else if (Z != "NULL") {
          cat("\n", paste("  w", X,i, " ~ pXX", i,i-j, "*w", X,i-j, " + pYX", i,i-j, "*w", Y,i-j, " + pZX", i,i-j, "*w", Z,i-j, sep=""))
          cat("\n", paste("  w", Y,i, " ~ pXY", i,i-j, "*w", X,i-j, " + pYY", i,i-j, "*w", Y,i-j, " + pZY", i,i-j, "*w", Z,i-j, sep=""))
          cat("\n", paste("  w", Z,i, " ~ pXZ", i,i-j, "*w", X,i-j, " + pYZ", i,i-j, "*w", Y,i-j, " + pZZ", i,i-j, "*w", Z,i-j, sep=""))
        } # end (if Z)
      } # end (for i)
    } # end (for j)

    cat(rep("\n",2), "  # -- Estimate grand means (intercepts) of indicators -- #")
    for (i in 1:no.waves) {
      cat("\n", paste("  ", X, i, " ~ Mw", X, i, "*1", sep=""))
      cat("\n", paste("  ", Y, i, " ~ Mw", Y, i, "*1", sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  ", Z, i, " ~ Mw", Z, i, "*1", sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  ", W, i, " ~ Mw", W, i, "*1", sep=""))
      } # end (if W)
    } # end (for i)


    cat(rep("\n",2), "  ##########################################")
    cat("\n", "  # Regression of indicators on C1 #")
    cat("\n", "  ##########################################")
    for (i in 1:no.waves) {
      cat("\n", paste("  #  ", X, i, " ~ sx", i,"*C1", sep=""))
      cat("\n", paste("  #  ", Y, i, " ~ sy", i,"*C1", sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  #  ", Z, i, " ~ sz", i,"*C1", sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  #  ", W, i, " ~ sz", i,"*C1", sep=""))
      } # end (if W)
    } # end (for i)

    cat(rep("\n",2), "  ################################################")
    cat("\n", "  # Regression of outcome D1 on latent variables #")
    cat("\n", "  ################################################")
    for (i in 1:no.waves) {
      cat("\n", paste("  #  D1 ~ bx", i, "*w", X, i, sep=""))
      cat("\n", paste("  #  D1 ~ by", i, "*w", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  #  D1 ~ bz", i, "*w", Z, i, sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  #  D1 ~ bw", i, "*w", W, i, sep=""))
      } # end (if W)
    } # end (for i)
    cat("\n", "  #  D1 ~~ D1 # Residual variance of D1")

    cat(rep("\n",2), "  '")

    # -- Run CLPMMLR -- #
    cat(rep("\n",2), "  CLPMMLR.fit <- lavaan::sem(CLPM, ")
    cat("\n", "  ", data.source, ",")
    cat("\n", "   missing = 'fiml',")
    cat("\n", "   meanstructure = TRUE,")
    cat("\n", "   information = 'observed',")
    cat("\n", "   estimator = 'MLR'")
    cat("\n", ")")

    # Request summary outputs
    cat(rep("\n",2), "  print(lavaan::summary(CLPMMLR.fit, fit.measure = TRUE, standardized = TRUE, rsq = TRUE))", rep("\n",3))

  sink() # Stop writing to file

  ## -- Execute CLPM.txt and request summary outputs-- ##
  source('CLPM.txt')
  if (lavaan::lavInspect(CLPMMLR.fit, what ="post.check") == FALSE) stop("The lavaan solution is non-admissible.")
  ## ------------------------------- ##


  ## -- Monte Carlo Simulation -- ##

  parEst <- lavaan::parameterEstimates(CLPMMLR.fit, remove.nonfree = TRUE)
  pest2 <- parEst[, "est"]  # Estimated Parameters
  pest3 <- lavaan::lavTech(CLPMMLR.fit, what = "vcov", add.labels = TRUE)  # Estimated Variance-Covariance of Estimated Parameters

  Invariance(parEst, pest2, pest3, no.path, no.waves, lag, Type1, Type1Adj, X, Y, Z, W)

  cat(rep("\n", 2))

}  # end (Function CLPM)

## ========================================================================================== ##





# ==================== Creating Function "RICLPM" ==================== #
#' Function RICLPM (Random Intercept Cross-Lagged Panel Model)
#'
#' Random Intercept Cross-Lagged Panel Model (RI-CLPM)
#'
#' @param data.source name of data.frame
#' @param no.waves number of waves (minimum = 3, must be grater than lag)
#' @param lag number of waves between two lags (minimum = 1, maximum = 4)
#' @param Type1 Overall Type I error rate (default is 0.05) for comparing estimated parameters in the List and Delete method.
#' @param Type1Adj Adjustment of Type I error rate for multiple tests for each estimated parameter. Default is "BON" (Bonferroni adjustment -- Type I error rate/no. of pairwise comparisons), can also be "NULL" (without adjustment).
#' @param X name of variable X.
#' @param Y name of variable Y.
#' @param Z name of variable Z.
#' @param W name of variable W (Z must come before W).
#'
#' @return Random Intercept Cross-Lagged Panel Model (RI-CLPM) outputs.
#' @export
#' @examples
#'
#' ## -- Example -- ##
#'
#' RICLPM(data.source="Data_A", 7, 2, Type1=0.05, Type1Adj="BON", X="EXPOSE.", Y="INTENS.")
#'

RICLPM <- function(data.source, no.waves, lag=1, Type1=0.05, Type1Adj="BON", X, Y, Z="NULL", W = "NULL") {

  ## -- Check inputs -- ##

  if (no.waves < 3) stop("Minimum number of waves is 3")

  if (lag < 1) stop("Minimum number of lag is 1")
  if (lag > 4) stop("Maximum number of lag is 4")

  if ((no.waves - lag) < 1) stop("Number of waves must be greater than lag plus 1")

  if (Type1 > 0.05) stop("Type I error rate > 0.05 is not recommended")
  if (Type1 < 0.0001) stop("Type I error rate < 0.0001 is not recommended")

  if (Z == "NULL" & W != "NULL") stop("Z must be defined before W")

  Type1Adj <- toupper(Type1Adj)
  match.arg(Type1Adj, c("BON", "NULL"))

  ## ----- Creating Model RICLPM ----- ###
  sink('RICLPM.txt')
    cat("\n", "# Specify the model (RICLPM)", "\n")
    cat("\n", "RICLPM <- '")

    # -- Create Between Components (Random Intercepts) -- #
    cat(rep("\n",2), "  # -- Create between components (random intercepts) -- #")
    BX <- paste("  RI", X, " =~ 1*", X, "1", sep="")
    BY <- paste("  RI", Y, " =~ 1*", Y, "1", sep="")
    for (i in 2:no.waves) {
      BX <- paste(BX, " +1*", X, i, sep="")
      BY <- paste(BY, " +1*", Y, i, sep="")
    } # end (for i)
    cat("\n", BX)
    cat("\n", BY)
    if (Z != "NULL") {
      BZ <- paste("  RI", Z, " =~ 1*", Z, "1", sep="")
      for (i in 2:no.waves) {
        BZ <- paste(BZ, " +1*", Z, i, sep="")
      } # end (for i)
      cat("\n", BZ)
    } # end (if Z)
    if (W != "NULL") {
      BW <- paste("  RI", W, " =~ 1*", W, "1", sep="")
      for (i in 2:no.waves) {
        BW <- paste(BW, " +1*", W, i, sep="")
      } # end (for i)
      cat("\n", BW)
    } # end (if W)

    # -- Constrain Residual Variance of Indicators to Zero -- #
    cat(rep("\n",2), "  # -- Constrain residual variance of indicators to zero -- #")
    for (i in 1:no.waves) {
      cat("\n", paste("  ", X, i, " ~~ 0*", X, i, sep=""))
      cat("\n", paste("  ", Y, i, " ~~ 0*", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  ", Z, i, " ~~ 0*", Z, i, sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  ", W, i, " ~~ 0*", W, i, sep=""))
      } # end (if W)
    } # end (for i)

    # -- Constrain Intercepts of Indicators to Zero -- #
    cat(rep("\n",2), "  # -- Constrain intercepts of indicators to zero -- #")
    for (i in 1:no.waves) {
      cat("\n", paste("  ", X, i, " ~ 0*1", sep=""))
      cat("\n", paste("  ", Y, i, " ~ 0*1", sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  ", Z, i, " ~ 0*1", sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  ", W, i, " ~ 0*1", sep=""))
      } # (if W)
    } # end (for i)

    # -- Constrain Means (Intercepts) of Random Intercepts to Zero -- #
    cat(rep("\n",2), "  # -- Constrain means (intercepts) of random intercepts to zero -- #")
    cat("\n", paste("  RI", X, " ~ 0*1", sep=""))
    cat("\n", paste("  RI", Y, " ~ 0*1", sep=""))
    if (Z != "NULL") {
      cat("\n", paste("  RI", Z, " ~ 0*1", sep=""))
    } # end (if Z)
    if (W != "NULL") {
      cat("\n", paste("  RI", W, " ~ 0*1", sep=""))
    } # (if W)

    # -- Create Latent Variables from Indicators -- #
    cat(rep("\n",2), "  # -- Create latent variables -- #")
    for (i in 1:no.waves) {
      cat("\n", paste("  w", X, i, " =~ 1*", X, i, sep=""))
      cat("\n", paste("  w", Y, i, " =~ 1*", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  w", Z, i, " =~ 1*", Z, i, sep=""))
      }  # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  w", W, i, " =~ 1*", W, i, sep=""))
      }  # end (if W)
    } # end (for i)

    # -- Estimate Means (Intercepts) of Latent Variables -- #
    cat(rep("\n",2), "  # -- Estimate means (intercepts) of latent variables -- #")
    for (i in 1:no.waves) {
      cat("\n", paste("  w", X, i, " ~ MwX", i, "*1", sep=""))
      cat("\n", paste("  w", Y, i, " ~ MwY", i, "*1", sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  w", Z, i, " ~ MwZ", i, "*1", sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  w", W, i, " ~ MwW", i, "*1", sep=""))
      } # (if W)
    } # end (for i)

    # -- Estimate (Residual) Variance of Latent Variables -- #
    cat(rep("\n",2), "  # -- Estimate (residual) variance of latent variables -- #")
    for (i in 1:no.waves) {
      cat("\n", paste("  w", X, i, " ~~ eXX", i, "*w", X, i, sep=""))
      cat("\n", paste("  w", Y, i, " ~~ eYY", i, "*w", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  w", Z, i, " ~~ eZZ", i, "*w", Z, i, sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  w", W, i, " ~~ eWW", i, "*w", W, i, sep=""))
      } # end (if W)
    } # end (for i)

    # -- Estimate Covariances Between Residuals of Latent Variables -- #
    cat(rep("\n",2), "  # -- Estimate covariances between residuals of latent variables -- #")
    cat("\n", "  ###################################################################")
    cat("\n", "  # Remove the subscripts for eXY for invariant residual covariance #")
    cat("\n", "  ###################################################################")
    for (i in 2:no.waves) {
      cat("\n", paste("  w", X, i, " ~~ eXY", i, "*w", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  w", X, i, " ~~ eXZ", i, "*w", Z, i, sep=""))
        cat("\n", paste("  w", Y, i, " ~~ eYZ", i, "*w", Z, i, sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  w", X, i, " ~~ eXW", i, "*w", W, i, sep=""))
        cat("\n", paste("  w", Y, i, " ~~ eYW", i, "*w", W, i, sep=""))
        cat("\n", paste("  w", Z, i, " ~~ eZW", i, "*w", W, i, sep=""))
      } # end (if W)
    } # end ((for i)

    # -- Estimate Covariance Between Latent Variables at First Wave -- #
    cat(rep("\n",2), "  # -- Estimate covariance between latent variables at first wave -- #")
    cat("\n", "   w", X, "1 ~~ w", Y, "1", sep="")
    if (Z != "NULL") {
      cat("\n", "    w", X, "1 ~~ w", Z, "1", sep="")
      cat("\n", "    w", Y, "1 ~~ w", Z, "1", sep="")
    } # end (if Z)
    if (W != "NULL") {
      cat("\n", "    w", X, "1 ~~ w", W, "1", sep="")
      cat("\n", "    w", Y, "1 ~~ w", W, "1", sep="")
      cat("\n", "    w", Z, "1 ~~ w", W, "1", sep="")
    } # end (if W)

    # -- Estimate Variance and Covariance of Random Intercepts -- #
    cat(rep("\n",2), "  # -- Estimate variance and covariance of random intercepts -- #")
    cat("\n", "   RI", X, " ~~ RI", X, sep="")
    cat("\n", "   RI", Y, " ~~ RI", Y, sep="")
    cat("\n", "   RI", X, " ~~ RI", Y, sep="")
    if (Z != "NULL") {
      cat("\n", "   RI", Z, " ~~ RI", Z, sep="")
      cat("\n", "   RI", X, " ~~ RI", Z, sep="")
      cat("\n", "   RI", Y, " ~~ RI", Z, sep="")
    } # end (if Z != "NULL")
    if (W != "NULL") {
      cat("\n", "   RI", W, " ~~ RI", W, sep="")
      cat("\n", "   RI", X, " ~~ RI", W, sep="")
      cat("\n", "   RI", Y, " ~~ RI", W, sep="")
      cat("\n", "   RI", Z, " ~~ RI", W, sep="")
    } # end (if W != "NULL")

    # -- Constrain covariance among RIX RIY, RIZ, RIW and wx1, wy1, wz1, ww1 -- #
    cat(rep("\n",2), "  # -- Constrain covariance among RIx, RIy, RIz, RIw and wx1, wy1, wz1, ww1 -- #")
    cat("\n", "   RI", X, " ~~ 0*w", X, "1", sep="")
    cat("\n", "   RI", X, " ~~ 0*w", Y, "1", sep="")
    cat("\n", "   RI", Y, " ~~ 0*w", X, "1", sep="")
    cat("\n", "   RI", Y, " ~~ 0*w", Y, "1", sep="")
    if (Z != "NULL") {
      cat("\n", "   RI", X, " ~~ 0*w", Z, "1", sep="")
      cat("\n", "   RI", Y, " ~~ 0*w", Z, "1", sep="")
      cat("\n", "   RI", Z, " ~~ 0*w", X, "1", sep="")
      cat("\n", "   RI", Z, " ~~ 0*w", Y, "1", sep="")
      cat("\n", "   RI", Z, " ~~ 0*w", Z, "1", sep="")
    } # end (if Z)
    if (W != "NULL") {
      cat("\n", "   RI", X, " ~~ 0*w", W, "1", sep="")
      cat("\n", "   RI", Y, " ~~ 0*w", W, "1", sep="")
      cat("\n", "   RI", Z, " ~~ 0*w", W, "1", sep="")
      cat("\n", "   RI", W, " ~~ 0*w", X, "1", sep="")
      cat("\n", "   RI", W, " ~~ 0*w", Y, "1", sep="")
      cat("\n", "   RI", W, " ~~ 0*w", Z, "1", sep="")
      cat("\n", "   RI", W, " ~~ 0*w", W, "1", sep="")
    } # end (if W)

    for (j in 1:lag) {
      cat(rep("\n",2), paste0("  # -- Estimate lagged effects between latent variables (Lag = ", j, " wave) -- #", sep =""))
      cat("\n", "  #############################################")
      cat("\n", "  # Remove the subscripts for invariant paths #")
      cat("\n", "  #############################################")
      for (i in (j+1):no.waves) {
        if (Z == "NULL") {
          cat("\n", paste("  w", X,i, " ~ pXX", i,i-j, "*w", X,i-j, " + pYX", i,i-j, "*w", Y,i-j, sep=""))
          cat("\n", paste("  w", Y,i, " ~ pXY", i,i-j, "*w", X,i-j, " + pYY", i,i-j, "*w", Y,i-j, sep=""))
        } else if (W != "NULL") {
          cat("\n", paste("  w", X,i, " ~ pXX", i,i-j, "*w", X,i-j, " + pYX", i,i-j, "*w", Y,i-j, " + pZX", i,i-j, "*w", Z,i-j, " + pWX", i,i-j, "*w", W,i-j, sep=""))
          cat("\n", paste("  w", Y,i, " ~ pXY", i,i-j, "*w", X,i-j, " + pYY", i,i-j, "*w", Y,i-j, " + pZY", i,i-j, "*w", Z,i-j, " + pWY", i,i-j, "*w", W,i-j, sep=""))
          cat("\n", paste("  w", Z,i, " ~ pXZ", i,i-j, "*w", X,i-j, " + pYZ", i,i-j, "*w", Y,i-j, " + pZZ", i,i-j, "*w", Z,i-j, " + pWZ", i,i-j, "*w", W,i-j, sep=""))
          cat("\n", paste("  w", W,i, " ~ pXW", i,i-j, "*w", X,i-j, " + pYW", i,i-j, "*w", Y,i-j, " + pZW", i,i-j, "*w", Z,i-j, " + pWW", i,i-j, "*w", W,i-j, sep=""))
        } else if (Z != "NULL") {
          cat("\n", paste("  w", X,i, " ~ pXX", i,i-j, "*w", X,i-j, " + pYX", i,i-j, "*w", Y,i-j, " + pZX", i,i-j, "*w", Z,i-j, sep=""))
          cat("\n", paste("  w", Y,i, " ~ pXY", i,i-j, "*w", X,i-j, " + pYY", i,i-j, "*w", Y,i-j, " + pZY", i,i-j, "*w", Z,i-j, sep=""))
          cat("\n", paste("  w", Z,i, " ~ pXZ", i,i-j, "*w", X,i-j, " + pYZ", i,i-j, "*w", Y,i-j, " + pZZ", i,i-j, "*w", Z,i-j, sep=""))
        } # end (if Z)
      } # end (for i)
    } # end (for j)

    cat(rep("\n",2), "  ##########################################")
    cat("\n", "  # Regression of indicators on C1 #")
    cat("\n", "  ##########################################")
    for (i in 1:no.waves) {
      cat("\n", paste("  #  ", X, i, " ~ sx", i,"*C1", sep=""))
      cat("\n", paste("  #  ", Y, i, " ~ sy", i,"*C1", sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  #  ", Z, i, " ~ sz", i,"*C1", sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  #  ", W, i, " ~ sz", i,"*C1", sep=""))
      } # end (if W)
    } # end (for i)

    cat(rep("\n",2), "  #########################################")
    cat("\n", "  # Regression of random intercepts on C1 #")
    cat("\n", "  #########################################")
    if (Z == "NULL") {
      cat("\n", "   #  RI", X, " + RI", Y, " ~ C1", sep="")
    } else if (Z != "NULL") {
      cat("\n", "   #  RI", X, " + RI", Y, " + RI", Z, " ~ C1", sep="")
    } else if (W != "NULL") {
      cat("\n", "   #  RI", X, " + RI", Y, " + RI", Z, " + RI", W, " ~ C1", sep="")
    } # end (if Z/W)

    cat(rep("\n",2), "  ################################################################")
    cat("\n", "  # Regression of time-invariant outcome D1 on random intercepts #")
    cat("\n", "  ################################################################")
    if (Z == "NULL") {
      cat("\n", "   #  D1 ~ RI", X, " + RI", Y, sep="")
    } else if (Z != "NULL") {
      cat("\n", "   #  D1 ~ RI", X, " + RI", Y, " + RI", Z, sep="")
    } else if (W != "NULL") {
      cat("\n", "   #  D1 ~ RI", X, " + RI", Y, " + RI", Z, " + RI", W, sep="")
    } # end (if Z/W)
    cat("\n", "  #  D1 ~~ D1 # Residual variance of D1")

    cat(rep("\n",2), "  ################################################")
    cat("\n", "  # Regression of outcome D1 on latent variables #")
    cat("\n", "  ################################################")
    for (i in 1:no.waves) {
      cat("\n", paste("  #  D1 ~ bx", i, "*w", X, i, sep=""))
      cat("\n", paste("  #  D1 ~ by", i, "*w", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  #  D1 ~ bz", i, "*w", Z, i, sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  #  D1 ~ bw", i, "*w", W, i, sep=""))
      } # end (if W)
    } # end (for i)
    cat("\n", "  #  D1 ~~ D1 # Residual variance of D1")

    cat(rep("\n",2), "  '")

    # -- Run RICLPMMLR -- #
    cat(rep("\n",2), "  RICLPMMLR.fit <- lavaan::sem(RICLPM, ")
    cat("\n", "   ", data.source, ",")
    cat("\n", "   missing = 'fiml',")
    cat("\n", "   meanstructure = TRUE,")
    cat("\n", "   information = 'observed',")
    cat("\n", "   estimator = 'MLR'")
    cat("\n", ")")

    # Request summary outputs
    cat(rep("\n",2), "  print(lavaan::summary(RICLPMMLR.fit, fit.measure = TRUE, standardized = TRUE, rsq = TRUE))", rep("\n", 3))

  sink() # Stop writing to file


  ## -- Execute RICLPM.txt and request summary outputs-- ##
  source('RICLPM.txt')
  if (lavaan::lavInspect(RICLPMMLR.fit, what ="post.check") == FALSE) stop("The lavaan solution is non-admissible.")
  ## ------------------------------- ##


  ## -- Monte Carlo Simulation -- ##

  parEst <- lavaan::parameterEstimates(RICLPMMLR.fit, remove.nonfree = TRUE)
  pest2 <- parEst[, "est"]  # Estimated Parameters
  pest3 <- lavaan::lavTech(RICLPMMLR.fit, what = "vcov", add.labels = TRUE)  # Estimated Variance-Covariance of Estimated Parameters

  Invariance(parEst, pest2, pest3, no.path, no.waves, lag, Type1, Type1Adj, X, Y, Z, W)

  cat(rep("\n", 2))

}  # end (Function RICLPM)

## ========================================================================================== ##





# ==================== Creating Function "LGCMSR" ==================== #
#' Function LGCMSR (Latent Growth Curve Model with Structural Residuals)
#'
#' Latent Growth Curve Model with Structural Residuals (LGCM-SR)
#'
#' @param data.source name of data.frame
#' @param no.waves number of waves (minimum = 3, must be grater than lag)
#' @param lag number of waves between two lags (minimum = 1, maximum = 4)
#' @param Type1 Overall Type I error rate (default is 0.05) for comparing estimated parameters in the List and Delete method.
#' @param Type1Adj Adjustment of Type I error rate for multiple tests for each estimated parameter. Default is "BON" (Bonferroni adjustment -- Type I error rate/no. of pairwise comparisons), can also be "NULL" (without adjustment).
#' @param X name of variable X.
#' @param Y name of variable Y.
#' @param Z name of variable Z.
#' @param W name of variable W (Z must come before W).
#'
#' @return Latent Growth Curve Model with Structural Residuals (LGCM-SR) outputs.
#' @export
#' @examples
#'
#' ## -- Example -- ##
#'
#' LGCMSR(data.source="Data_A", 7, 2, Type1=0.05, Type1Adj="BON", X="EXPOSE.", Y="INTENS.")
#'

LGCMSR <- function(data.source, no.waves, lag=1, Type1=0.05, Type1Adj="BON", X, Y, Z="NULL", W = "NULL") {

  ## -- Check inputs -- ##

  if (no.waves < 3) stop("Minimum number of waves is 3")

  if (lag < 1) stop("Minimum number of lag is 1")
  if (lag > 4) stop("Maximum number of lag is 4")

  if ((no.waves - lag) < 1) stop("Number of waves must be greater than lag plus 1")

  if (Type1 > 0.05) stop("Type I error rate > 0.05 is not recommended")
  if (Type1 < 0.0001) stop("Type I error rate < 0.0001 is not recommended")

  if (Z == "NULL" & W != "NULL") stop("Z must be defined before W")

  Type1Adj <- toupper(Type1Adj)
  match.arg(Type1Adj, c("BON", "NULL"))


  ## ----- Creating Model LGCMSR ----- ###
  sink('LGCMSR.txt')
    cat("\n", "# Specify the model (LGCMSR)", "\n")
    cat("\n", "LGCMSR <- '")

    # -- Create Between Components (Random Intercepts) -- #
    cat(rep("\n",2), "  # -- Create between components (random intercepts) -- #")
    BX <- paste("  RI", X, " =~ 1*", X, "1", sep="")
    BY <- paste("  RI", Y, " =~ 1*", Y, "1", sep="")
    for (i in 2:no.waves) {
      BX <- paste(BX, " +1*", X, i, sep="")
      BY <- paste(BY, " +1*", Y, i, sep="")
    } # end (for i)
    cat("\n", BX)
    cat("\n", BY)
    if (Z != "NULL") {
      BZ <- paste("  RI", Z, " =~ 1*", Z, "1", sep="")
      for (i in 2:no.waves) {
        BZ <- paste(BZ, " +1*", Z, i, sep="")
      } # end (for i)
      cat("\n", BZ)
    } # end (if Z)
    if (W != "NULL") {
      BW <- paste("  RI", W, " =~ 1*", W, "1", sep="")
      for (i in 2:no.waves) {
        BW <- paste(BW, " +1*", W, i, sep="")
      } # end (for i)
      cat("\n", BW)
    } # end (if W)

    # -- Create Between Components (Random Slopes) -- #
    cat(rep("\n",2), "  # -- Create between components (random slopes) -- #")
    BX <- paste("  RS", X, " =~ 0*", X, "1", sep="")
    BY <- paste("  RS", Y, " =~ 0*", Y, "1", sep="")
    for (i in 2:no.waves) {
      BX <- paste(BX, " +", (i-1), "*", X, i, sep="")
      BY <- paste(BY, " +", (i-1), "*", Y, i, sep="")
    } # end (for i)
    cat("\n", BX)
    cat("\n", BY)
    if (Z != "NULL") {
      BZ <- paste("  RS", Z, " =~ 0*", Z, "1", sep="")
      for (i in 2:no.waves) {
        BZ <- paste(BZ, " +", (i-1), "*", Z, i, sep="")
      } # end (for i)
      cat("\n", BZ)
    } # end (if Z)
    if (W != "NULL") {
      BW <- paste("  RS", W, " =~ 0*", W, "1", sep="")
      for (i in 2:no.waves) {
        BW <- paste(BW, " +", (i-1), "*", W, i, sep="")
      } # end (for i)
      cat("\n", BW)
    } # end (if W)

    # -- Constrain Residual Variance of Indicators to Zero -- #
    cat(rep("\n",2), "  # -- Constrain residual variance of indicators to zero -- #")
    for (i in 1:no.waves) {
      cat("\n", paste("  ", X, i, " ~~ 0*", X, i, sep=""))
      cat("\n", paste("  ", Y, i, " ~~ 0*", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  ", Z, i, " ~~ 0*", Z, i, sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  ", W, i, " ~~ 0*", W, i, sep=""))
      } # end (if W)
    } # end (for i)###   ###

    # -- Constrain Means (Intercepts) of Indicators to zero -- #
    cat(rep("\n",2), "  # -- Constrain means (intercepts) of indicators to zero -- #")
    for (i in 1:no.waves) {
      cat("\n", paste("  ", X, i, " ~ 0*1", sep=""))
      cat("\n", paste("  ", Y, i, " ~ 0*1", sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  ", Z, i, " ~ 0*1", sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  ", W, i, " ~ 0*1", sep=""))
      } # (if W)
    } # end (for i)

    # -- Estimate Means (Intercepts) of Random Intercepts -- #
    cat(rep("\n",2), "  # -- Estimate means (intercepts) of random intercepts -- #")
    cat("\n", paste("  RI", X, " ~ MRI", X, "*1", sep=""))
    cat("\n", paste("  RI", Y, " ~ MRI", Y, "*1", sep=""))
    if (Z != "NULL") {
      cat("\n", paste("  RI", Z, " ~ MRI", Z, "*1", sep=""))
    } # end (if Z)
    if (W != "NULL") {
      cat("\n", paste("  RI", W, " ~ MRI", W, "*1", sep=""))
    } # (if W)

    # -- Estimate Means (Intercepts) of Random Slopes -- #
    cat(rep("\n",2), "  # -- Estimate means (intercepts) of random slopes -- #")
    cat("\n", paste("  RS", X, " ~ MRS", X, "*1", sep=""))
    cat("\n", paste("  RS", Y, " ~ MRS", Y, "*1", sep=""))
    if (Z != "NULL") {
      cat("\n", paste("  RS", Z, " ~ MRS", Z, "*1", sep=""))
    } # end (if Z)
    if (W != "NULL") {
      cat("\n", paste("  RS", W, " ~ MRS", W, "*1", sep=""))
    } # (if W)

    # -- Estimate Variance and Covariance of Random Intercepts and Random Slopes -- #
    cat(rep("\n",2), "  # -- Estimate variance and covariance of random intercepts and random slopes -- #")
    cat("\n", "   RI", X, " ~~ RI", X, sep="")
    cat("\n", "   RS", X, " ~~ RS", X, sep="")
    cat("\n", "   RI", Y, " ~~ RI", Y, sep="")
    cat("\n", "   RS", Y, " ~~ RS", Y, sep="")
    if (Z != "NULL") {
      cat("\n", "   RI", Z, " ~~ RI", Z, sep="")
      cat("\n", "   RS", Z, " ~~ RS", Z, sep="")
    } # end (if Z != "NULL")
    if (W != "NULL") {
      cat("\n", "   RI", W, " ~~ RI", W, sep="")
      cat("\n", "   RS", W, " ~~ RS", W, sep="")
    } # end (if W != "NULL")

    cat("\n", "   RI", X, " ~~ RS", X, sep="")
    cat("\n", "   RI", X, " ~~ RI", Y, sep="")
    cat("\n", "   RI", X, " ~~ RS", Y, sep="")
    cat("\n", "   RS", X, " ~~ RI", Y, sep="")
    cat("\n", "   RS", X, " ~~ RS", Y, sep="")
    cat("\n", "   RI", Y, " ~~ RS", Y, sep="")

    if (Z != "NULL") {
      cat("\n", "   RI", X, " ~~ RI", Z, sep="")
      cat("\n", "   RI", X, " ~~ RS", Z, sep="")
      cat("\n", "   RS", X, " ~~ RI", Z, sep="")
      cat("\n", "   RS", X, " ~~ RS", Z, sep="")
      cat("\n", "   RI", Y, " ~~ RI", Z, sep="")
      cat("\n", "   RI", Y, " ~~ RS", Z, sep="")
      cat("\n", "   RS", Y, " ~~ RI", Z, sep="")
      cat("\n", "   RS", Y, " ~~ RS", Z, sep="")
      cat("\n", "   RI", Z, " ~~ RS", Z, sep="")
    } # end (if Z != "NULL")
    if (W != "NULL") {
      cat("\n", "   RI", X, " ~~ RI", W, sep="")
      cat("\n", "   RI", X, " ~~ RS", W, sep="")
      cat("\n", "   RS", X, " ~~ RI", W, sep="")
      cat("\n", "   RS", X, " ~~ RS", W, sep="")
      cat("\n", "   RI", Y, " ~~ RI", W, sep="")
      cat("\n", "   RI", Y, " ~~ RS", W, sep="")
      cat("\n", "   RS", Y, " ~~ RI", W, sep="")
      cat("\n", "   RS", Y, " ~~ RS", W, sep="")
      cat("\n", "   RI", Z, " ~~ RI", W, sep="")
      cat("\n", "   RI", Z, " ~~ RS", W, sep="")
      cat("\n", "   RS", Z, " ~~ RI", W, sep="")
      cat("\n", "   RS", Z, " ~~ RS", W, sep="")
      cat("\n", "   RI", W, " ~~ RS", W, sep="")
    } # end (if W != "NULL")

    # -- Create Latent Variables from Indicators -- #
    cat(rep("\n",2), "  # -- Create latent variables -- #")
    for (i in 1:no.waves) {
      cat("\n", paste("  w", X, i, " =~ 1*", X, i, sep=""))
      cat("\n", paste("  w", Y, i, " =~ 1*", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  w", Z, i, " =~ 1*", Z, i, sep=""))
      }  # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  w", W, i, " =~ 1*", W, i, sep=""))
      }  # end (if W)
    } # end (for i)

    # -- Constrain Means (Intercepts) of Latent Variables to Zero -- #
    cat(rep("\n",2), "  # -- Constrain means (intercepts) of latent variables to zero -- #")
    for (i in 1:no.waves) {
      cat("\n", paste("  w", X, i, " ~ 0*1", sep=""))
      cat("\n", paste("  w", Y, i, " ~ 0*1", sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  w", Z, i, " ~ 0*1", sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  w", W, i, " ~ 0*1", sep=""))
      } # (if W)
    } # end (for i)

    # -- Estimate (Residual) Variance of Latent Variables -- #
    cat(rep("\n",2), "  # -- Estimate (residual) variance of latent variables -- #")
    for (i in 1:no.waves) {
      cat("\n", paste("  w", X, i, " ~~ eXX", i, "*w", X, i, sep=""))
      cat("\n", paste("  w", Y, i, " ~~ eYY", i, "*w", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  w", Z, i, " ~~ eZZ", i, "*w", Z, i, sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  w", W, i, " ~~ eWW", i, "*w", W, i, sep=""))
      } # end (if W)
    } # end (for i)

    # -- Estimate Covariance Between Latent Variables at First Wave -- #
    cat(rep("\n",2), "  # -- Estimate covariance between latent variables at first wave -- #")
    cat("\n", "    w", X, "1 ~~ w", Y, "1", sep="")
    if (Z != "NULL") {
      cat("\n", "    w", X, "1 ~~ w", Z, "1", sep="")
      cat("\n", "    w", Y, "1 ~~ w", Z, "1", sep="")
    }  # end (if Z)
    if (W != "NULL") {
      cat("\n", "    w", X, "1 ~~ w", W, "1", sep="")
      cat("\n", "    w", Y, "1 ~~ w", W, "1", sep="")
      cat("\n", "    w", Z, "1 ~~ w", W, "1", sep="")
    }  # end (if W)

    # -- Estimate Covariances Between Residuals of Latent Variables -- #
    cat(rep("\n",2), "  # -- Estimate covariances between residuals of latent variables -- #")
    cat("\n", "  ###################################################################")
    cat("\n", "  # Remove the subscripts for eXY for invariant residual covariance #")
    cat("\n", "  ###################################################################")
    for (i in 2:no.waves) {
      cat("\n", paste("  w", X, i, " ~~ eXY", i, "*w", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  w", X, i, " ~~ eXZ", i, "*w", Z, i, sep=""))
        cat("\n", paste("  w", Y, i, " ~~ eYZ", i, "*w", Z, i, sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  w", X, i, " ~~ eXW", i, "*w", W, i, sep=""))
        cat("\n", paste("  w", Y, i, " ~~ eYW", i, "*w", W, i, sep=""))
        cat("\n", paste("  w", Z, i, " ~~ eZW", i, "*w", W, i, sep=""))
      } # end (if W)
    } # end ((for i)

    # -- Constrain covariance among RIX RIY, RIZ, RIW and wx1, wy1, wz1, ww1 -- #
    cat(rep("\n",2), "  # -- Constrain covariance among RIx, RIy, RIz, RIw and wx1, wy1, wz1, ww1 -- #")
    cat("\n", "   RI", X, " ~~ 0*w", X, "1", sep="")
    cat("\n", "   RI", X, " ~~ 0*w", Y, "1", sep="")
    cat("\n", "   RI", Y, " ~~ 0*w", X, "1", sep="")
    cat("\n", "   RI", Y, " ~~ 0*w", Y, "1", sep="")
    if (Z != "NULL") {
      cat("\n", "   RI", X, " ~~ 0*w", Z, "1", sep="")
      cat("\n", "   RI", Y, " ~~ 0*w", Z, "1", sep="")
      cat("\n", "   RI", Z, " ~~ 0*w", X, "1", sep="")
      cat("\n", "   RI", Z, " ~~ 0*w", Y, "1", sep="")
      cat("\n", "   RI", Z, " ~~ 0*w", Z, "1", sep="")
    } # end (if Z)
    if (W != "NULL") {
      cat("\n", "   RI", X, " ~~ 0*w", W, "1", sep="")
      cat("\n", "   RI", Y, " ~~ 0*w", W, "1", sep="")
      cat("\n", "   RI", Z, " ~~ 0*w", W, "1", sep="")
      cat("\n", "   RI", W, " ~~ 0*w", X, "1", sep="")
      cat("\n", "   RI", W, " ~~ 0*w", Y, "1", sep="")
      cat("\n", "   RI", W, " ~~ 0*w", Z, "1", sep="")
      cat("\n", "   RI", W, " ~~ 0*w", W, "1", sep="")
    } # end (if W)

    # -- Constrain covariance among RSX RSY, RSZ, RSW and wx1, wy1, wz1, ww1 -- #
    cat(rep("\n",2), "  # -- Constrain covariance among RSx, RSy, RSz, RSw and wx1, wy1, wz1, ww1 -- #")
    cat("\n", "   RS", X, " ~~ 0*w", X, "1", sep="")
    cat("\n", "   RS", X, " ~~ 0*w", Y, "1", sep="")
    cat("\n", "   RS", Y, " ~~ 0*w", X, "1", sep="")
    cat("\n", "   RS", Y, " ~~ 0*w", Y, "1", sep="")
    if (Z != "NULL") {
      cat("\n", "   RS", X, " ~~ 0*w", Z, "1", sep="")
      cat("\n", "   RS", Y, " ~~ 0*w", Z, "1", sep="")
      cat("\n", "   RS", Z, " ~~ 0*w", X, "1", sep="")
      cat("\n", "   RS", Z, " ~~ 0*w", Y, "1", sep="")
      cat("\n", "   RS", Z, " ~~ 0*w", Z, "1", sep="")
    } # end (if Z != "NULL")
    if (W != "NULL") {
      cat("\n", "   RS", X, " ~~ 0*w", W, "1", sep="")
      cat("\n", "   RS", Y, " ~~ 0*w", W, "1", sep="")
      cat("\n", "   RS", Z, " ~~ 0*w", W, "1", sep="")
      cat("\n", "   RS", W, " ~~ 0*w", X, "1", sep="")
      cat("\n", "   RS", W, " ~~ 0*w", Y, "1", sep="")
      cat("\n", "   RS", W, " ~~ 0*w", Z, "1", sep="")
      cat("\n", "   RS", W, " ~~ 0*w", W, "1", sep="")
    } # end (if W != "NULL")

    for (j in 1:lag) {
      cat(rep("\n",2), paste0("  # -- Estimate lagged effects between latent variables (Lag = ", j, " wave) -- #", sep =""))
      cat("\n", "  #############################################")
      cat("\n", "  # Remove the subscripts for invariant paths #")
      cat("\n", "  #############################################")
      for (i in (j+1):no.waves) {
        if (Z == "NULL") {
          cat("\n", paste("  w", X,i, " ~ pXX", i,i-j, "*w", X,i-j, " + pYX", i,i-j, "*w", Y,i-j, sep=""))
          cat("\n", paste("  w", Y,i, " ~ pXY", i,i-j, "*w", X,i-j, " + pYY", i,i-j, "*w", Y,i-j, sep=""))
        } else if (W != "NULL") {
          cat("\n", paste("  w", X,i, " ~ pXX", i,i-j, "*w", X,i-j, " + pYX", i,i-j, "*w", Y,i-j, " + pZX", i,i-j, "*w", Z,i-j, " + pWX", i,i-j, "*w", W,i-j, sep=""))
          cat("\n", paste("  w", Y,i, " ~ pXY", i,i-j, "*w", X,i-j, " + pYY", i,i-j, "*w", Y,i-j, " + pZY", i,i-j, "*w", Z,i-j, " + pWY", i,i-j, "*w", W,i-j, sep=""))
          cat("\n", paste("  w", Z,i, " ~ pXZ", i,i-j, "*w", X,i-j, " + pYZ", i,i-j, "*w", Y,i-j, " + pZZ", i,i-j, "*w", Z,i-j, " + pWZ", i,i-j, "*w", W,i-j, sep=""))
          cat("\n", paste("  w", W,i, " ~ pXW", i,i-j, "*w", X,i-j, " + pYW", i,i-j, "*w", Y,i-j, " + pZW", i,i-j, "*w", Z,i-j, " + pWW", i,i-j, "*w", W,i-j, sep=""))
        } else if (Z != "NULL") {
          cat("\n", paste("  w", X,i, " ~ pXX", i,i-j, "*w", X,i-j, " + pYX", i,i-j, "*w", Y,i-j, " + pZX", i,i-j, "*w", Z,i-j, sep=""))
          cat("\n", paste("  w", Y,i, " ~ pXY", i,i-j, "*w", X,i-j, " + pYY", i,i-j, "*w", Y,i-j, " + pZY", i,i-j, "*w", Z,i-j, sep=""))
          cat("\n", paste("  w", Z,i, " ~ pXZ", i,i-j, "*w", X,i-j, " + pYZ", i,i-j, "*w", Y,i-j, " + pZZ", i,i-j, "*w", Z,i-j, sep=""))
        } # end (if Z)
      } # end (for i)
    } # end (for j)

    cat(rep("\n",2), "  ##########################################")
    cat("\n", "  # Regression of indicators on C1 #")
    cat("\n", "  ##########################################")
    for (i in 1:no.waves) {
      cat("\n", paste("  #  ", X, i, " ~ sx", i,"*C1", sep=""))
      cat("\n", paste("  #  ", Y, i, " ~ sy", i,"*C1", sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  #  ", Z, i, " ~ sz", i,"*C1", sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  #  ", W, i, " ~ sz", i,"*C1", sep=""))
      } # end (if W)
    } # end (for i)

    cat(rep("\n",2), "  #########################################")
    cat("\n", "  # Regression of random intercepts on C1 #")
    cat("\n", "  #########################################")
    if (Z == "NULL") {
      cat("\n", "   #  RI", X, " + RI", Y, " ~ C1", sep="")
    } else if (Z != "NULL") {
      cat("\n", "   #  RI", X, " + RI", Y, " + RI", Z, " ~ C1", sep="")
    } else if (W != "NULL") {
      cat("\n", "   #  RI", X, " + RI", Y, " + RI", Z, " + RI", W, " ~ C1", sep="")
    } # end (if Z/W)

    cat(rep("\n",2), "  ################################################################")
    cat("\n", "  # Regression of time-invariant outcome D1 on random intercepts #")
    cat("\n", "  ################################################################")
    if (Z == "NULL") {
      cat("\n", "   #  D1 ~ RI", X, " + RI", Y, sep="")
    } else if (Z != "NULL") {
      cat("\n", "   #  D1 ~ RI", X, " + RI", Y, " + RI", Z, sep="")
    } else if (W != "NULL") {
      cat("\n", "   #  D1 ~ RI", X, " + RI", Y, " + RI", Z, " + RI", W, sep="")
    } # end (if Z/W)
    cat("\n", "  #  D1 ~~ D1 # Residual variance of D1")

    cat(rep("\n",2), "  ################################################")
    cat("\n", "  # Regression of outcome D1 on latent variables #")
    cat("\n", "  ################################################")
    for (i in 1:no.waves) {
      cat("\n", paste("  #  D1 ~ bx", i, "*w", X, i, sep=""))
      cat("\n", paste("  #  D1 ~ by", i, "*w", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  #  D1 ~ bz", i, "*w", Z, i, sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  #  D1 ~ bw", i, "*w", W, i, sep=""))
      } # end (if W)
    } # end (for i)
    cat("\n", "  #  D1 ~~ D1 # Residual variance of D1")

    cat(rep("\n",2), "  '")

    # -- Run LGCMSRMLR -- #
    cat(rep("\n",2), "  LGCMSRMLR.fit <- lavaan::sem(LGCMSR, ")
    cat("\n", "   ", data.source, ",")
    cat("\n", "   missing = 'fiml',")
    cat("\n", "   meanstructure = TRUE,")
    cat("\n", "   information = 'observed',")
    cat("\n", "   estimator = 'MLR'")
    cat("\n", ")")

    # Request summary outputs
    cat(rep("\n",2), "  print(lavaan::summary(LGCMSRMLR.fit, fit.measure = TRUE, standardized = TRUE, rsq = TRUE))", rep("\n", 3))

  sink() # Stop writing to file


  ## -- Execute LGCMSR.txt and request summary outputs-- ##
  source('LGCMSR.txt')
  if (lavaan::lavInspect(LGCMSRMLR.fit, what ="post.check") == FALSE) stop("The lavaan solution is non-admissible.")
  ## ------------------------------- ##


  ## -- Monte Carlo Simulation -- ##

  parEst <- lavaan::parameterEstimates(LGCMSRMLR.fit, remove.nonfree = TRUE)
  pest2 <- parEst[, "est"]  # Estimated Parameters
  pest3 <- lavaan::lavTech(LGCMSRMLR.fit, what = "vcov", add.labels = TRUE)  # Estimated Variance-Covariance of Estimated Parameters

  Invariance(parEst, pest2, pest3, no.path, no.waves, lag, Type1, Type1Adj, X, Y, Z, W)

  cat(rep("\n", 2))

}  # end (Function LGCMSR)

## ========================================================================================== ##





# ==================== Creating Function "STARTS" ==================== #
#' Function STARTS (Stable Trait Autoregressive Trait and State Model)
#'
#' Stable Trait Autoregressive Trait and State Model (STARTS)
#'
#' @param data.source name of data.frame.
#' @param no.waves number of waves (minimum = 3, must be grater than lag).
#' @param lag number of waves between two lags (minimum = 1, maximum = 4).
#' @param varI.eq indicator residual variances are invariant across waves if TRUE (default is FALSE).
#' @param Type1 Overall Type I error rate (default is 0.05) for comparing estimated parameters in the List and Delete method.
#' @param Type1Adj Adjustment of Type I error rate for multiple tests for each estimated parameter. Default is "BON" (Bonferroni adjustment -- Type I error rate/no. of pairwise comparisons), can also be "NULL" (without adjustment).
#' @param X name of variable X.
#' @param Y name of variable Y.
#' @param Z name of variable Z.
#' @param W name of variable W (Z must come before W).
#'
#' @return Stable Trait Autoregressive Trait and State Model (STARTS) outputs.
#' @export
#' @examples
#'
#' ## -- Example -- ##
#'
#' STARTS(data.source="Data_A", 8, 1, varI.eq=TRUE, Type1=0.05, Type1Adj="BON", X="EXPOSE.",
#'  Y="INTENS.")
#'

STARTS <- function(data.source, no.waves, lag=1, varI.eq = FALSE, Type1=0.05, Type1Adj="BON", X, Y, Z="NULL", W = "NULL") {

  ## -- Check inputs -- ##

  if (no.waves < 3) stop("Minimum number of waves is 3")

  if (lag < 1) stop("Minimum number of lag is 1")
  if (lag > 4) stop("Maximum number of lag is 4")

  if ((no.waves - lag) < 1) stop("Number of waves must be greater than lag plus 1")

  if (Type1 > 0.05) stop("Type I error rate > 0.05 is not recommended")
  if (Type1 < 0.0001) stop("Type I error rate < 0.0001 is not recommended")

  if (Z == "NULL" & W != "NULL") stop("Z must be defined before W")

  if (is.logical(varI.eq) == FALSE) stop("varI.eq (equivalence of indicator residual variance) can only be TRUE or FALSE")

  Type1Adj <- toupper(Type1Adj)
  match.arg(Type1Adj, c("BON", "NULL"))

  ## ----- Creating Model STARTS ----- ###
  sink('STARTS.txt')
    cat("\n", "# Specify the model (STARTS)", "\n")
    cat("\n", "STARTS <- '")

    # -- Create Between Components (Random Intercepts) -- #
    cat(rep("\n",2), "  # -- Create between components (random intercepts) -- #")
    BX <- paste("  RI", X, " =~ 1*", X, "1", sep="")
    BY <- paste("  RI", Y, " =~ 1*", Y, "1", sep="")
    for (i in 2:no.waves) {
      BX <- paste(BX, " +1*", X, i, sep="")
      BY <- paste(BY, " +1*", Y, i, sep="")
    } # end (for i)
    cat("\n", BX)
    cat("\n", BY)
    if (Z != "NULL") {
      BZ <- paste("  RI", Z, " =~ 1*", Z, "1", sep="")
      for (i in 2:no.waves) {
        BZ <- paste(BZ, " +1*", Z, i, sep="")
      } # end (for i)
      cat("\n", BZ)
    } # end (if Z)
    if (W != "NULL") {
      BW <- paste("  RI", W, " =~ 1*", W, "1", sep="")
      for (i in 2:no.waves) {
        BW <- paste(BW, " +1*", W, i, sep="")
      } # end (for i)
      cat("\n", BW)
    } # end (if W)

    # -- Create Residual Variance of Indicators -- #
    if (isFALSE(varI.eq)) {
      cat(rep("\n",2), "  # -- Create residual variance of indicators -- #")
      cat("\n", "  ##########################################################")
      cat("\n", "  # Remove the subscripts for invariant indicator variance #")
      cat("\n", "  ##########################################################")
      for (i in 1:no.waves) {
        cat("\n", paste("  ", X, i, " ~~ eIXX", i, "*", X, i, sep=""))
        cat("\n", paste("  ", Y, i, " ~~ eIYY", i, "*", Y, i, sep=""))
        if (Z != "NULL") {
          cat("\n", paste("  ", Z, i, " ~~ eIZZ", i, "*", Z, i, sep=""))
        } # end (if Z)
        if (W != "NULL") {
          cat("\n", paste("  ", W, i, " ~~ eIWW", i, "*", W, i, sep=""))
        } # end (if W)
      } # end (for i)###   ###
    } else {
      cat(rep("\n",2), "  # -- Create residual variance of indicators -- #")
      cat("\n", "  ##########################################################")
      cat("\n", "  # Remove the subscripts for invariant indicator variance #")
      cat("\n", "  ##########################################################")
      for (i in 1:no.waves) {
        cat("\n", paste("  ", X, i, " ~~ eIXX*", X, i, sep=""))
        cat("\n", paste("  ", Y, i, " ~~ eIYY*", Y, i, sep=""))
        if (Z != "NULL") {
          cat("\n", paste("  ", Z, i, " ~~ eIZZ*", Z, i, sep=""))
        } # end (if Z)
        if (W != "NULL") {
          cat("\n", paste("  ", W, i, " ~~ eIWW*", W, i, sep=""))
        } # end (if W)
      } # end (for i)###   ###
    } # end (if varI.eq == FALSE)


    # -- Create Covariance of Indicator Residuals -- #
    cat(rep("\n",2), "  # -- Create covariance of indicator residuals -- #")
    cat("\n", "  ############################################################")
    cat("\n", "  # Remove the subscripts for invariant indicator covariance #")
    cat("\n", "  ############################################################")
    for (i in 1:no.waves) {
      cat("\n", paste("  ", X, i, " ~~ eIXY", i, "*", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  ", X, i, " ~~ eIXZ", i, "*", Z, i, sep=""))
        cat("\n", paste("  ", Y, i, " ~~ eIYZ", i, "*", Z, i, sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  ", X, i, " ~~ eIXW", i, "*", W, i, sep=""))
        cat("\n", paste("  ", Y, i, " ~~ eIYW", i, "*", W, i, sep=""))
        cat("\n", paste("  ", Z, i, " ~~ eIZW", i, "*", W, i, sep=""))
      } # end (if W)
    } # end (for i)

    # -- Create Latent Variables from Indicators -- #
    cat(rep("\n",2), "  # -- Create latent variables -- #")
    for (i in 1:no.waves) {
      cat("\n", paste("  w", X, i, " =~ 1*", X, i, sep=""))
      cat("\n", paste("  w", Y, i, " =~ 1*", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  w", Z, i, " =~ 1*", Z, i, sep=""))
      }  # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  w", W, i, " =~ 1*", W, i, sep=""))
      }  # end (if W)
    } # end (for i)

    cat(rep("\n",2), "  # -- Estimate covariances between residuals of latent variables -- #")
    cat("\n", "  ###################################################################")
    cat("\n", "  # Remove the subscripts for eXY for invariant residual covariance #")
    cat("\n", "  ###################################################################")
    for (i in 2:no.waves) {
      cat("\n", paste("  w", X, i, " ~~ eXY", i, "*w", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  w", X, i, " ~~ eXZ", i, "*w", Z, i, sep=""))
        cat("\n", paste("  w", Y, i, " ~~ eYZ", i, "*w", Z, i, sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  w", X, i, " ~~ eXW", i, "*w", W, i, sep=""))
        cat("\n", paste("  w", Y, i, " ~~ eYW", i, "*w", W, i, sep=""))
        cat("\n", paste("  w", Z, i, " ~~ eZW", i, "*w", W, i, sep=""))
      } # end (if W)
    } # end ((for i)


    cat(rep("\n",2), "  # -- Estimate covariance between latent variables at first wave -- #")
    cat("\n", "   w", X, "1 ~~ w", Y, "1", sep="")
    if (Z != "NULL") {
      cat("\n", "    w", X, "1 ~~ w", Z, "1", sep="")
      cat("\n", "    w", Y, "1 ~~ w", Z, "1", sep="")
    } # end (if Z)
    if (W != "NULL") {
      cat("\n", "    w", X, "1 ~~ w", W, "1", sep="")
      cat("\n", "    w", Y, "1 ~~ w", W, "1", sep="")
      cat("\n", "    w", Z, "1 ~~ w", W, "1", sep="")
    } # end (if W)

    # -- Estimate Variance and Covariance of Random Intercepts -- #
    cat(rep("\n",2), "  # -- Estimate variance and covariance of random intercepts -- #")
    cat("\n", "   RI", X, " ~~ RI", X, sep="")
    cat("\n", "   RI", Y, " ~~ RI", Y, sep="")
    cat("\n", "   RI", X, " ~~ RI", Y, sep="")
    if (Z != "NULL") {
      cat("\n", "   RI", Z, " ~~ RI", Z, sep="")
      cat("\n", "   RI", X, " ~~ RI", Z, sep="")
      cat("\n", "   RI", Y, " ~~ RI", Z, sep="")
    } # end (if Z != "NULL")
    if (W != "NULL") {
      cat("\n", "   RI", W, " ~~ RI", W, sep="")
      cat("\n", "   RI", X, " ~~ RI", W, sep="")
      cat("\n", "   RI", Y, " ~~ RI", W, sep="")
      cat("\n", "   RI", Z, " ~~ RI", W, sep="")
    } # end (if W != "NULL")


    # -- Estimate (Residual) Variance of Latent Variables -- #
    cat(rep("\n",2), "  # -- Estimate (residual) variance of latent variables -- #")
    for (i in 1:no.waves) {
      cat("\n", paste("  w", X, i, " ~~ eXX", i, "*w", X, i, sep=""))
      cat("\n", paste("  w", Y, i, " ~~ eYY", i,  "*w", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  w", Z, i, " ~~ eZZ", i, "*w", Z, i, sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  w", W, i, " ~~ eWW", i, "*w", W, i, sep=""))
      } # end (if W)
    } # end (for i)


    # -- Estimate Grand Means (Intercepts) of Indicators -- #
    cat(rep("\n",2), "  # -- Estimate grand means (intercepts) of indicators -- #")
    for (i in 1:no.waves) {
      cat("\n", paste("  ", X, i, " ~ Mw", X, i, "*1", sep=""))
      cat("\n", paste("  ", Y, i, " ~ Mw", Y, i, "*1", sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  ", Z, i, " ~ Mw", Z, i, "*1", sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  ", W, i, " ~ Mw", W, i, "*1", sep=""))
      } # (if W)
    } # end (for i)


    # -- Constrain covariance among RIX RIY, RIZ, RIW and wx1, wy1, wz1, ww1 -- #
    cat(rep("\n",2), "  # -- Constrain covariance among RIx, RIy, RIz, RIw and wx1, wy1, wz1, ww1 -- #")
    cat("\n", "   RI", X, " ~~ 0*w", X, "1", sep="")
    cat("\n", "   RI", X, " ~~ 0*w", Y, "1", sep="")
    cat("\n", "   RI", Y, " ~~ 0*w", X, "1", sep="")
    cat("\n", "   RI", Y, " ~~ 0*w", Y, "1", sep="")
    if (Z != "NULL") {
      cat("\n", "   RI", X, " ~~ 0*w", Z, "1", sep="")
      cat("\n", "   RI", Y, " ~~ 0*w", Z, "1", sep="")
      cat("\n", "   RI", Z, " ~~ 0*w", X, "1", sep="")
      cat("\n", "   RI", Z, " ~~ 0*w", Y, "1", sep="")
      cat("\n", "   RI", Z, " ~~ 0*w", Z, "1", sep="")
    } # end (if Z)
    if (W != "NULL") {
      cat("\n", "   RI", X, " ~~ 0*w", W, "1", sep="")
      cat("\n", "   RI", Y, " ~~ 0*w", W, "1", sep="")
      cat("\n", "   RI", Z, " ~~ 0*w", W, "1", sep="")
      cat("\n", "   RI", W, " ~~ 0*w", X, "1", sep="")
      cat("\n", "   RI", W, " ~~ 0*w", Y, "1", sep="")
      cat("\n", "   RI", W, " ~~ 0*w", Z, "1", sep="")
      cat("\n", "   RI", W, " ~~ 0*w", W, "1", sep="")
    } # end (if W)

    for (j in 1:lag) {
      cat(rep("\n",2), paste0("  # -- Estimate lagged effects between latent variables (Lag = ", j, " wave) -- #", sep =""))
      cat("\n", "  #############################################")
      cat("\n", "  # Remove the subscripts for invariant paths #")
      cat("\n", "  #############################################")
      for (i in (j+1):no.waves) {
        if (Z == "NULL") {
          cat("\n", paste("  w", X,i, " ~ pXX", i,i-j, "*w", X,i-j, " + pYX", i,i-j, "*w", Y,i-j, sep=""))
          cat("\n", paste("  w", Y,i, " ~ pXY", i,i-j, "*w", X,i-j, " + pYY", i,i-j, "*w", Y,i-j, sep=""))
        } else if (W != "NULL") {
          cat("\n", paste("  w", X,i, " ~ pXX", i,i-j, "*w", X,i-j, " + pYX", i,i-j, "*w", Y,i-j, " + pZX", i,i-j, "*w", Z,i-j, " + pWX", i,i-j, "*w", W,i-j, sep=""))
          cat("\n", paste("  w", Y,i, " ~ pXY", i,i-j, "*w", X,i-j, " + pYY", i,i-j, "*w", Y,i-j, " + pZY", i,i-j, "*w", Z,i-j, " + pWY", i,i-j, "*w", W,i-j, sep=""))
          cat("\n", paste("  w", Z,i, " ~ pXZ", i,i-j, "*w", X,i-j, " + pYZ", i,i-j, "*w", Y,i-j, " + pZZ", i,i-j, "*w", Z,i-j, " + pWZ", i,i-j, "*w", W,i-j, sep=""))
          cat("\n", paste("  w", W,i, " ~ pXW", i,i-j, "*w", X,i-j, " + pYW", i,i-j, "*w", Y,i-j, " + pZW", i,i-j, "*w", Z,i-j, " + pWW", i,i-j, "*w", W,i-j, sep=""))
        } else if (Z != "NULL") {
          cat("\n", paste("  w", X,i, " ~ pXX", i,i-j, "*w", X,i-j, " + pYX", i,i-j, "*w", Y,i-j, " + pZX", i,i-j, "*w", Z,i-j, sep=""))
          cat("\n", paste("  w", Y,i, " ~ pXY", i,i-j, "*w", X,i-j, " + pYY", i,i-j, "*w", Y,i-j, " + pZY", i,i-j, "*w", Z,i-j, sep=""))
          cat("\n", paste("  w", Z,i, " ~ pXZ", i,i-j, "*w", X,i-j, " + pYZ", i,i-j, "*w", Y,i-j, " + pZZ", i,i-j, "*w", Z,i-j, sep=""))
        } # end (if Z)
      } # end (for i)
    } # end (for j)

    cat(rep("\n",2), "  ##########################################")
    cat("\n", "  # Regression of indicators on C1 #")
    cat("\n", "  ##########################################")
    for (i in 1:no.waves) {
      cat("\n", paste("  #  ", X, i, " ~ sx", i,"*C1", sep=""))
      cat("\n", paste("  #  ", Y, i, " ~ sy", i,"*C1", sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  #  ", Z, i, " ~ sz", i,"*C1", sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  #  ", W, i, " ~ sz", i,"*C1", sep=""))
      } # end (if W)
    } # end (for i)

    cat(rep("\n",2), "  #########################################")
    cat("\n", "  # Regression of random intercepts on C1 #")
    cat("\n", "  #########################################")
    if (Z == "NULL") {
      cat("\n", "   #  RI", X, " + RI", Y, " ~ C1", sep="")
    } else if (Z != "NULL") {
      cat("\n", "   #  RI", X, " + RI", Y, " + RI", Z, " ~ C1", sep="")
    } else if (W != "NULL") {
      cat("\n", "   #  RI", X, " + RI", Y, " + RI", Z, " + RI", W, " ~ C1", sep="")
    } # end (if Z/W)

    cat(rep("\n",2), "  ################################################################")
    cat("\n", "  # Regression of time-invariant outcome D1 on random intercepts #")
    cat("\n", "  ################################################################")
    if (Z == "NULL") {
      cat("\n", "   #  D1 ~ RI", X, " + RI", Y, sep="")
    } else if (Z != "NULL") {Type1Adj
      cat("\n", "   #  D1 ~ RI", X, " + RI", Y, " + RI", Z, sep="")
    } else if (W != "NULL") {
      cat("\n", "   #  D1 ~ RI", X, " + RI", Y, " + RI", Z, " + RI", W, sep="")
    } # end (if Z/W)
    cat("\n", "  #  D1 ~~ D1 # Residual variance of D1")

    cat(rep("\n",2), "  ################################################")
    cat("\n", "  # Regression of outcome D1 on latent variables #")
    cat("\n", "  ################################################")
    for (i in 1:no.waves) {
      cat("\n", paste("  #  D1 ~ bx", i, "*w", X, i, sep=""))
      cat("\n", paste("  #  D1 ~ by", i, "*w", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  #  D1 ~ bz", i, "*w", Z, i, sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  #  D1 ~ bw", i, "*w", W, i, sep=""))
      } # end (if W)
    } # end (for i)
    cat("\n", "  #  D1 ~~ D1 # Residual variance of D1")

    cat(rep("\n",2), "  '")

    # -- Run STARTSMLR -- #
    cat(rep("\n",2), "  STARTSMLR.fit <- lavaan::sem(STARTS, ")
    cat("\n", "   ", data.source, ",")
    cat("\n", "   missing = 'fiml',")
    cat("\n", "   meanstructure = TRUE,")
    cat("\n", "   information = 'observed',")
    cat("\n", "   estimator = 'MLR'")
    cat("\n", ")")

    # Request summary outputs
    cat(rep("\n",2), "  print(lavaan::summary(STARTSMLR.fit, fit.measure = TRUE, standardized = TRUE, rsq = TRUE))", rep("\n",3))

  sink() # Stop writing to file


  ## -- Execute STARTS.txt and request summary outputs-- ##
  source('STARTS.txt')
  if (lavaan::lavInspect(STARTSMLR.fit, what ="post.check") == FALSE) stop("The lavaan solution is non-admissible.")
  ## ------------------------------- ##


  ## -- Monte Carlo Simulation -- ##

  parEst <- lavaan::parameterEstimates(STARTSMLR.fit, remove.nonfree = TRUE)
  pest2 <- parEst[, "est"]  # Estimated Parameters
  pest3 <- lavaan::lavTech(STARTSMLR.fit, what = "vcov", add.labels = TRUE)  # Estimated Variance-Covariance of Estimated Parameters

  Invariance(parEst, pest2, pest3, no.path, no.waves, lag, Type1, Type1Adj, X, Y, Z, W)

  cat(rep("\n", 2))


}  # end (Function STARTS)

## ========================================================================================== ##





# ==================== Creating Function "ALT" ==================== #
#' Function ALT (Autoregressive latent trajectory model)
#'
#' Autoregressive Latent Trajectory Model (ALT)
#'
#' @param data.source name of data.frame
#' @param no.waves number of waves (minimum = 3, must be grater than lag)
#' @param lag number of waves between two lags (minimum = 1, maximum = 4)
#' @param Type1 Overall Type I error rate (default is 0.05) for comparing estimated parameters in the List and Delete method.
#' @param Type1Adj Adjustment of Type I error rate for multiple tests for each estimated parameter. Default is "BON" (Bonferroni adjustment -- Type I error rate/no. of pairwise comparisons), can also be "NULL" (without adjustment).
#' @param X name of variable X.
#' @param Y name of variable Y.
#' @param Z name of variable Z.
#' @param W name of variable W (Z must come before W).
#'
#' @return Autoregressive Latent Trajectory Model (ALT) outputs.
#' @export
#' @examples
#'
#' ## -- Example -- ##
#'
#' ALT(data.source="Data_A", 7, 1, Type1=0.05, Type1Adj="BON", X="EXPOSE.", Y="INTENS.")
#'

ALT <- function(data.source, no.waves, lag=1, Type1=0.05, Type1Adj="BON", X, Y, Z="NULL", W = "NULL") {

  ## -- Check inputs -- ##

  if (no.waves < 3) stop("Minimum number of waves is 3")

  if (lag < 1) stop("Minimum number of lag is 1")
  if (lag > 4) stop("Maximum number of lag is 4")

  if ((no.waves - lag) < 1) stop("Number of waves must be greater than lag plus 1")

  if (Type1 > 0.05) stop("Type I error rate > 0.05 is not recommended")
  if (Type1 < 0.0001) stop("Type I error rate < 0.0001 is not recommended")

  if (Z == "NULL" & W != "NULL") stop("Z must be defined before W")

  Type1Adj <- toupper(Type1Adj)
  match.arg(Type1Adj, c("BON", "NULL"))

  ## ----- Creating Model ALT ----- ###
  sink('ALT.txt')
    cat("\n", "# Specify the model (ALT)", "\n")
    cat("\n", "ALT <- '")

    # -- Create Latent Variables from Indicators -- #
    cat(rep("\n",2), "  # -- Create latent variables -- #")
    for (i in 1:no.waves) {
      cat("\n", paste("  w", X, i, " =~ 1*", X, i, sep=""))
      cat("\n", paste("  w", Y, i, " =~ 1*", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  w", Z, i, " =~ 1*", Z, i, sep=""))
      }  # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  w", W, i, " =~ 1*", W, i, sep=""))
      }  # end (if W)
    } # end (for i)

    # -- Constrain Residual Variance of Indicators to Zero -- #
    cat(rep("\n",2), "  # -- Constrain residual variance of indicators to zero -- #")
    for (i in 1:no.waves) {
      cat("\n", paste("  ", X, i, " ~~ 0*", X, i, sep=""))
      cat("\n", paste("  ", Y, i, " ~~ 0*", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  ", Z, i, " ~~ 0*", Z, i, sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  ", W, i, " ~~ 0*", W, i, sep=""))
      } # end (if W)
    } # end (for i)

    # -- Create Between Components (Random Intercepts) from Latent Variables -- #
    cat(rep("\n",2), "  # -- Create between components (random intercepts) -- #")
    BX <- paste("  RI", X, " =~ 1*w", X, "2", sep="")
    BY <- paste("  RI", Y, " =~ 1*w", Y, "2", sep="")
    for (i in 3:no.waves) {
      BX <- paste(BX, " +1*w", X, i, sep="")
      BY <- paste(BY, " +1*w", Y, i, sep="")
    } # end (for i)
    cat("\n", BX)
    cat("\n", BY)
    if (Z != "NULL") {
      BZ <- paste("  RI", Z, " =~ 1*w", Z, "2", sep="")
      for (i in 3:no.waves) {
        BZ <- paste(BZ, " +1*w", Z, i, sep="")
      } # end (for i)
      cat("\n", BZ)
    } # end (if Z)
    if (W != "NULL") {
      BW <- paste("  RI", W, " =~ 1*w", W, "2", sep="")
      for (i in 3:no.waves) {
        BW <- paste(BW, " +1*w", W, i, sep="")
      } # end (for i)
      cat("\n", BW)
    } # end (if W)

    # -- Create Between Components (Random Slopes) from Latent Variables -- #
    cat(rep("\n",2), "  # -- Create between components (random slopes) -- #")
    BX <- paste("  RS", X, " =~ 1*w", X, "2", sep="")
    BY <- paste("  RS", Y, " =~ 1*w", Y, "2", sep="")
    for (i in 3:no.waves) {
      BX <- paste(BX, " +", (i-1), "*w", X, i, sep="")
      BY <- paste(BY, " +", (i-1), "*w", Y, i, sep="")
    } # end (for i)
    cat("\n", BX)
    cat("\n", BY)
    if (Z != "NULL") {
      BZ <- paste("  RS", Z, " =~ 1*w", Z, "2", sep="")
      for (i in 3:no.waves) {
        BZ <- paste(BZ, " +", (i-1), "*w", Z, i, sep="")
      } # end (for i)
      cat("\n", BZ)
    } # end (if Z)
    if (W != "NULL") {
      BW <- paste("  RS", W, " =~ 1*w", W, "2", sep="")
      for (i in 3:no.waves) {
        BW <- paste(BW, " +", (i-1), "*w", W, i, sep="")
      } # end (for i)
      cat("\n", BW)
    } # end (if W)

    # -- Constrain Means (Intercepts) of Indicators to Zero -- #
    cat(rep("\n",2), "  # -- Constrain means (intercepts) of indicators to zero -- #")
    for (i in 1:no.waves) {
      cat("\n", paste("  ", X, i, " ~ 0*1", sep=""))
      cat("\n", paste("  ", Y, i, " ~ 0*1", sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  ", Z, i, " ~ 0*1", sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  ", W, i, " ~ 0*1", sep=""))
      } # (if W)
    } # end (for i)

    # -- Constrain Intercepts of Latent Variables to Zero -- #
    cat(rep("\n",2), "  # -- Constrain intercepts of latent variables to zero -- #")
    for (i in 2:no.waves) {
      cat("\n", paste("  w", X, i, " ~ 0*1", sep=""))
      cat("\n", paste("  w", Y, i, " ~ 0*1", sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  w", Z, i, " ~ 0*1", sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  w", W, i, " ~ 0*1", sep=""))
      } # (if W)
    } # end (for i)

    # -- Estimate Means of First Indicator -- #
    cat(rep("\n",2), "  # -- Estimate means of first indicator -- #")
    cat("\n", paste("  w", X, "1 ~ MwX1*1", sep=""))
    cat("\n", paste("  w", Y, "1 ~ MwY1*1", sep=""))
    if (Z != "NULL") {
      cat("\n", paste("  w", Z, "1 ~ MwZ1*1", sep=""))
    } # end (if Z)
    if (W != "NULL") {
      cat("\n", paste("  w", W, "1 ~ MwW1*1", sep=""))
    } # (if W)

    # -- Estimate Means (Intercepts) of Random Intercepts -- #
    cat(rep("\n",2), "  # -- Estimate means (intercepts) of random intercepts -- #")
    cat("\n", paste("  RI", X, " ~ MRI", X, "*1", sep=""))
    cat("\n", paste("  RI", Y, " ~ MRI", Y, "*1", sep=""))
    if (Z != "NULL") {
      cat("\n", paste("  RI", Z, " ~ MRI", Z, "*1", sep=""))
    } # end (if Z)
    if (W != "NULL") {
      cat("\n", paste("  RI", W, " ~ MRI", W, "*1", sep=""))
    } # (if W)

    # -- Estimate Means (Intercepts) of Random Slopes -- #
    cat(rep("\n",2), "  # -- Estimate means (intercepts) of random slopes -- #")
    cat("\n", paste("  RS", X, " ~ MRS", X, "*1", sep=""))
    cat("\n", paste("  RS", Y, " ~ MRS", Y, "*1", sep=""))
    if (Z != "NULL") {
      cat("\n", paste("  RS", Z, " ~ MRS", Z, "*1", sep=""))
    } # end (if Z)
    if (W != "NULL") {
      cat("\n", paste("  RS", W, " ~ MRS", W, "*1", sep=""))
    } # (if W)

    # -- Estimate Variance and Covariance of Random Intercepts and Random Slopes -- #
    cat(rep("\n",2), "  # -- Estimate variance and covariance of random intercepts and random slopes -- #")
    cat("\n", "   RI", X, " ~~ RI", X, sep="")
    cat("\n", "   RS", X, " ~~ RS", X, sep="")
    cat("\n", "   RI", Y, " ~~ RI", Y, sep="")
    cat("\n", "   RS", Y, " ~~ RS", Y, sep="")
    if (Z != "NULL") {
      cat("\n", "   RI", Z, " ~~ RI", Z, sep="")
      cat("\n", "   RS", Z, " ~~ RS", Z, sep="")
    } # end (if Z != "NULL")
    if (W != "NULL") {
      cat("\n", "   RI", W, " ~~ RI", W, sep="")
      cat("\n", "   RS", W, " ~~ RS", W, sep="")
    } # end (if W != "NULL")

    cat("\n", "   RI", X, " ~~ RS", X, sep="")
    cat("\n", "   RI", X, " ~~ RI", Y, sep="")
    cat("\n", "   RI", X, " ~~ RS", Y, sep="")
    cat("\n", "   RS", X, " ~~ RI", Y, sep="")
    cat("\n", "   RS", X, " ~~ RS", Y, sep="")
    cat("\n", "   RI", Y, " ~~ RS", Y, sep="")

    if (Z != "NULL") {
      cat("\n", "   RI", X, " ~~ RI", Z, sep="")
      cat("\n", "   RI", X, " ~~ RS", Z, sep="")
      cat("\n", "   RS", X, " ~~ RI", Z, sep="")
      cat("\n", "   RS", X, " ~~ RS", Z, sep="")
      cat("\n", "   RI", Y, " ~~ RI", Z, sep="")
      cat("\n", "   RI", Y, " ~~ RS", Z, sep="")
      cat("\n", "   RS", Y, " ~~ RI", Z, sep="")
      cat("\n", "   RS", Y, " ~~ RS", Z, sep="")
      cat("\n", "   RI", Z, " ~~ RS", Z, sep="")
    } # end (if Z != "NULL")
    if (W != "NULL") {
      cat("\n", "   RI", X, " ~~ RI", W, sep="")
      cat("\n", "   RI", X, " ~~ RS", W, sep="")
      cat("\n", "   RS", X, " ~~ RI", W, sep="")
      cat("\n", "   RS", X, " ~~ RS", W, sep="")
      cat("\n", "   RI", Y, " ~~ RI", W, sep="")
      cat("\n", "   RI", Y, " ~~ RS", W, sep="")
      cat("\n", "   RS", Y, " ~~ RI", W, sep="")
      cat("\n", "   RS", Y, " ~~ RS", W, sep="")
      cat("\n", "   RI", Z, " ~~ RI", W, sep="")
      cat("\n", "   RI", Z, " ~~ RS", W, sep="")
      cat("\n", "   RS", Z, " ~~ RI", W, sep="")
      cat("\n", "   RS", Z, " ~~ RS", W, sep="")
      cat("\n", "   RI", W, " ~~ RS", W, sep="")
    } # end (if W != "NULL")

    # -- Estimate (residual) variance of latent variables -- #
    cat(rep("\n",2), "  # -- Estimate (residual) variance of latent variables -- #")
    for (i in 1:no.waves) {
      cat("\n", paste("  w", X, i, " ~~ eXX", i, "*w", X, i, sep=""))
      cat("\n", paste("  w", Y, i, " ~~ eYY", i, "*w", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  w", Z, i, " ~~ eZZ", i, "*w", Z, i, sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  w", W, i, " ~~ eWW", i, "*w", W, i, sep=""))
      } # end (if W)
    } # end (for i)

    # -- Estimate Covariances Between Residuals of Latent Variables -- #
    cat(rep("\n",2), "  # -- Estimate covariances between residuals of latent variables -- #")
    cat("\n", "  ###################################################################")
    cat("\n", "  # Remove the subscripts for eXY for invariant residual covariance #")
    cat("\n", "  ###################################################################")
    for (i in 2:no.waves) {
      cat("\n", paste("  w", X, i, " ~~ eXY", i, "*w", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  w", X, i, " ~~ eXZ", i, "*w", Z, i, sep=""))
        cat("\n", paste("  w", Y, i, " ~~ eYZ", i, "*w", Z, i, sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  w", X, i, " ~~ eXW", i, "*w", W, i, sep=""))
        cat("\n", paste("  w", Y, i, " ~~ eYW", i, "*w", W, i, sep=""))
        cat("\n", paste("  w", Z, i, " ~~ eZW", i, "*w", W, i, sep=""))
      } # end (if W)
    } # end ((for i)

    # -- Estimate Covariance Between Latent Variables at First Wave -- #
    cat(rep("\n",2), "  # -- Estimate covariance between latent variables at first wave -- #")
    cat("\n", "    w", X, "1 ~~ w", Y, "1", sep="")
    if (Z != "NULL") {
      cat("\n", "    w", X, "1 ~~ w", Z, "1", sep="")
      cat("\n", "    w", Y, "1 ~~ w", Z, "1", sep="")
    }  # end (if Z)
    if (W != "NULL") {
      cat("\n", "    w", X, "1 ~~ w", W, "1", sep="")
      cat("\n", "    w", Y, "1 ~~ w", W, "1", sep="")
      cat("\n", "    w", Z, "1 ~~ w", W, "1", sep="")
    }  # end (if W)

    # -- Estimate covariance among RIX RIY, RIZ, RIW and wx1, wy1, wz1, ww1 -- #
    cat(rep("\n",2), "  # -- Estimate covariance among RIx, RIy, RIz, RIw and wx1, wy1, wz1, ww1 -- #")
    cat("\n", "   RI", X, " ~~ w", X, "1", sep="")
    cat("\n", "   RI", X, " ~~ w", Y, "1", sep="")
    cat("\n", "   RI", Y, " ~~ w", X, "1", sep="")
    cat("\n", "   RI", Y, " ~~ w", Y, "1", sep="")

    if (Z != "NULL") {
      cat("\n", "   RI", X, " ~~ w", Z, "1", sep="")
      cat("\n", "   RI", Y, " ~~ w", Z, "1", sep="")
      cat("\n", "   RI", Z, " ~~ w", X, "1", sep="")
      cat("\n", "   RI", Z, " ~~ w", Y, "1", sep="")
      cat("\n", "   RI", Z, " ~~ w", Z, "1", sep="")
    } # end (if Z != "NULL")
    if (W != "NULL") {
      cat("\n", "   RI", X, " ~~ w", W, "1", sep="")
      cat("\n", "   RI", Y, " ~~ w", W, "1", sep="")
      cat("\n", "   RI", Z, " ~~ w", W, "1", sep="")
      cat("\n", "   RI", W, " ~~ w", X, "1", sep="")
      cat("\n", "   RI", W, " ~~ w", Y, "1", sep="")
      cat("\n", "   RI", W, " ~~ w", Z, "1", sep="")
      cat("\n", "   RI", W, " ~~ w", W, "1", sep="")
    } # end (if W != "NULL")

    # -- Estimate covariance among RSX RSY, RSZ, RSW and wx1, wy1, wz1, ww1 -- #
    cat(rep("\n",2), "  # -- Estimate covariance among RSx, RSy, RSz, RSw and wx1, wy1, wz1, ww1 -- #")
    cat("\n", "   RS", X, " ~~ w", X, "1", sep="")
    cat("\n", "   RS", X, " ~~ w", Y, "1", sep="")
    cat("\n", "   RS", Y, " ~~ w", X, "1", sep="")
    cat("\n", "   RS", Y, " ~~ w", Y, "1", sep="")

    if (Z != "NULL") {
      cat("\n", "   RS", X, " ~~ w", Z, "1", sep="")
      cat("\n", "   RS", Y, " ~~ w", Z, "1", sep="")
      cat("\n", "   RS", Z, " ~~ w", X, "1", sep="")
      cat("\n", "   RS", Z, " ~~ w", Y, "1", sep="")
      cat("\n", "   RS", Z, " ~~ w", Z, "1", sep="")
    } # end (if Z != "NULL")
    if (W != "NULL") {
      cat("\n", "   RS", X, " ~~ w", W, "1", sep="")
      cat("\n", "   RS", Y, " ~~ w", W, "1", sep="")
      cat("\n", "   RS", Z, " ~~ w", W, "1", sep="")
      cat("\n", "   RS", W, " ~~ w", X, "1", sep="")
      cat("\n", "   RS", W, " ~~ w", Y, "1", sep="")
      cat("\n", "   RS", W, " ~~ w", Z, "1", sep="")
      cat("\n", "   RS", W, " ~~ w", W, "1", sep="")
    } # end (if W != "NULL")

    for (j in 1:lag) {
      cat(rep("\n",2), paste0("  # -- Estimate lagged effects between latent variables (Lag = ", j, " wave) -- #", sep =""))
      cat("\n", "  #############################################")
      cat("\n", "  # Remove the subscripts for invariant paths #")
      cat("\n", "  #############################################")
      for (i in (j+1):no.waves) {
        if (Z == "NULL") {
          cat("\n", paste("  w", X,i, " ~ pXX", i,i-j, "*w", X,i-j, " + pYX", i,i-j, "*w", Y,i-j, sep=""))
          cat("\n", paste("  w", Y,i, " ~ pXY", i,i-j, "*w", X,i-j, " + pYY", i,i-j, "*w", Y,i-j, sep=""))
        } else if (W != "NULL") {
          cat("\n", paste("  w", X,i, " ~ pXX", i,i-j, "*w", X,i-j, " + pYX", i,i-j, "*w", Y,i-j, " + pZX", i,i-j, "*w", Z,i-j, " + pWX", i,i-j, "*w", W,i-j, sep=""))
          cat("\n", paste("  w", Y,i, " ~ pXY", i,i-j, "*w", X,i-j, " + pYY", i,i-j, "*w", Y,i-j, " + pZY", i,i-j, "*w", Z,i-j, " + pWY", i,i-j, "*w", W,i-j, sep=""))
          cat("\n", paste("  w", Z,i, " ~ pXZ", i,i-j, "*w", X,i-j, " + pYZ", i,i-j, "*w", Y,i-j, " + pZZ", i,i-j, "*w", Z,i-j, " + pWZ", i,i-j, "*w", W,i-j, sep=""))
          cat("\n", paste("  w", W,i, " ~ pXW", i,i-j, "*w", X,i-j, " + pYW", i,i-j, "*w", Y,i-j, " + pZW", i,i-j, "*w", Z,i-j, " + pWW", i,i-j, "*w", W,i-j, sep=""))
        } else if (Z != "NULL") {
          cat("\n", paste("  w", X,i, " ~ pXX", i,i-j, "*w", X,i-j, " + pYX", i,i-j, "*w", Y,i-j, " + pZX", i,i-j, "*w", Z,i-j, sep=""))
          cat("\n", paste("  w", Y,i, " ~ pXY", i,i-j, "*w", X,i-j, " + pYY", i,i-j, "*w", Y,i-j, " + pZY", i,i-j, "*w", Z,i-j, sep=""))
          cat("\n", paste("  w", Z,i, " ~ pXZ", i,i-j, "*w", X,i-j, " + pYZ", i,i-j, "*w", Y,i-j, " + pZZ", i,i-j, "*w", Z,i-j, sep=""))
        } # end (if Z)
      } # end (for i)
    } # end (for j)

    cat(rep("\n",2), "  ##########################################")
    cat("\n", "  # Regression of indicators on C1 #")
    cat("\n", "  ##########################################")
    for (i in 1:no.waves) {
      cat("\n", paste("  #  ", X, i, " ~ sx", i,"*C1", sep=""))
      cat("\n", paste("  #  ", Y, i, " ~ sy", i,"*C1", sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  #  ", Z, i, " ~ sz", i,"*C1", sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  #  ", W, i, " ~ sz", i,"*C1", sep=""))
      } # end (if W)
    } # end (for i)

    cat(rep("\n",2), "  #########################################")
    cat("\n", "  # Regression of random intercepts on C1 #")
    cat("\n", "  #########################################")
    if (Z == "NULL") {
      cat("\n", "   #  RI", X, " + RI", Y, " ~ C1", sep="")
    } else if (Z != "NULL") {
      cat("\n", "   #  RI", X, " + RI", Y, " + RI", Z, " ~ C1", sep="")
    } else if (W != "NULL") {
      cat("\n", "   #  RI", X, " + RI", Y, " + RI", Z, " + RI", W, " ~ C1", sep="")
    } # end (if Z/W)

    cat(rep("\n",2), "  ################################################################")
    cat("\n", "  # Regression of time-invariant outcome D1 on random intercepts #")
    cat("\n", "  ################################################################")
    if (Z == "NULL") {
      cat("\n", "   #  D1 ~ RI", X, " + RI", Y, sep="")
    } else if (Z != "NULL") {
      cat("\n", "   #  D1 ~ RI", X, " + RI", Y, " + RI", Z, sep="")
    } else if (W != "NULL") {
      cat("\n", "   #  D1 ~ RI", X, " + RI", Y, " + RI", Z, " + RI", W, sep="")
    } # end (if Z/W)
    cat("\n", "  #  D1 ~~ D1 # Residual variance of D1")

    cat(rep("\n",2), "  ################################################")
    cat("\n", "  # Regression of outcome D1 on latent variables #")
    cat("\n", "  ################################################")
    for (i in 1:no.waves) {
      cat("\n", paste("  #  D1 ~ bx", i, "*w", X, i, sep=""))
      cat("\n", paste("  #  D1 ~ by", i, "*w", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  #  D1 ~ bz", i, "*w", Z, i, sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  #  D1 ~ bw", i, "*w", W, i, sep=""))
      } # end (if W)
    } # end (for i)
    cat("\n", "  #  D1 ~~ D1 # Residual variance of D1")

    cat(rep("\n",2), "  '")

    # -- Run ALTMLR -- #
    cat(rep("\n",2), "  ALTMLR.fit <- lavaan::sem(ALT, ")
    cat("\n", "   ", data.source, ",")
    cat("\n", "   missing = 'fiml',")
    cat("\n", "   meanstructure = TRUE,")
    cat("\n", "   information = 'observed',")
    cat("\n", "   estimator = 'MLR'")
    cat("\n", ")")

    # Request summary outputs
    cat(rep("\n",2), "  print(lavaan::summary(ALTMLR.fit, fit.measure = TRUE, standardized = TRUE, rsq = TRUE))", rep("\n",3))

  sink() # Stop writing to file


  ## -- Execute ALT.txt and request summary outputs-- ##
  source('ALT.txt')
  if (lavaan::lavInspect(ALTMLR.fit, what ="post.check") == FALSE) stop("The lavaan solution is non-admissible.")
  ## ------------------------------- ##


  ## -- Monte Carlo Simulation -- ##

  parEst <- lavaan::parameterEstimates(ALTMLR.fit, remove.nonfree = TRUE)
  pest2 <- parEst[, "est"]  # Estimated Parameters
  pest3 <- lavaan::lavTech(ALTMLR.fit, what = "vcov", add.labels = TRUE)  # Estimated Variance-Covariance of Estimated Parameters

  Invariance(parEst, pest2, pest3, no.path, no.waves, lag, Type1, Type1Adj, X, Y, Z, W)

  cat(rep("\n", 2))


}  # end (Function ALT)

## ========================================================================================== ##





# ==================== Creating Function "LGCM" ==================== #
#' Function LGCM (Latent Growth Curve Model)
#'
#' Latent Growth Curve Model (LGCM)
#'
#' @param data.source name of data.frame
#' @param no.waves number of waves (minimum = 3)
#' @param Type1 Overall Type I error rate (default is 0.05) for comparing estimated parameters in the List and Delete method.
#' @param Type1Adj Adjustment of Type I error rate for multiple tests for each estimated parameter. Default is "BON" (Bonferroni adjustment -- Type I error rate/no. of pairwise comparisons), can also be "NULL" (without adjustment).
#' @param X name of variable X.
#' @param Y name of variable Y.
#' @param Z name of variable Z.
#' @param W name of variable W (Z must come before W).
#'
#' @return Latent Growth Curve Model (LGCM) outputs.
#' @export
#' @examples
#'
#' ## -- Example -- ##
#'
#' LGCM(data.source="Data_A", 7, Type1=0.05, Type1Adj="BON", X="EXPOSE.", Y="INTENS.")
#'

LGCM <- function(data.source, no.waves, Type1=0.05, Type1Adj="BON", X, Y, Z="NULL", W = "NULL") {

  ## -- Check inputs -- ##

  if (no.waves < 3) stop("Minimum number of waves is 3")

  if (Type1 > 0.05) stop("Type I error rate > 0.05 is not recommended")
  if (Type1 < 0.0001) stop("Type I error rate < 0.0001 is not recommended")

  if (Z == "NULL" & W != "NULL") stop("Z must be defined before W")

  Type1Adj <- toupper(Type1Adj)
  match.arg(Type1Adj, c("BON", "NULL"))

  ## ----- Creating Model LGCM ----- ###
  sink('LGCM.txt') # Start writing script to LGCM.txt

    cat("\n", "## ----- Specify the model (LGCM) ----- ##", "\n")
    cat("\n", "LGCM <- '")

    # -- Create Between Components (Random Intercepts) -- #
    cat(rep("\n",2), "  # -- Create between components (random intercepts) -- #")
    BX <- paste("  RI", X, " =~ 1*", X, "1", sep="")
    BY <- paste("  RI", Y, " =~ 1*", Y, "1", sep="")
    for (i in 2:no.waves) {
      BX <- paste(BX, " +1*", X, i, sep="")
      BY <- paste(BY, " +1*", Y, i, sep="")
    } # end (for i)
    cat("\n", BX)
    cat("\n", BY)
    if (Z != "NULL") {
      BZ <- paste("  RI", Z, " =~ 1*", Z, "1", sep="")
      for (i in 2:no.waves) {
        BZ <- paste(BZ, " +1*", Z, i, sep="")
      } # end (for i)
      cat("\n", BZ)
    } # end (if Z)
    if (W != "NULL") {
      BW <- paste("  RI", W, " =~ 1*", W, "1", sep="")
      for (i in 2:no.waves) {
        BW <- paste(BW, " +1*", W, i, sep="")
      } # end (for i)
      cat("\n", BW)
    } # end (if W)

    # -- Create Between Components (Random Slopes) -- #
    cat(rep("\n",2), "  # -- Create between components (random slopes) -- #")
    BX <- paste("  RS", X, " =~ 0*", X, "1", sep="")
    BY <- paste("  RS", Y, " =~ 0*", Y, "1", sep="")
    for (i in 2:no.waves) {
      BX <- paste(BX, " +", (i-1), "*", X, i, sep="")
      BY <- paste(BY, " +", (i-1), "*", Y, i, sep="")
    } # end (for i)
    cat("\n", BX)
    cat("\n", BY)
    if (Z != "NULL") {
      BZ <- paste("  RS", Z, " =~ 0*", Z, "1", sep="")
      for (i in 2:no.waves) {
        BZ <- paste(BZ, " +", (i-1), "*", Z, i, sep="")
      } # end (for i)
      cat("\n", BZ)
    } # end (if Z)
    if (W != "NULL") {
      BW <- paste("  RS", W, " =~ 0*", W, "1", sep="")
      for (i in 2:no.waves) {
        BW <- paste(BW, " +", (i-1), "*", W, i, sep="")
      } # end (for i)
      cat("\n", BW)
    } # end (if W)

    # -- Constrain Means (Intercepts) of Indicators to zero -- #
    cat(rep("\n",2), "  # -- Constrain means (intercepts) of indicators to zero -- #")
    for (i in 1:no.waves) {
      cat("\n", paste("  ", X, i, " ~ 0*1", sep=""))
      cat("\n", paste("  ", Y, i, " ~ 0*1", sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  ", Z, i, " ~ 0*1", sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  ", W, i, " ~ 0*1", sep=""))
      } # (if W)
    } # end (for i)

    # -- Estimate Residual Variance of Indicators -- #
    cat(rep("\n",2), "  # -- Estimate residual variance of indicators -- #")
    for (i in 1:no.waves) {
      cat("\n", paste("  ", X, i, " ~~ eIXX", i, "*", X, i, sep=""))
      cat("\n", paste("  ", Y, i, " ~~ eIYY", i, "*", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  ", Z, i, " ~~ eIZZ", i, "*", Z, i, sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  ", W, i, " ~~ eIWW", i, "*", W, i, sep=""))
      } # end (if W)
    } # end (for i)

    # -- Estimate Means (Intercepts) of Random Intercepts -- #
    cat(rep("\n",2), "  # -- Estimate means (intercepts) of random intercepts -- #")
    cat("\n", paste("  RI", X, " ~ MRI", X, "*1", sep=""))
    cat("\n", paste("  RI", Y, " ~ MRI", Y, "*1", sep=""))
    if (Z != "NULL") {
      cat("\n", paste("  RI", Z, " ~ MRI", Z, "*1", sep=""))
    } # end (if Z)
    if (W != "NULL") {
      cat("\n", paste("  RI", W, " ~ MRI", W, "*1", sep=""))
    } # (if W)

    # -- Estimate Means (Intercepts) of Random Slopes -- #
    cat(rep("\n",2), "  # -- Estimate means (intercepts) of random slopes -- #")
    cat("\n", paste("  RS", X, " ~ MRS", X, "*1", sep=""))
    cat("\n", paste("  RS", Y, " ~ MRS", Y, "*1", sep=""))
    if (Z != "NULL") {
      cat("\n", paste("  RS", Z, " ~ MRS", Z, "*1", sep=""))
    } # end (if Z)
    if (W != "NULL") {
      cat("\n", paste("  RS", W, " ~ MRS", W, "*1", sep=""))
    } # (if W)

    # -- Estimate Variance and Covariance of Random Intercepts and Random Slopes -- #
    cat(rep("\n",2), "  # -- Estimate variance and covariance of random intercepts and random slopes -- #")
    cat("\n", "   RI", X, " ~~ RI", X, sep="")
    cat("\n", "   RS", X, " ~~ RS", X, sep="")
    cat("\n", "   RI", Y, " ~~ RI", Y, sep="")
    cat("\n", "   RS", Y, " ~~ RS", Y, sep="")
    if (Z != "NULL") {
      cat("\n", "   RI", Z, " ~~ RI", Z, sep="")
      cat("\n", "   RS", Z, " ~~ RS", Z, sep="")
    } # end (if Z != "NULL")
    if (W != "NULL") {
      cat("\n", "   RI", W, " ~~ RI", W, sep="")
      cat("\n", "   RS", W, " ~~ RS", W, sep="")
    } # end (if W != "NULL")

    cat("\n", "   RI", X, " ~~ RS", X, sep="")
    cat("\n", "   RI", X, " ~~ RI", Y, sep="")
    cat("\n", "   RI", X, " ~~ RS", Y, sep="")
    cat("\n", "   RS", X, " ~~ RI", Y, sep="")
    cat("\n", "   RS", X, " ~~ RS", Y, sep="")
    cat("\n", "   RI", Y, " ~~ RS", Y, sep="")

    if (Z != "NULL") {
      cat("\n", "   RI", X, " ~~ RI", Z, sep="")
      cat("\n", "   RI", X, " ~~ RS", Z, sep="")
      cat("\n", "   RS", X, " ~~ RI", Z, sep="")
      cat("\n", "   RS", X, " ~~ RS", Z, sep="")
      cat("\n", "   RI", Y, " ~~ RI", Z, sep="")
      cat("\n", "   RI", Y, " ~~ RS", Z, sep="")
      cat("\n", "   RS", Y, " ~~ RI", Z, sep="")
      cat("\n", "   RS", Y, " ~~ RS", Z, sep="")
      cat("\n", "   RI", Z, " ~~ RS", Z, sep="")
    } # end (if Z != "NULL")
    if (W != "NULL") {
      cat("\n", "   RI", X, " ~~ RI", W, sep="")
      cat("\n", "   RI", X, " ~~ RS", W, sep="")
      cat("\n", "   RS", X, " ~~ RI", W, sep="")
      cat("\n", "   RS", X, " ~~ RS", W, sep="")
      cat("\n", "   RI", Y, " ~~ RI", W, sep="")
      cat("\n", "   RI", Y, " ~~ RS", W, sep="")
      cat("\n", "   RS", Y, " ~~ RI", W, sep="")
      cat("\n", "   RS", Y, " ~~ RS", W, sep="")
      cat("\n", "   RI", Z, " ~~ RI", W, sep="")
      cat("\n", "   RI", Z, " ~~ RS", W, sep="")
      cat("\n", "   RS", Z, " ~~ RI", W, sep="")
      cat("\n", "   RS", Z, " ~~ RS", W, sep="")
      cat("\n", "   RI", W, " ~~ RS", W, sep="")
    } # end (if W != "NULL")

    cat(rep("\n",2), "  '")

    # -- Run Model LGCMMLR -- #
    cat(rep("\n",2), "# -- Run Model LGCMMLR -- #")
    cat(rep("\n",2), "  LGCMMLR.fit <- suppressWarnings(lavaan::sem(LGCM,")
    cat("\n", "   ", data.source, ",")
    cat("\n", "   missing = 'fiml',")
    cat("\n", "   meanstructure = TRUE,")
    cat("\n", "   information = 'observed',")
    cat("\n", "   estimator = 'MLR'")
    cat("\n", "))")

  sink()  # Stop writing to file LGCM.txt
  ## -------------------------------------- ##


  ## -- Execute LGCM.txt and request summary outputs-- ##
  source('LGCM.txt')
  print(lavaan::summary(LGCMMLR.fit, fit.measure = TRUE, standardized = TRUE, rsq = TRUE))
  if (lavaan::lavInspect(LGCMMLR.fit, what ="post.check") == FALSE) stop("The lavaan solution is non-admissible.")
  ## ------------------------------- ##



  ## -- Monte Carlo Simulation -- ##

  parEst <- lavaan::parameterEstimates(LGCMMLR.fit, remove.nonfree = TRUE)
  pest2 <- parEst[, "est"]  # Estimated Parameters
  pest3 <- lavaan::lavTech(LGCMMLR.fit, what = "vcov", add.labels = TRUE)  # Estimated Variance-Covariance of Estimated Parameters

  Invariance(parEst, pest2, pest3, no.path, no.waves, lag, Type1, Type1Adj, X, Y, Z, W)

  cat(rep("\n", 2))

  ## ----- Start writing script to LGCM.txt ----- ##
  sink('LGCM.txt')
    cat("\n", "# Specify the model (LGCM)", "\n")
    cat("\n", "LGCMSR <- '")

    # -- Create Between Components (Random Intercepts) -- #
    cat(rep("\n",2), "  # -- Create between components (random intercepts) -- #")
    BX <- paste("  RI", X, " =~ 1*", X, "1", sep="")
    BY <- paste("  RI", Y, " =~ 1*", Y, "1", sep="")
    for (i in 2:no.waves) {
      BX <- paste(BX, " +1*", X, i, sep="")
      BY <- paste(BY, " +1*", Y, i, sep="")
    } # end (for i)
    cat("\n", BX)
    cat("\n", BY)
    if (Z != "NULL") {
      BZ <- paste("  RI", Z, " =~ 1*", Z, "1", sep="")
      for (i in 2:no.waves) {
        BZ <- paste(BZ, " +1*", Z, i, sep="")
      } # end (for i)
      cat("\n", BZ)
    } # end (if Z)
    if (W != "NULL") {
      BW <- paste("  RI", W, " =~ 1*", W, "1", sep="")
      for (i in 2:no.waves) {
        BW <- paste(BW, " +1*", W, i, sep="")
      } # end (for i)
      cat("\n", BW)
    } # end (if W)

    # -- Create Between Components (Random Slopes) -- #
    cat(rep("\n",2), "  # -- Create between components (random slopes) -- #")
    BX <- paste("  RS", X, " =~ 0*", X, "1", sep="")
    BY <- paste("  RS", Y, " =~ 0*", Y, "1", sep="")
    for (i in 2:no.waves) {
      BX <- paste(BX, " +", (i-1), "*", X, i, sep="")
      BY <- paste(BY, " +", (i-1), "*", Y, i, sep="")
    } # end (for i)
    cat("\n", BX)
    cat("\n", BY)
    if (Z != "NULL") {
      BZ <- paste("  RS", Z, " =~ 0*", Z, "1", sep="")
      for (i in 2:no.waves) {
        BZ <- paste(BZ, " +", (i-1), "*", Z, i, sep="")
      } # end (for i)
      cat("\n", BZ)
    } # end (if Z)
    if (W != "NULL") {
      BW <- paste("  RS", W, " =~ 0*", W, "1", sep="")
      for (i in 2:no.waves) {
        BW <- paste(BW, " +", (i-1), "*", W, i, sep="")
      } # end (for i)
      cat("\n", BW)
    } # end (if W)

    # -- Constrain Means (Intercepts) of Indicators to zero -- #
    cat(rep("\n",2), "  # -- Constrain means (intercepts) of indicators to zero -- #")
    for (i in 1:no.waves) {
      cat("\n", paste("  ", X, i, " ~ 0*1", sep=""))
      cat("\n", paste("  ", Y, i, " ~ 0*1", sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  ", Z, i, " ~ 0*1", sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  ", W, i, " ~ 0*1", sep=""))
      } # (if W)
    } # end (for i)

    # -- Estimate Residual Variance of Indicators -- #
    cat(rep("\n",2), "  # -- Estimate residual variance of indicators -- #")
    cat("\n", "  ##########################################################")
    cat("\n", "  # Remove the subscripts for invariant indicator variance #")
    cat("\n", "  ##########################################################")
    for (i in 1:no.waves) {
      cat("\n", paste("  ", X, i, " ~~ eIXX", i, "*", X, i, sep=""))
      cat("\n", paste("  ", Y, i, " ~~ eIYY", i, "*", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  ", Z, i, " ~~ eIZZ", i, "*", Z, i, sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  ", W, i, " ~~ eIWW", i, "*", W, i, sep=""))
      } # end (if W)
    } # end (for i)###   ###

    # -- Estimate Means (Intercepts) of Random Intercepts -- #
    cat(rep("\n",2), "  # -- Estimate means (intercepts) of random intercepts -- #")
    cat("\n", paste("  RI", X, " ~ MRI", X, "*1", sep=""))
    cat("\n", paste("  RI", Y, " ~ MRI", Y, "*1", sep=""))
    if (Z != "NULL") {
      cat("\n", paste("  RI", Z, " ~ MRI", Z, "*1", sep=""))
    } # end (if Z)
    if (W != "NULL") {
      cat("\n", paste("  RI", W, " ~ MRI", W, "*1", sep=""))
    } # (if W)

    # -- Estimate Means (Intercepts) of Random Slopes -- #
    cat(rep("\n",2), "  # -- Estimate means (intercepts) of random slopes -- #")
    cat("\n", paste("  RS", X, " ~ MRS", X, "*1", sep=""))
    cat("\n", paste("  RS", Y, " ~ MRS", Y, "*1", sep=""))
    if (Z != "NULL") {
      cat("\n", paste("  RS", Z, " ~ MRS", Z, "*1", sep=""))
    } # end (if Z)
    if (W != "NULL") {
      cat("\n", paste("  RS", W, " ~ MRS", W, "*1", sep=""))
    } # (if W)

    # -- Estimate Variance and Covariance of Random Intercepts and Random Slopes -- #
    cat(rep("\n",2), "  # -- Estimate variance and covariance of random intercepts and random slopes -- #")
    cat("\n", "   RI", X, " ~~ RI", X, sep="")
    cat("\n", "   RS", X, " ~~ RS", X, sep="")
    cat("\n", "   RI", Y, " ~~ RI", Y, sep="")
    cat("\n", "   RS", Y, " ~~ RS", Y, sep="")
    if (Z != "NULL") {
      cat("\n", "   RI", Z, " ~~ RI", Z, sep="")
      cat("\n", "   RS", Z, " ~~ RS", Z, sep="")
    } # end (if Z != "NULL")
    if (W != "NULL") {
      cat("\n", "   RI", W, " ~~ RI", W, sep="")
      cat("\n", "   RS", W, " ~~ RS", W, sep="")
    } # end (if W != "NULL")

    cat("\n", "   RI", X, " ~~ RS", X, sep="")
    cat("\n", "   RI", X, " ~~ RI", Y, sep="")
    cat("\n", "   RI", X, " ~~ RS", Y, sep="")
    cat("\n", "   RS", X, " ~~ RI", Y, sep="")
    cat("\n", "   RS", X, " ~~ RS", Y, sep="")
    cat("\n", "   RI", Y, " ~~ RS", Y, sep="")

    if (Z != "NULL") {
      cat("\n", "   RI", X, " ~~ RI", Z, sep="")
      cat("\n", "   RI", X, " ~~ RS", Z, sep="")
      cat("\n", "   RS", X, " ~~ RI", Z, sep="")
      cat("\n", "   RS", X, " ~~ RS", Z, sep="")
      cat("\n", "   RI", Y, " ~~ RI", Z, sep="")
      cat("\n", "   RI", Y, " ~~ RS", Z, sep="")
      cat("\n", "   RS", Y, " ~~ RI", Z, sep="")
      cat("\n", "   RS", Y, " ~~ RS", Z, sep="")
      cat("\n", "   RI", Z, " ~~ RS", Z, sep="")
    } # end (if Z != "NULL")
    if (W != "NULL") {
      cat("\n", "   RI", X, " ~~ RI", W, sep="")
      cat("\n", "   RI", X, " ~~ RS", W, sep="")
      cat("\n", "   RS", X, " ~~ RI", W, sep="")
      cat("\n", "   RS", X, " ~~ RS", W, sep="")
      cat("\n", "   RI", Y, " ~~ RI", W, sep="")
      cat("\n", "   RI", Y, " ~~ RS", W, sep="")
      cat("\n", "   RS", Y, " ~~ RI", W, sep="")
      cat("\n", "   RS", Y, " ~~ RS", W, sep="")
      cat("\n", "   RI", Z, " ~~ RI", W, sep="")
      cat("\n", "   RI", Z, " ~~ RS", W, sep="")
      cat("\n", "   RS", Z, " ~~ RI", W, sep="")
      cat("\n", "   RS", Z, " ~~ RS", W, sep="")
      cat("\n", "   RI", W, " ~~ RS", W, sep="")
    } # end (if W != "NULL")

    cat(rep("\n",2), "  ##########################################")
    cat("\n", "  # Regression of indicators on C1 #")
    cat("\n", "  ##########################################")
    for (i in 1:no.waves) {
      cat("\n", paste("  #  ", X, i, " ~ sx", i,"*C1", sep=""))
      cat("\n", paste("  #  ", Y, i, " ~ sy", i,"*C1", sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  #  ", Z, i, " ~ sz", i,"*C1", sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  #  ", W, i, " ~ sz", i,"*C1", sep=""))
      } # end (if W)
    } # end (for i)

    cat(rep("\n",2), "  #########################################")
    cat("\n", "  # Regression of random intercepts on C1 #")
    cat("\n", "  #########################################")
    if (Z == "NULL") {
      cat("\n", "   #  RI", X, " + RI", Y, " ~ C1", sep="")
    } else if (Z != "NULL") {
      cat("\n", "   #  RI", X, " + RI", Y, " + RI", Z, " ~ C1", sep="")
    } else if (W != "NULL") {
      cat("\n", "   #  RI", X, " + RI", Y, " + RI", Z, " + RI", W, " ~ C1", sep="")
    } # end (if Z/W)

    cat(rep("\n",2), "  ################################################################")
    cat("\n", "  # Regression of time-invariant outcome D1 on random intercepts #")
    cat("\n", "  ################################################################")
    if (Z == "NULL") {
      cat("\n", "   #  D1 ~ RI", X, " + RI", Y, sep="")
    } else if (Z != "NULL") {
      cat("\n", "   #  D1 ~ RI", X, " + RI", Y, " + RI", Z, sep="")
    } else if (W != "NULL") {
      cat("\n", "   #  D1 ~ RI", X, " + RI", Y, " + RI", Z, " + RI", W, sep="")
    } # end (if Z/W)
    cat("\n", "  #  D1 ~~ D1 # Residual variance of D1")

    cat(rep("\n",2), "  ################################################")
    cat("\n", "  # Regression of outcome D1 on latent variables #")
    cat("\n", "  ################################################")
    for (i in 1:no.waves) {
      cat("\n", paste("  #  D1 ~ bx", i, "*w", X, i, sep=""))
      cat("\n", paste("  #  D1 ~ by", i, "*w", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  #  D1 ~ bz", i, "*w", Z, i, sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  #  D1 ~ bw", i, "*w", W, i, sep=""))
      } # end (if W)
    } # end (for i)
    cat("\n", "  #  D1 ~~ D1 # Residual variance of D1")

    cat(rep("\n",2), "  '")

    # -- Run LGCMMLR -- #
    cat(rep("\n",2), "  LGCMMLR.fit <- lavaan::sem(LGCM, ")
    cat("\n", "   ", data.source, ",")
    cat("\n", "   missing = 'fiml',")
    cat("\n", "   meanstructure = TRUE,")
    cat("\n", "   information = 'observed',")
    cat("\n", "   estimator = 'MLR'")
    cat("\n", ")")

    # Request summary outputs
    cat(rep("\n",2), "  print(lavaan::summary(LGCMMLR.fit, fit.measure = TRUE, standardized = TRUE, rsq = TRUE))", rep("\n", 3))

  sink() # Stop writing to file

}  # end (Function LGCM)

## ========================================================================================== ##




# ==================== Creating Function "LCS" ==================== #
#' Function LCS (Latent Change Score Model)
#'
#' Latent Change Score Model (LCS)
#'
#' @param data.source name of data.frame.
#' @param no.waves number of waves (minimum = 3).
#' @param Type1 Overall Type I error rate (default is 0.05) for comparing estimated parameters in the List and Delete method.
#' @param Type1Adj Adjustment of Type I error rate for multiple tests for each estimated parameter. Default is "BON" (Bonferroni adjustment -- Type I error rate/no. of pairwise comparisons), can also be "NULL" (without adjustment).
#' @param X name of variable X.
#' @param Y name of variable Y.
#' @param Z name of variable Z.
#' @param W name of variable W (Z must come before W).
#'
#' @return Latent Change Score Model (LCS) outputs.
#' @export
#' @examples
#'
#' ## -- Example -- ##
#'
#' LCS(data.source="Data_A", 7, Type1=0.05, Type1Adj="BON", X="EXPOSE.", Y="INTENS.")
#'

LCS <- function(data.source, no.waves, Type1=0.05, Type1Adj="BON", X, Y, Z="NULL", W = "NULL") {

  lag <- 1

  ## -- Check inputs -- ##

  if (no.waves < 3) stop("Minimum number of waves is 3")

  if (Type1 > 0.05) stop("Type I error rate > 0.05 is not recommended")
  if (Type1 < 0.0001) stop("Type I error rate < 0.0001 is not recommended")

  if (Z == "NULL" & W != "NULL") stop("Z must be defined before W")

  Type1Adj <- toupper(Type1Adj)
  match.arg(Type1Adj, c("BON", "NULL"))

  ## ----- Creating Model LCS ----- ###

  sink('LCS.txt')
    cat("\n", "# Specify the model (LCS)", "\n")
    cat("\n", "LCS <- '")

    # -- Create Latent Variables from Indicators -- #
    cat(rep("\n",2), "  # -- Create latent variables -- #")
    for (i in 1:no.waves) {
      cat("\n", paste("  w", X, i, " =~ 1*", X, i, sep=""))
      cat("\n", paste("  w", Y, i, " =~ 1*", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  w", Z, i, " =~ 1*", Z, i, sep=""))
      }  # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  w", W, i, " =~ 1*", W, i, sep=""))
      }  # end (if W)
    } # end (for i)

    # -- Create Latent Change Score Latent Variables -- #
    cat(rep("\n",2), "  # -- Create latent change score latent variables -- #")
    for (i in 2:no.waves) {
      cat("\n", paste("  cs", X, i, " =~ 1*w", X, i, sep=""))
      cat("\n", paste("  cs", Y, i, " =~ 1*w", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  cs", Z, i, " =~ 1*w", Z, i, sep=""))
      }  # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  cs", W, i, " =~ 1*w", W, i, sep=""))
      }  # end (if W)
    } # end (for i)

    # -- Estimate Proportional Change Between Adjacent Latent Variables -- #
    cat(rep("\n",2), "  # -- Estimate proportional change between adjacent latent variables (Lag = 1 wave) -- #")
    cat("\n", "  ###########################################################")
    cat("\n", "  # Remove the subscripts for invariant proportional change #")
    cat("\n", "  ###########################################################")
    for (i in 2:no.waves) {
      cat("\n", paste("  cs", X, i, " ~ dXX", i,i-1, "*w", X,i-1, sep=""))
      cat("\n", paste("  cs", Y, i, " ~ dYY", i,i-1, "*w", Y,i-1, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  cs", Z, i, " ~ dZZ", i,i-1, "*w", Z,i-1, sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  cs", W, i, " ~ dWW", i,i-1, "*w", W,i-1, sep=""))
      } # end (if W)
    } # end (for i)

    # -- Define Change Score Between Adjacent Latent Variables -- #
    cat(rep("\n",2), "  # -- Define change score between adjacent latent variables (Lag = 1 wave) -- #")
    for (i in 2:no.waves) {
      cat("\n", paste("  w", X, i, " ~ 1*w", X,i-1, sep=""))
      cat("\n", paste("  w", Y, i, " ~ 1*w", Y,i-1, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  w", Z, i, " ~ 1*w", Z,i-1, sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  w", W, i, " ~ 1*w", W,i-1, sep=""))
      } # end (if W)
    } # end (for i)

    # -- Estimate Cross-lagged Effects Between Change Scores -- #
    cat(rep("\n",2), "  # -- Estimate cross-lagged effects between change scores (Lag = 1 wave) -- #")
    cat("\n", "  ############################################################")
    cat("\n", "  # Remove the subscripts for invariant cross-lagged effects #")
    cat("\n", "  ############################################################")
    for (i in 2:no.waves) {
      cat("\n", paste("  cs", X, i, " ~ dYX", i,i-1, "*w", Y,i-1, sep=""))
      cat("\n", paste("  cs", Y, i, " ~ dXY", i,i-1, "*w", X,i-1, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  cs", X, i, " ~ dYX", i,i-1, "*w", Y,i-1, " + dZX", i,i-1, "*w", Z,i-1, sep=""))
        cat("\n", paste("  cs", Y, i, " ~ dXY", i,i-1, "*w", X,i-1, " + dZY", i,i-1, "*w", Z,i-1, sep=""))
        cat("\n", paste("  cs", Z, i, " ~ dXZ", i,i-1, "*w", X,i-1, " + dYZ", i,i-1, "*w", Y,i-1, sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  cs", X, i, " ~ dYX", i,i-1, "*w", Y,i-1, " + dZX", i,i-1, "*w", Z,i-1, " + dWX", i,i-1, "*w", W,i-1, sep=""))
        cat("\n", paste("  cs", Y, i, " ~ dXY", i,i-1, "*w", X,i-1, " + dZY", i,i-1, "*w", Z,i-1, " + dWY", i,i-1, "*w", W,i-1, sep=""))
        cat("\n", paste("  cs", Z, i, " ~ dXZ", i,i-1, "*w", X,i-1, " + dYZ", i,i-1, "*w", Y,i-1, " + dWZ", i,i-1, "*w", W,i-1, sep=""))
        cat("\n", paste("  cs", W, i, " ~ dXW", i,i-1, "*w", X,i-1, " + dYW", i,i-1, "*w", Y,i-1, " + dZW", i,i-1, "*w", Z,i-1, sep=""))
      } # end (if W)
    } # end (for i)

    # -- Constrain Intercepts of Indicators to zero -- #
    cat(rep("\n",2), "  # -- Constrain intercepts of indicators to zero -- #")
    for (i in 1:no.waves) {
      cat("\n", paste("  ", X, i, " ~ 0*1", sep=""))
      cat("\n", paste("  ", Y, i, " ~ 0*1", sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  ", Z, i, " ~ 0*1", sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  ", W, i, " ~ 0*1", sep=""))
      } # (if W)
    } # end (for i)

    # -- Create Variance of Indicator Residuals -- #
    cat(rep("\n",2), "  # -- Create variance of indicator residuals -- #")
    cat("\n", "  ###################################################################")
    cat("\n", "  # Remove the subscripts for invariant indicator residual variance #")
    cat("\n", "  ###################################################################")
    for (i in 2:no.waves) {
      cat("\n", paste("  ", X, i, " ~~ eIXX", i, "*", X, i, sep=""))
      cat("\n", paste("  ", Y, i, " ~~ eIYY", i, "*", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  ", Z, i, " ~~ eIZZ", i, "*", Z, i, sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  ", W, i, " ~~ eIWW",  i, "*", W, i, sep=""))
      } # end (if W)
    } # end (for i)###   ###

    # -- Create Covariance of Indicator Residuals -- #
    cat(rep("\n",2), "  # -- Create covariance of indicator residuals -- #")
    cat("\n", "  #####################################################################")
    cat("\n", "  # Remove the subscripts for invariant indicator residual covariance #")
    cat("\n", "  #####################################################################")
    for (i in 2:no.waves) {
      cat("\n", paste("  ", X, i, " ~~ eIXY", i, "*", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  ", X, i, " ~~ eIXZ", i, "*", Z, i, sep=""))
        cat("\n", paste("  ", Y, i, " ~~ eIYZ", i, "*", Z, i, sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  ", X, i, " ~~ eIXW", i, "*", W, i, sep=""))
        cat("\n", paste("  ", Y, i, " ~~ eIYW", i, "*", W, i, sep=""))
        cat("\n", paste("  ", Z, i, " ~~ eIZW", i, "*", W, i, sep=""))
      } # end (if W)
    } # end (for i)###   ###

    # -- Constrain Intercepts of Latent Variables to zero -- #
    cat(rep("\n",2), "  # -- Constrain intercepts of latent variables to zero -- #")
    for (i in 1:no.waves) {
      cat("\n", paste("  w", X, i, " ~ 0*1", sep=""))
      cat("\n", paste("  w", Y, i, " ~ 0*1", sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  w", Z, i, " ~ 0*1", sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  w", W, i, " ~ 0*1", sep=""))
      } # (if W)
    } # end (for i)

    # -- Constrain Residual Variance of Latent Variables to Zero -- #
    cat(rep("\n",2), "  # -- Constrain residual variance of latent variables to zero -- #")
    for (i in 1:no.waves) {
      cat("\n", paste("  w", X, i, " ~~ 0*w", X, i, sep=""))
      cat("\n", paste("  w", Y, i, " ~~ 0*w", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  w", Z, i, " ~~ 0*w", Z, i, sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  w", W, i, " ~~ 0*w", W, i, sep=""))
      } # end (if W)
    } # end (for i)

    # -- Constrain Intercepts of Change Scores to zero -- #
    cat(rep("\n",2), "  # -- Constrain intercepts of changes scores to zero -- #")
    for (i in 2:no.waves) {
      cat("\n", paste("  cs", X, i, " ~ 0*1", sep=""))
      cat("\n", paste("  cs", Y, i, " ~ 0*1", sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  cs", Z, i, " ~ 0*1", sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  cs", W, i, " ~ 0*1", sep=""))
      } # (if W)
    } # end (for i)

    # -- Constrain Residual Variance of Change Scores to Zero -- #
    cat(rep("\n",2), "  # -- Constrain residual variance of change scores to zero -- #")
    for (i in 2:no.waves) {
      cat("\n", paste("  cs", X, i, " ~~ 0*cs", X, i, sep=""))
      cat("\n", paste("  cs", Y, i, " ~~ 0*cs", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  cs", Z, i, " ~~ 0*cs", Z, i, sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  cs", W, i, " ~~ 0*cs", W, i, sep=""))
      } # end (if W)
    } # end (for i)

    # -- Create Initial True Score (Intercept) -- #
    cat(rep("\n",2), "  # -- Create initial true score (intercept) -- #")
    cat("\n", paste("  RI", X, " =~ 1*w", X, "1", sep=""))
    cat("\n", paste("  RI", Y, " =~ 1*w", Y, "1", sep=""))
    if (Z != "NULL") {
      cat("\n", paste("  RI", Z, " =~ 1*w", Z, "1", sep=""))
    } # end (if Z)
    if (W != "NULL") {
      cat("\n", paste("  RI", W, " =~ 1*w", W, "1", sep=""))
    } # end (if W)

    # -- Create Constant Change Estimate -- #
    cat(rep("\n",2), "  # -- Create constant change estimate -- #")
    BX <- paste("  RS", X, " =~ 1*cs", X, "2", sep="")
    BY <- paste("  RS", Y, " =~ 1*cs", Y, "2", sep="")
    for (i in 3:no.waves) {
      BX <- paste(BX, " + 1*cs", X, i, sep="")
      BY <- paste(BY, " + 1*cs", Y, i, sep="")
    } # end (for i)
    cat("\n", BX)
    cat("\n", BY)
    if (Z != "NULL") {
      BZ <- paste("  RS", Z, " =~ 1*cs", Z, "2", sep="")
      for (i in 3:no.waves) {
        BZ <- paste(BZ, " + 1*cs", Z, i, sep="")
      } # end (for i)
      cat("\n", BZ)
    } # end (if Z)
    if (W != "NULL") {
      BW <- paste("  RS", W, " =~ 1*cs", W, "2", sep="")
      for (i in 3:no.waves) {
        BW <- paste(BW, " + 1*cs", W, i, sep="")
      } # end (for i)
      cat("\n", BW)
    } # end (if W)

    # -- Estimate Means of Initial True Score (Intercept) -- #
    cat(rep("\n",2), "  # -- Estimate means of initial true score (intercept) -- #")
    cat("\n", paste("  RI", X, " ~ MRI", X, "*1", sep=""))
    cat("\n", paste("  RI", Y, " ~ MRI", Y, "*1", sep=""))
    if (Z != "NULL") {
      cat("\n", paste("  RI", Z, " ~ MRI", Z, "*1", sep=""))
    } # end (if Z)
    if (W != "NULL") {
      cat("\n", paste("  RI", W, " ~ MRI", W, "*1", sep=""))
    } # (if W)

    # -- Estimate Means of Constant Change Estimate -- #
    cat(rep("\n",2), "  # -- Estimate means of constant change estimates -- #")
    cat("\n", paste("  RS", X, " ~ MRS", X, "*1", sep=""))
    cat("\n", paste("  RS", Y, " ~ MRS", Y, "*1", sep=""))
    if (Z != "NULL") {
      cat("\n", paste("  RS", Z, " ~ MRS", Z, "*1", sep=""))
    } # end (if Z)
    if (W != "NULL") {
      cat("\n", paste("  RS", W, " ~ MRS", W, "*1", sep=""))
    } # (if W)

    # -- Estimate Variance and Covariance of Initial True Score and Constant Change Estimate -- #
    cat(rep("\n",2), "  # -- Estimate variance and covariance of initial true score and constant change estimate -- #")
    cat("\n", "   RI", X, " ~~ RI", X, sep="")
    cat("\n", "   RI", X, " ~~ RS", X, sep="")
    cat("\n", "   RI", X, " ~~ RI", Y, sep="")
    cat("\n", "   RI", X, " ~~ RS", Y, sep="")
    cat("\n", "   RS", X, " ~~ RS", X, sep="")
    cat("\n", "   RS", X, " ~~ RI", Y, sep="")
    cat("\n", "   RS", X, " ~~ RS", Y, sep="")
    cat("\n", "   RI", Y, " ~~ RI", Y, sep="")
    cat("\n", "   RI", Y, " ~~ RS", Y, sep="")
    cat("\n", "   RS", Y, " ~~ RS", Y, sep="")

    if (Z != "NULL") {
      cat("\n", "   RI", X, " ~~ RI", Z, sep="")
      cat("\n", "   RI", X, " ~~ RS", Z, sep="")
      cat("\n", "   RS", X, " ~~ RI", Z, sep="")
      cat("\n", "   RS", X, " ~~ RS", Z, sep="")
      cat("\n", "   RI", Y, " ~~ RI", Z, sep="")
      cat("\n", "   RI", Y, " ~~ RS", Z, sep="")
      cat("\n", "   RS", Y, " ~~ RI", Z, sep="")
      cat("\n", "   RS", Y, " ~~ RS", Z, sep="")
      cat("\n", "   RI", Z, " ~~ RI", Z, sep="")
      cat("\n", "   RI", Z, " ~~ RS", Z, sep="")
      cat("\n", "   RS", Z, " ~~ RS", Z, sep="")
    } # end (if Z != "NULL")
    if (W != "NULL") {
      cat("\n", "   RI", X, " ~~ RI", W, sep="")
      cat("\n", "   RI", X, " ~~ RS", W, sep="")
      cat("\n", "   RS", X, " ~~ RI", W, sep="")
      cat("\n", "   RS", X, " ~~ RS", W, sep="")
      cat("\n", "   RI", Y, " ~~ RI", W, sep="")
      cat("\n", "   RI", Y, " ~~ RS", W, sep="")
      cat("\n", "   RS", Y, " ~~ RI", W, sep="")
      cat("\n", "   RS", Y, " ~~ RS", W, sep="")
      cat("\n", "   RI", Z, " ~~ RI", W, sep="")
      cat("\n", "   RI", Z, " ~~ RS", W, sep="")
      cat("\n", "   RS", Z, " ~~ RI", W, sep="")
      cat("\n", "   RS", Z, " ~~ RS", W, sep="")
      cat("\n", "   RI", W, " ~~ RI", W, sep="")
      cat("\n", "   RI", W, " ~~ RS", W, sep="")
      cat("\n", "   RS", W, " ~~ RS", W, sep="")
    } # end (if W != "NULL")

    cat(rep("\n",2), "  ##########################################")
    cat("\n", "  # Regression of indicators on C1 #")
    cat("\n", "  ##########################################")
    for (i in 1:no.waves) {
      cat("\n", paste("  #  ", X, i, " ~ sx", i,"*C1", sep=""))
      cat("\n", paste("  #  ", Y, i, " ~ sy", i,"*C1", sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  #  ", Z, i, " ~ sz", i,"*C1", sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  #  ", W, i, " ~ sz", i,"*C1", sep=""))
      } # end (if W)
    } # end (for i)

    cat(rep("\n",2), "  ##########################################")
    cat("\n", "  # Regression of initial true score on C1 #")
    cat("\n", "  ##########################################")
    if (Z == "NULL") {
      cat("\n", "   #  RI", X, " + RI", Y, " ~ C1", sep="")
    } else if (Z != "NULL") {
      cat("\n", "   #  RI", X, " + RI", Y, " + RI", Z, " ~ C1", sep="")
    } else if (W != "NULL") {
      cat("\n", "   #  RI", X, " + RI", Y, " + RI", Z, " + RI", W, " ~ C1", sep="")
    } # end (if Z/W)

    cat(rep("\n",2), "  #################################################################")
    cat("\n", "  # Regression of time-invariant outcome D1 on initial true score #")
    cat("\n", "  #################################################################")
    if (Z == "NULL") {
      cat("\n", "   #  D1 ~ RI", X, " + RI", Y, sep="")
    } else if (Z != "NULL") {
      cat("\n", "   #  D1 ~ RI", X, " + RI", Y, " + RI", Z, sep="")
    } else if (W != "NULL") {
      cat("\n", "   #  D1 ~ RI", X, " + RI", Y, " + RI", Z, " + RI", W, sep="")
    } # end (if Z/W)
    cat("\n", "  #  D1 ~~ D1 # Residual variance of D1")

    cat(rep("\n",2), "  ################################################")
    cat("\n", "  # Regression of outcome D1 on latent variables #")
    cat("\n", "  ################################################")
    for (i in 1:no.waves) {
      cat("\n", paste("  #  D1 ~ bx", i, "*w", X, i, sep=""))
      cat("\n", paste("  #  D1 ~ by", i, "*w", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  #  D1 ~ bz", i, "*w", Z, i, sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  #  D1 ~ bw", i, "*w", W, i, sep=""))
      } # end (if W)
    } # end (for i)
    cat("\n", "  #  D1 ~~ D1 # Residual variance of D1")

    cat(rep("\n",2), "  '")

    # -- Run LCSMLR -- #
    cat(rep("\n",2), "  LCSMLR.fit <- lavaan::sem(LCS, ")
    cat("\n", "   ", data.source, ",")
    cat("\n", "   missing = 'fiml',")
    cat("\n", "   meanstructure = TRUE,")
    cat("\n", "   information = 'observed',")
    cat("\n", "   estimator = 'MLR'")
    cat("\n", ")")

    # Request summary outputs
    cat(rep("\n",2), "  print(lavaan::summary(LCSMLR.fit, fit.measure = TRUE, standardized = TRUE, rsq = TRUE))", rep("\n",3))

  sink() # Stop writing to file

  ## -- Execute LCS.txt and request summary outputs-- ##
  source('LCS.txt')

  if (lavaan::lavInspect(LCSMLR.fit, what ="post.check") == FALSE) stop("The lavaan solution is non-admissible.")
  ## ------------------------------- ##



  ## -- Monte Carlo Simulation -- ##

  parEst <- lavaan::parameterEstimates(LCSMLR.fit, remove.nonfree = TRUE)
  pest2 <- parEst[, "est"]  # Estimated Parameters
  pest3 <- lavaan::lavTech(LCSMLR.fit, what = "vcov", add.labels = TRUE)  # Estimated Variance-Covariance of Estimated Parameters

  Invariance(parEst, pest2, pest3, no.path, no.waves, lag, Type1, Type1Adj, X, Y, Z, W)

  cat(rep("\n", 2))

}  # end (Function LCS)

## ========================================================================================== ##





# ==================== Creating Function "LCSCC" ==================== #
#' Function LCSCC (Latent Change Score Model with Changes-To-Changes Extension)
#'
#' Latent Change Score Model with Changes-To-Changes Extension (LCS-CC)
#'
#' @param data.source name of data.frame.
#' @param no.waves number of waves (minimum = 3).
#' @param varI.eq whether indicator residual variances are constrained to be equal (default is false).
#' @param Type1 Overall Type I error rate (default is 0.05) for comparing estimated parameters in the List and Delete method.
#' @param Type1Adj Adjustment of Type I error rate for multiple tests for each estimated parameter. Default is "BON" (Bonferroni adjustment -- Type I error rate/no. of pairwise comparisons), can also be "NULL" (without adjustment).
#' @param X name of variable X.
#' @param Y name of variable Y.
#' @param Z name of variable Z.
#' @param W name of variable W (Z must come before W).
#'
#' @return Latent Change Score Model with Changes-To-Changes Extension (LCSCC) outputs.
#' @export
#' @examples
#'
#' ## -- Example -- ##
#'
#' LCSCC(data.source="Data_A", 7, varI.eq=TRUE, Type1=0.05, Type1Adj="BON", X="EXPOSE.", Y="INTENS.")
#'

LCSCC <- function(data.source, no.waves, varI.eq=FALSE, Type1=0.05, Type1Adj="BON", X, Y, Z="NULL", W = "NULL") {

  lag <- 1

  ## -- Check inputs -- ##

  if (no.waves < 3) stop("Minimum number of waves is 3")

  if (Type1 > 0.05) stop("Type I error rate > 0.05 is not recommended")
  if (Type1 < 0.0001) stop("Type I error rate < 0.0001 is not recommended")

  if (Z == "NULL" & W != "NULL") stop("Z must be defined before W")

  if (is.logical(varI.eq) == FALSE) stop("varI.eq (equivalence of indicator residual variance) can only be TRUE or FALSE")

  Type1Adj <- toupper(Type1Adj)
  match.arg(Type1Adj, c("BON", "NULL"))

  ## ----- Creating Model LCSCC ----- ###

  sink('LCSCC.txt')
    cat("\n", "# Specify the model (LCSCC)", "\n")
    cat("\n", "LCSCC <- '")

    # -- Create Latent Variables from Indicators -- #
    cat(rep("\n",2), "  # -- Create latent variables -- #")
    for (i in 1:no.waves) {
      cat("\n", paste("  w", X, i, " =~ 1*", X, i, sep=""))
      cat("\n", paste("  w", Y, i, " =~ 1*", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  w", Z, i, " =~ 1*", Z, i, sep=""))
      }  # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  w", W, i, " =~ 1*", W, i, sep=""))
      }  # end (if W)
    } # end (for i)

    # -- Create Latent Change Score Latent Variables -- #
    cat(rep("\n",2), "  # -- Create latent change score latent variables -- #")
    for (i in 2:no.waves) {
      cat("\n", paste("  cs", X, i, " =~ 1*w", X, i, sep=""))
      cat("\n", paste("  cs", Y, i, " =~ 1*w", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  cs", Z, i, " =~ 1*w", Z, i, sep=""))
      }  # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  cs", W, i, " =~ 1*w", W, i, sep=""))
      }  # end (if W)
    } # end (for i)

    # -- Estimate Proportional Change Between Adjacent Latent Variables -- #
    cat(rep("\n",2), "  # -- Estimate proportional change between adjacent latent variables (Lag = 1 wave) -- #")
    cat("\n", "  ###########################################################")
    cat("\n", "  # Remove the subscripts for invariant proportional change #")
    cat("\n", "  ###########################################################")
    for (i in 2:no.waves) {
      cat("\n", paste("  cs", X, i, " ~ dXX", i,i-1, "*w", X,i-1, sep=""))
      cat("\n", paste("  cs", Y, i, " ~ dYY", i,i-1, "*w", Y,i-1, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  cs", Z, i, " ~ dZZ", i,i-1, "*w", Z,i-1, sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  cs", W, i, " ~ dWW", i,i-1, "*w", W,i-1, sep=""))
      } # end (if W)
    } # end (for i)

    # -- Define Change Score Between Adjacent Latent Variables -- #
    cat(rep("\n",2), "  # -- Define change score between adjacent latent variables (Lag = 1 wave) -- #")
    for (i in 2:no.waves) {
      cat("\n", paste("  w", X, i, " ~ 1*w", X,i-1, sep=""))
      cat("\n", paste("  w", Y, i, " ~ 1*w", Y,i-1, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  w", Z, i, " ~ 1*w", Z,i-1, sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  w", W, i, " ~ 1*w", W,i-1, sep=""))
      } # end (if W)
    } # end (for i)

    # -- Estimate Cross-lagged Effects Between Change Scores -- #
    cat(rep("\n",2), "  # -- Estimate cross-lagged effects between change scores (Lag = 1 wave) -- #")
    cat("\n", "  ############################################################")
    cat("\n", "  # Remove the subscripts for invariant cross-lagged effects #")
    cat("\n", "  ############################################################")
    for (i in 2:no.waves) {
      cat("\n", paste("  cs", X, i, " ~ dYX", i,i-1, "*w", Y,i-1, sep=""))
      cat("\n", paste("  cs", Y, i, " ~ dXY", i,i-1, "*w", X,i-1, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  cs", X, i, " ~ dYX", i,i-1, "*w", Y,i-1, " + dZX", i,i-1, "*w", Z,i-1, sep=""))
        cat("\n", paste("  cs", Y, i, " ~ dXY", i,i-1, "*w", X,i-1, " + dZY", i,i-1, "*w", Z,i-1, sep=""))
        cat("\n", paste("  cs", Z, i, " ~ dXZ", i,i-1, "*w", X,i-1, " + dYZ", i,i-1, "*w", Y,i-1, sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  cs", X, i, " ~ dYX", i,i-1, "*w", Y,i-1, " + dZX", i,i-1, "*w", Z,i-1, " + dWX", i,i-1, "*w", W,i-1, sep=""))
        cat("\n", paste("  cs", Y, i, " ~ dXY", i,i-1, "*w", X,i-1, " + dZY", i,i-1, "*w", Z,i-1, " + dWY", i,i-1, "*w", W,i-1, sep=""))
        cat("\n", paste("  cs", Z, i, " ~ dXZ", i,i-1, "*w", X,i-1, " + dYZ", i,i-1, "*w", Y,i-1, " + dWZ", i,i-1, "*w", W,i-1, sep=""))
        cat("\n", paste("  cs", W, i, " ~ dXW", i,i-1, "*w", X,i-1, " + dYW", i,i-1, "*w", Y,i-1, " + dZW", i,i-1, "*w", Z,i-1, sep=""))
      } # end (if W)
    } # end (for i)

    # -- Estimate Changes-to-Changes -- #
    cat(rep("\n",2), "  # -- Estimate changes-to-changes (Lag = 1 wave) -- #")
    cat("\n", "  ##########################################################")
    cat("\n", "  # Remove the subscripts for invariant changes-to-changes #")
    cat("\n", "  ##########################################################")
    for (i in 3:no.waves) {
      if (Z == "NULL") {
        cat("\n", paste("  cs", X, i, " ~ ccXX", i,i-1, "*cs", X,i-1, " + ccYX", i,i-1, "*cs", Y,i-1, sep=""))
        cat("\n", paste("  cs", Y, i, " ~ ccYY", i,i-1, "*cs", Y,i-1, " + ccXY", i,i-1, "*cs", X,i-1, sep=""))
      }
      if (Z != "NULL") {
        cat("\n", paste("  cs", X, i, " ~ ccXX", i,i-1, "*cs", X,i-1, " + ccYX", i,i-1, "*cs", Y,i-1, " + ccZX", i,i-1, "*cs", Z,i-1, sep=""))
        cat("\n", paste("  cs", Y, i, " ~ ccXY", i,i-1, "*cs", X,i-1, " + ccYY", i,i-1, "*cs", Y,i-1, " + ccZY", i,i-1, "*cs", Z,i-1, sep=""))
        cat("\n", paste("  cs", Z, i, " ~ ccXZ", i,i-1, "*cs", X,i-1, " + ccYZ", i,i-1, "*cs", Y,i-1, " + ccZZ", i,i-1, "*cs", Z,i-1, sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  cs", X,i, " ~ ccXX",i,i-1, "*cs", X,i-1, " + ccYX",i,i-1, "*cs",Y,i-1, " + ccZX",i,i-1, "*cs",Z,i-1," + ccWX", i,i-1, "*cs",W,i-1,sep=""))
        cat("\n", paste("  cs", Y,i, " ~ ccXY",i,i-1, "*cs", X,i-1, " + ccYY",i,i-1, "*cs",Y,i-1, " + ccZY",i,i-1, "*cs",Z,i-1," + ccWY", i,i-1, "*cs",W,i-1,sep=""))
        cat("\n", paste("  cs", Z,i, " ~ ccXZ",i,i-1, "*cs", X,i-1, " + ccYZ",i,i-1, "*cs",Y,i-1, " + ccZZ",i,i-1, "*cs",Z,i-1," + ccWZ", i,i-1, "*cs",W,i-1,sep=""))
        cat("\n", paste("  cs", W,i, " ~ ccXW",i,i-1, "*cs", X,i-1, " + ccYW",i,i-1, "*cs",Y,i-1, " + ccZW",i,i-1, "*cs",Z,i-1," + ccWW", i,i-1, "*cs",W,i-1,sep=""))
      } # end (if W)
    } # end (for i)

    # -- Constrain Intercepts of Indicators to zero -- #
    cat(rep("\n",2), "  # -- Constrain intercepts of indicators to zero -- #")
    for (i in 1:no.waves) {
      cat("\n", paste("  ", X, i, " ~ 0*1", sep=""))
      cat("\n", paste("  ", Y, i, " ~ 0*1", sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  ", Z, i, " ~ 0*1", sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  ", W, i, " ~ 0*1", sep=""))
      } # (if W)
    } # end (for i)

    # -- Create Variance of Indicator Residuals -- #
    cat(rep("\n",2), "  # -- Create variance of indicator residuals -- #")
    cat("\n", "  ###################################################################")
    cat("\n", "  # Remove the subscripts for invariant indicator residual variance #")
    cat("\n", "  ###################################################################")
    if (isFALSE(varI.eq)) {
      for (i in 2:no.waves) {
        cat("\n", paste("  ", X, i, " ~~ eIXX", i, "*", X, i, sep=""))
        cat("\n", paste("  ", Y, i, " ~~ eIYY", i, "*", Y, i, sep=""))
        if (Z != "NULL") {
          cat("\n", paste("  ", Z, i, " ~~ eIZZ", i, "*", Z, i, sep=""))
        } # end (if Z)
        if (W != "NULL") {
          cat("\n", paste("  ", W, i, " ~~ eIWW", i, "*", W, i, sep=""))
        } # end (if W)
      } # end (for i)
    } else {
      for (i in 2:no.waves) {
        cat("\n", paste("  ", X, i, " ~~ eIXX*", X, i, sep=""))
        cat("\n", paste("  ", Y, i, " ~~ eIYY*", Y, i, sep=""))
        if (Z != "NULL") {
          cat("\n", paste("  ", Z, i, " ~~ eIZZ*", Z, i, sep=""))
        } # end (if Z)
        if (W != "NULL") {
          cat("\n", paste("  ", W, i, " ~~ eIWW*", W, i, sep=""))
        } # end (if W)
      } # end (for i)
    } # end (if varI.eq)

    # -- Create Covariance of Indicator Residuals -- #
    cat(rep("\n",2), "  # -- Create covariance of indicator residuals -- #")
    cat("\n", "  #####################################################################")
    cat("\n", "  # Remove the subscripts for invariant indicator residual covariance #")
    cat("\n", "  #####################################################################")
    for (i in 2:no.waves) {
      cat("\n", paste("  ", X, i, " ~~ eIXY", i, "*", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  ", X, i, " ~~ eIXZ", i, "*", Z, i, sep=""))
        cat("\n", paste("  ", Y, i, " ~~ eIYZ", i, "*", Z, i, sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  ", X, i, " ~~ eIXW", i, "*", W, i, sep=""))
        cat("\n", paste("  ", Y, i, " ~~ eIYW", i, "*", W, i, sep=""))
        cat("\n", paste("  ", Z, i, " ~~ eIZW", i, "*", W, i, sep=""))
      } # end (if W)
    } # end (for i)

    # -- Constrain Intercepts of Latent Variables to zero -- #
    cat(rep("\n",2), "  # -- Constrain intercepts of latent variables to zero -- #")
    for (i in 1:no.waves) {
      cat("\n", paste("  w", X, i, " ~ 0*1", sep=""))
      cat("\n", paste("  w", Y, i, " ~ 0*1", sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  w", Z, i, " ~ 0*1", sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  w", W, i, " ~ 0*1", sep=""))
      } # (if W)
    } # end (for i)

    # -- Constrain Residual Variance of Latent Variables to Zero -- #
    cat(rep("\n",2), "  # -- Constrain residual variance of latent variables to zero -- #")
    for (i in 1:no.waves) {
      cat("\n", paste("  w", X, i, " ~~ 0*w", X, i, sep=""))
      cat("\n", paste("  w", Y, i, " ~~ 0*w", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  w", Z, i, " ~~ 0*w", Z, i, sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  w", W, i, " ~~ 0*w", W, i, sep=""))
      } # end (if W)
    } # end (for i)

    # -- Constrain Intercepts of Change Scores to zero -- #
    cat(rep("\n",2), "  # -- Constrain intercepts of changes scores to zero -- #")
    for (i in 2:no.waves) {
      cat("\n", paste("  cs", X, i, " ~ 0*1", sep=""))
      cat("\n", paste("  cs", Y, i, " ~ 0*1", sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  cs", Z, i, " ~ 0*1", sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  cs", W, i, " ~ 0*1", sep=""))
      } # (if W)
    } # end (for i)

    # -- Constrain Residual Variance of Change Scores to Zero -- #
    cat(rep("\n",2), "  # -- Constrain residual variance of change scores to zero -- #")
    for (i in 2:no.waves) {
      cat("\n", paste("  cs", X, i, " ~~ 0*cs", X, i, sep=""))
      cat("\n", paste("  cs", Y, i, " ~~ 0*cs", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  cs", Z, i, " ~~ 0*cs", Z, i, sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  cs", W, i, " ~~ 0*cs", W, i, sep=""))
      } # end (if W)
    } # end (for i)

    # -- Create Initial True Score (Intercept) -- #
    cat(rep("\n",2), "  # -- Create initial true score (intercept) -- #")
    cat("\n", paste("  RI", X, " =~ 1*w", X, "1", sep=""))
    cat("\n", paste("  RI", Y, " =~ 1*w", Y, "1", sep=""))
    if (Z != "NULL") {
      cat("\n", paste("  RI", Z, " =~ 1*w", Z, "1", sep=""))
    } # end (if Z)
    if (W != "NULL") {
      cat("\n", paste("  RI", W, " =~ 1*w", W, "1", sep=""))
    } # end (if W)

    # -- Create Constant Change Estimate -- #
    cat(rep("\n",2), "  # -- Create constant change estimate -- #")
    BX <- paste("  RS", X, " =~ 1*cs", X, "2", sep="")
    BY <- paste("  RS", Y, " =~ 1*cs", Y, "2", sep="")
    for (i in 3:no.waves) {
      BX <- paste(BX, " + 1*cs", X, i, sep="")
      BY <- paste(BY, " + 1*cs", Y, i, sep="")
    } # end (for i)
    cat("\n", BX)
    cat("\n", BY)
    if (Z != "NULL") {
      BZ <- paste("  RS", Z, " =~ 1*cs", Z, "2", sep="")
      for (i in 3:no.waves) {
        BZ <- paste(BZ, " + 1*cs", Z, i, sep="")
      } # end (for i)
      cat("\n", BZ)
    } # end (if Z)
    if (W != "NULL") {
      BW <- paste("  RS", W, " =~ 1*cs", W, "2", sep="")
      for (i in 3:no.waves) {
        BW <- paste(BW, " + 1*cs", W, i, sep="")
      } # end (for i)
      cat("\n", BW)
    } # end (if W)

    # -- Estimate Means of Initial True Score (Intercept) -- #
    cat(rep("\n",2), "  # -- Estimate means of initial true score (intercept) -- #")
    cat("\n", paste("  RI", X, " ~ MRI", X, "*1", sep=""))
    cat("\n", paste("  RI", Y, " ~ MRI", Y, "*1", sep=""))
    if (Z != "NULL") {
      cat("\n", paste("  RI", Z, " ~ MRI", Z, "*1", sep=""))
    } # end (if Z)
    if (W != "NULL") {
      cat("\n", paste("  RI", W, " ~ MRI", W, "*1", sep=""))
    } # (if W)

    # -- Estimate Means of Constant Change Estimate -- #
    cat(rep("\n",2), "  # -- Estimate means of constant change estimates -- #")
    cat("\n", paste("  RS", X, " ~ MRS", X, "*1", sep=""))
    cat("\n", paste("  RS", Y, " ~ MRS", Y, "*1", sep=""))
    if (Z != "NULL") {
      cat("\n", paste("  RS", Z, " ~ MRS", Z, "*1", sep=""))
    } # end (if Z)
    if (W != "NULL") {
      cat("\n", paste("  RS", W, " ~ MRS", W, "*1", sep=""))
    } # (if W)

    # -- Estimate Variance and Covariance of Initial True Score and Constant Change Estimate -- #
    cat(rep("\n",2), "  # -- Estimate variance and covariance of initial true score and constant change estimate -- #")
    cat("\n", "   RI", X, " ~~ RI", X, sep="")
    cat("\n", "   RI", X, " ~~ RS", X, sep="")
    cat("\n", "   RI", X, " ~~ RI", Y, sep="")
    cat("\n", "   RI", X, " ~~ RS", Y, sep="")
    cat("\n", "   RS", X, " ~~ RS", X, sep="")
    cat("\n", "   RS", X, " ~~ RI", Y, sep="")
    cat("\n", "   RS", X, " ~~ RS", Y, sep="")
    cat("\n", "   RI", Y, " ~~ RI", Y, sep="")
    cat("\n", "   RI", Y, " ~~ RS", Y, sep="")
    cat("\n", "   RS", Y, " ~~ RS", Y, sep="")

    if (Z != "NULL") {
      cat("\n", "   RI", X, " ~~ RI", Z, sep="")
      cat("\n", "   RI", X, " ~~ RS", Z, sep="")
      cat("\n", "   RS", X, " ~~ RI", Z, sep="")
      cat("\n", "   RS", X, " ~~ RS", Z, sep="")
      cat("\n", "   RI", Y, " ~~ RI", Z, sep="")
      cat("\n", "   RI", Y, " ~~ RS", Z, sep="")
      cat("\n", "   RS", Y, " ~~ RI", Z, sep="")
      cat("\n", "   RS", Y, " ~~ RS", Z, sep="")
      cat("\n", "   RI", Z, " ~~ RI", Z, sep="")
      cat("\n", "   RI", Z, " ~~ RS", Z, sep="")
      cat("\n", "   RS", Z, " ~~ RS", Z, sep="")
    } # end (if Z != "NULL")
    if (W != "NULL") {
      cat("\n", "   RI", X, " ~~ RI", W, sep="")
      cat("\n", "   RI", X, " ~~ RS", W, sep="")
      cat("\n", "   RS", X, " ~~ RI", W, sep="")
      cat("\n", "   RS", X, " ~~ RS", W, sep="")
      cat("\n", "   RI", Y, " ~~ RI", W, sep="")
      cat("\n", "   RI", Y, " ~~ RS", W, sep="")
      cat("\n", "   RS", Y, " ~~ RI", W, sep="")
      cat("\n", "   RS", Y, " ~~ RS", W, sep="")
      cat("\n", "   RI", Z, " ~~ RI", W, sep="")
      cat("\n", "   RI", Z, " ~~ RS", W, sep="")
      cat("\n", "   RS", Z, " ~~ RI", W, sep="")
      cat("\n", "   RS", Z, " ~~ RS", W, sep="")
      cat("\n", "   RI", W, " ~~ RI", W, sep="")
      cat("\n", "   RI", W, " ~~ RS", W, sep="")
      cat("\n", "   RS", W, " ~~ RS", W, sep="")
    } # end (if W != "NULL")

    cat(rep("\n",2), "  ##########################################")
    cat("\n", "  # Regression of indicators on C1 #")
    cat("\n", "  ##########################################")
    for (i in 1:no.waves) {
      cat("\n", paste("  #  ", X, i, " ~ sx", i,"*C1", sep=""))
      cat("\n", paste("  #  ", Y, i, " ~ sy", i,"*C1", sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  #  ", Z, i, " ~ sz", i,"*C1", sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  #  ", W, i, " ~ sz", i,"*C1", sep=""))
      } # end (if W)
    } # end (for i)

    cat(rep("\n",2), "  ##########################################")
    cat("\n", "  # Regression of initial true score on C1 #")
    cat("\n", "  ##########################################")
    if (Z == "NULL") {
      cat("\n", "   #  RI", X, " + RI", Y, " ~ C1", sep="")
    } else if (Z != "NULL") {
      cat("\n", "   #  RI", X, " + RI", Y, " + RI", Z, " ~ C1", sep="")
    } else if (W != "NULL") {
      cat("\n", "   #  RI", X, " + RI", Y, " + RI", Z, " + RI", W, " ~ C1", sep="")
    } # end (if Z/W)

    cat(rep("\n",2), "  #################################################################")
    cat("\n", "  # Regression of time-invariant outcome D1 on initial true score #")
    cat("\n", "  #################################################################")
    if (Z == "NULL") {
      cat("\n", "   #  D1 ~ RI", X, " + RI", Y, sep="")
    } else if (Z != "NULL") {
      cat("\n", "   #  D1 ~ RI", X, " + RI", Y, " + RI", Z, sep="")
    } else if (W != "NULL") {
      cat("\n", "   #  D1 ~ RI", X, " + RI", Y, " + RI", Z, " + RI", W, sep="")
    } # end (if Z/W)
    cat("\n", "  #  D1 ~~ D1 # Residual variance of D1")

    cat(rep("\n",2), "  ################################################")
    cat("\n", "  # Regression of outcome D1 on latent variables #")
    cat("\n", "  ################################################")
    for (i in 1:no.waves) {
      cat("\n", paste("  #  D1 ~ bx", i, "*w", X, i, sep=""))
      cat("\n", paste("  #  D1 ~ by", i, "*w", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  #  D1 ~ bz", i, "*w", Z, i, sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  #  D1 ~ bw", i, "*w", W, i, sep=""))
      } # end (if W)
    } # end (for i)
    cat("\n", "  #  D1 ~~ D1 # Residual variance of D1")

    cat(rep("\n",2), "  '")

    # -- Run LCSCCMLR -- #
    cat(rep("\n",2), "  LCSCCMLR.fit <- lavaan::sem(LCSCC, ")
    cat("\n", "   ", data.source, ",")
    cat("\n", "   missing = 'fiml',")
    cat("\n", "   meanstructure = TRUE,")
    cat("\n", "   information = 'observed',")
    cat("\n", "   estimator = 'MLR'")
    cat("\n", ")")

    # Request summary outputs
    cat(rep("\n",2), "  print(lavaan::summary(LCSCCMLR.fit, fit.measure = TRUE, standardized = TRUE, rsq = TRUE))", rep("\n",3))

  sink() # Stop writing to file

  ## -- Execute LCSCC.txt and request summary outputs-- ##
  source('LCSCC.txt')

  if (lavaan::lavInspect(LCSCCMLR.fit, what ="post.check") == FALSE) stop("The lavaan solution is non-admissible.")
  ## ------------------------------- ##



  ## -- Monte Carlo Simulation -- ##

  parEst <- lavaan::parameterEstimates(LCSCCMLR.fit, remove.nonfree = TRUE)
  pest2 <- parEst[, "est"]  # Estimated Parameters
  pest3 <- lavaan::lavTech(LCSCCMLR.fit, what = "vcov", add.labels = TRUE)  # Estimated Variance-Covariance of Estimated Parameters

  Invariance(parEst, pest2, pest3, no.path, no.waves, lag, Type1, Type1Adj, X, Y, Z, W)

  cat(rep("\n", 2))

}  # end (Function LCSCC)

## ========================================================================================== ##





# ==================== Creating Function "GCLM" ==================== #
#' Function GCLM (General Cross-Lagged Panel Model )
#'
#' General Cross-Lagged Panel Model (GCLM)
#'
#' @param data.source name of data.frame
#' @param no.waves number of waves (minimum = 3, must be grater than AR & MA)
#' @param lag number of lags in autoregressive effect (AR) and moving average (MA) (default = 1) AR = MA
#' @param Type1 Overall Type I error rate (default is 0.05) for comparing estimated parameters in the List and Delete method.
#' @param Type1Adj Adjustment of Type I error rate for multiple tests for each estimated parameter. Default is "BON" (Bonferroni adjustment -- Type I error rate/no. of pairwise comparisons), can also be "NULL" (without adjustment).
#' @param X name of variable X.
#' @param Y name of variable Y.
#' @param Z name of variable Z.
#' @param W name of variable W (Z must come before W).
#'
#' @return General Cross-Lagged Panel Model (GCLM) outputs.
#' @export
#' @examples
#'
#' ## -- Example -- ##
#'
#' GCLM(data.source="Data_A", 7, lag=1, Type1=0.05, Type1Adj="BON", X="EXPOSE.", Y="INTENS.")
#'

GCLM <- function(data.source, no.waves, lag=1, Type1=0.05, Type1Adj="BON", X, Y, Z="NULL", W = "NULL") {

  ## -- Check inputs -- ##

  if (no.waves < 3) stop("Minimum number of waves is 3")

  if (lag < 1) stop("Minimum number of lag is 1")
  if (lag > 4) stop("Maximum number of lag is 4")
  if ((no.waves - lag) < 1) stop("Number of waves must be greater than lag plus 1")

  if (Type1 > 0.05) stop("Type I error rate > 0.05 is not recommended")
  if (Type1 < 0.0001) stop("Type I error rate < 0.0001 is not recommended")

  if (Z == "NULL" & W != "NULL") stop("Z must be defined before W")

  Type1Adj <- toupper(Type1Adj)
  match.arg(Type1Adj, c("BON", "NULL"))

  ## ----- Creating Model GCLM ----- ###
  sink('GCLM.txt')
    cat("\n", "# Specify the model (GCLM)", "\n")
    cat("\n", "GCLM <- '")

    # -- Create Latent Variables from Indicators -- #
    cat(rep("\n",2), "  # -- Create latent variables -- #")
    for (i in 1:no.waves) {
      cat("\n", paste("  w", X, i, " =~ 1*", X, i, sep=""))
      cat("\n", paste("  w", Y, i, " =~ 1*", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  w", Z, i, " =~ 1*", Z, i, sep=""))
      }  # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  w", W, i, " =~ 1*", W, i, sep=""))
      }  # end (if W)
    } # end (for i)

    # -- Constrain Residual Variance of Indicators to Zero -- #
    cat(rep("\n",2), "  # -- Constrain residual variance of indicators to zero -- #")
    for (i in 1:no.waves) {
      cat("\n", paste("  ", X, i, " ~~ 0*", X, i, sep=""))
      cat("\n", paste("  ", Y, i, " ~~ 0*", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  ", Z, i, " ~~ 0*", Z, i, sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  ", W, i, " ~~ 0*", W, i, sep=""))
      } # end (if W)
    } # end (for i)

    # -- Create Unit Effects (Random Slopes) from Latent Variables -- #
    cat(rep("\n",2), "  # -- Create unit effects (random slopes) -- #")
    BX <- paste("  RS", X, " =~ 1*w", X, "1", sep="")
    BY <- paste("  RS", Y, " =~ 1*w", Y, "1", sep="")
    for (i in 2:no.waves) {
      BX <- paste(BX, " + uX", i, "*w", X, i, sep="")
      BY <- paste(BY, " + uY", i, "*w", Y, i, sep="")
    } # end (for i)
    cat("\n", BX)
    cat("\n", BY)
    if (Z != "NULL") {
      BZ <- paste("  RS", Z, " =~ 1*w", Z, "1", sep="")
      for (i in 2:no.waves) {
        BZ <- paste(BZ, " + uZ", i, "*w", Z, i, sep="")
      } # end (for i)
      cat("\n", BZ)
    } # end (if Z)
    if (W != "NULL") {
      BW <- paste("  RS", W, " =~ 1*w", W, "1", sep="")
      for (i in 2:no.waves) {
        BW <- paste(BW, " + uW", i, "*w", W, i, sep="")
      } # end (for i)
      cat("\n", BW)
    } # end (if W)

    # -- Constrain Means (Intercepts) of Indicators to Zero -- #
    cat(rep("\n",2), "  # -- Constrain means (intercepts) of indicators to zero -- #")
    for (i in 1:no.waves) {
      cat("\n", paste("  ", X, i, " ~ 0*1", sep=""))
      cat("\n", paste("  ", Y, i, " ~ 0*1", sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  ", Z, i, " ~ 0*1", sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  ", W, i, " ~ 0*1", sep=""))
      } # (if W)
    } # end (for i)

    # -- Create Impulses from Latent Variables -- #
    cat(rep("\n",2), "  # -- Create impulses from latent variables -- #")
    for (i in 1:no.waves) {
      cat("\n", paste("  d", X, i, " =~ 1*w", X, i, sep=""))
      cat("\n", paste("  d", Y, i, " =~ 1*w", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  d", Z, i, " =~ 1*w", Z, i, sep=""))
      }  # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  d", W, i, " =~ 1*w", W, i, sep=""))
      }  # end (if W)
    } # end (for i)

    # -- Constrain Residual Variance of Latent Variables to Zero -- #
    cat(rep("\n",2), "  # -- Constrain residual variance of latent variables to zero -- #")
    for (i in 1:no.waves) {
      cat("\n", paste("  w", X, i, " ~~ 0*w", X, i, sep=""))
      cat("\n", paste("  w", Y, i, " ~~ 0*w", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  w", Z, i, " ~~ 0*w", Z, i, sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  w", W, i, " ~~ 0*w", W, i, sep=""))
      } # end (if W)
    } # end (for i)

    # -- Estimate Intercepts of Latent Variables -- #
    cat(rep("\n",2), "  # -- Estimate intercepts of latent variables -- #")
    for (i in 1:no.waves) {
      cat("\n", paste("  w", X, i, " ~ MX", i, "*1", sep=""))
      cat("\n", paste("  w", Y, i, " ~ MY", i, "*1", sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  w", Z, i, " ~ MZ", i, "*1", sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  w", W, i, " ~ MW", i, "*1", sep=""))
      } # (if W)
    } # end (for i)

    # -- Constrain Means (Intercepts) of Unit Effects to zero -- #
    cat(rep("\n",2), "  # -- Estimate means (intercepts) of unit effects to zero -- #")
    cat("\n", paste("  RS", X, " ~ 0*1", sep=""))
    cat("\n", paste("  RS", Y, " ~ 0*1", sep=""))
    if (Z != "NULL") {
      cat("\n", paste("  RS", Z, " ~ 0*1", sep=""))
    } # end (if Z)
    if (W != "NULL") {
      cat("\n", paste("  RS", W, " ~ 0*1", sep=""))
    } # (if W)

    for (j in 1:lag) {
      cat(rep("\n",2), paste0("  # -- Estimate autoregression (AR) between latent variables (lag = ", j," wave) -- #", sep=""))
      cat("\n", "  ###########################################################")
      cat("\n", "  # Remove the subscripts for invariant autoregression (AR) #")
      cat("\n", "  ###########################################################")
      for (i in (lag+1):no.waves) {
        cat("\n", paste("  w", X,i, " ~ pXX", j, "*w", X,i-j, sep=""))
        cat("\n", paste("  w", Y,i, " ~ pYY", j, "*w", Y,i-j, sep=""))
        if (Z != "NULL") {
          cat("\n", paste("  w", Z,i, " ~ pZZ", j, "*w", Z,i-j, sep=""))
        } # end (if Z)
        if (W != "NULL") {
          cat("\n", paste("  w", W,i, " ~ pWW", j, "*w", W,i-j, sep=""))
        } # end (if W)
      } # end (for i)
    } # end (for j)
    if (lag > 1) {
      for (j in 1:(lag-1)) {
        for (i in (j+1):lag) {
          cat("\n", paste("  w", X,i, " ~ w", X, j, sep=""))
          cat("\n", paste("  w", Y,i, " ~ w", Y, j, sep=""))
          if (Z != "NULL") {
            cat("\n", paste("  w", Z,i, " ~ w", Z, j, sep=""))
          } # end (if Z)
          if (W != "NULL") {
            cat("\n", paste("  w", W,i, " ~ w", W, j, sep=""))
          } # end (if W)
        } # end (for i)
      } # end (for j)
    } # end (if lag)

    for (j in 1:lag) {
      cat(rep("\n",2), paste0("  # -- Estimate moving average (MA) between latent variables (lag = ", j," wave) -- #", sep=""))
      cat("\n", "  ###########################################################")
      cat("\n", "  # Remove the subscripts for invariant moving average (MA) #")
      cat("\n", "  ###########################################################")
      for (i in (lag+1):no.waves) {
        cat("\n", paste("  w", X,i, " ~ maXX", j, "*d", X,i-j, sep=""))
        cat("\n", paste("  w", Y,i, " ~ maYY", j, "*d", Y,i-j, sep=""))
        if (Z != "NULL") {
          cat("\n", paste("  w", Z,i, " ~ maZZ", j, "*d", Z,i-j, sep=""))
        } # end (if Z)
        if (W != "NULL") {
          cat("\n", paste("  w", W,i, " ~ maWW", j, "*d", W,i-j, sep=""))
        } # end (if W)
      } # end (for i)
    } # end (for j)

    cat(rep("\n",2), paste0("  # -- Estimate cross-lagged (CL) effects between latent variables (lag = ", j," wave) -- #", sep=""))
    cat("\n", "  ##########################################")
    cat("\n", "  # Remove the subscripts for invariant CL #")
    cat("\n", "  ##########################################")
    for (i in (lag+1):no.waves) {
      if (Z == "NULL") {
        cat("\n", paste("  w", X,i, " ~ pYX*w", Y,i-1, sep=""))
        cat("\n", paste("  w", Y,i, " ~ pXY*w", X,i-1, sep=""))
      } else if (W != "NULL") {
        cat("\n", paste("  w", X,i, " ~ pYX*w", Y,i-1, " + pZX*w", Z,i-1, " + pWX*w", W,i-1, sep=""))
        cat("\n", paste("  w", Y,i, " ~ pXY*w", X,i-1, " + pZY*w", Z,i-1, " + pWY*w", W,i-1, sep=""))
        cat("\n", paste("  w", Z,i, " ~ pXZ*w", X,i-1, " + pYZ*w", Y,i-1, " + pWZ*w", W,i-1, sep=""))
        cat("\n", paste("  w", W,i, " ~ pXW*w", X,i-1, " + pYW*w", Y,i-1, " + pZW*w", Z,i-1, sep=""))
      } else if (Z != "NULL") {
        cat("\n", paste("  w", X,i, " ~ pYX*w", Y,i-1, " + pZX*w", Z,i-1, sep=""))
        cat("\n", paste("  w", Y,i, " ~ pXY*w", X,i-1, " + pZY*w", Z,i-1, sep=""))
        cat("\n", paste("  w", Z,i, " ~ pXZ*w", X,i-1, " + pYZ*w", Y,i-1, sep=""))
      } # end (if Z)
    } # end (for i)
    if (lag > 1) {
      for (j in 1:(lag-1)) {
        for (i in (j+1):lag) {
          if (Z == "NULL") {
            cat("\n", paste("  w", X,i, " ~ w", Y,i-j, sep=""))
            cat("\n", paste("  w", Y,i, " ~ w", X,i-j, sep=""))
          } else if (W != "NULL") {
            cat("\n", paste("  w", X,i, " ~ w", Y,i-j, " + w", Z,i-j, " + w", W,i-j, sep=""))
            cat("\n", paste("  w", Y,i, " ~ w", X,i-j, " + w", Z,i-j, " + w", W,i-j, sep=""))
            cat("\n", paste("  w", Z,i, " ~ w", X,i-j, " + w", Y,i-j, " + w", W,i-j, sep=""))
            cat("\n", paste("  w", W,i, " ~ w", X,i-j, " + w", Y,i-j, " + w", Z,i-j, sep=""))
          } else if (Z != "NULL") {
            cat("\n", paste("  w", X,i, " ~ w", Y,i-j, " + w", Z,i-j, sep=""))
            cat("\n", paste("  w", Y,i, " ~ w", X,i-j, " + w", Z,i-j, sep=""))
            cat("\n", paste("  w", Z,i, " ~ w", X,i-j, " + w", Y,i-j, sep=""))
          } # end (if Z)
        } # end (for i)
      } # end (for j)
    } # end (if lag)

    cat(rep("\n",2), paste0("  # -- Estimate cross-lagged moving average (CLMA) (Lag = ", j, " wave) -- #", sep=""))
    cat("\n", "  ############################################")
    cat("\n", "  # Remove the subscripts for invariant CLMA #")
    cat("\n", "  ############################################")
    for (i in (lag+1):no.waves) {
      if (Z == "NULL") {
        cat("\n", paste("  w", X,i, " ~ maYX*d", Y, i-1, sep=""))
        cat("\n", paste("  w", Y,i, " ~ maXY*d", X, i-1, sep=""))
      } else if (W != "NULL") {
        cat("\n", paste("  w", X,i, " ~ maYX*d", Y, i-1, " + maZX*d", Z, i-1, " + maWX*d", W, i-1, sep=""))
        cat("\n", paste("  w", Y,i, " ~ maXY*d", X, i-1, " + maZY*d", Z, i-1, " + maWY*d", W, i-1, sep=""))
        cat("\n", paste("  w", Z,i, " ~ maXZ*d", X, i-1, " + maYZ*d", Y, i-1, " + maWZ*d", W, i-1, sep=""))
        cat("\n", paste("  w", W,i, " ~ maXW*d", X, i-1, " + maYW*d", Y, i-1, " + maZW*d", Z, i-1, sep=""))
      } else if (Z != "NULL") {
        cat("\n", paste("  w", X,i, " ~ maYX*d", Y, i-1, " + maZX*d", Z, i-1, sep=""))
        cat("\n", paste("  w", Y,i, " ~ maXY*d", X, i-1, " + maZY*d", Z, i-1, sep=""))
        cat("\n", paste("  w", Z,i, " ~ maXZ*d", X, i-1, " + maYZ*d", Y, i-1, sep=""))
      } # end (if Z)
    } # end (for i)

    # -- Estimate Variance and Covariance of Random Slopes -- #
    cat(rep("\n",2), "  # -- Estimate variance and covariance of random slopes -- #")
    cat("\n", "   RS", X, " ~~ RS", X, sep="")
    cat("\n", "   RS", Y, " ~~ RS", Y, sep="")
    cat("\n", "   RS", X, " ~~ RS", Y, sep="")
    if (Z != "NULL") {
      cat("\n", "   RS", Z, " ~~ RS", Z, sep="")
      cat("\n", "   RS", X, " ~~ RS", Z, sep="")
      cat("\n", "   RS", Y, " ~~ RS", Z, sep="")
    } # end (if Z != "NULL")
    if (W != "NULL") {
      cat("\n", "   RS", W, " ~~ RS", W, sep="")
      cat("\n", "   RS", X, " ~~ RS", W, sep="")
      cat("\n", "   RS", Y, " ~~ RS", W, sep="")
      cat("\n", "   RS", Z, " ~~ RS", W, sep="")
    } # end (if W != "NULL")

    # -- Constrain covariance among RSX RSY, RSZ, RSW and dx, dy, dz, dw -- #
    cat(rep("\n",2), "  # -- Estimate covariance among RSx, RSy, RSz, RSw and dx, dy, dz, dw -- #")
    for (i in 1:no.waves) {
      cat("\n", "   RS", X, " ~~ 0*d", X, i, sep="")
      cat("\n", "   RS", X, " ~~ 0*d", Y, i, sep="")
      cat("\n", "   RS", Y, " ~~ 0*d", X, i, sep="")
      cat("\n", "   RS", Y, " ~~ 0*d", Y, i, sep="")
      if (Z != "NULL") {
        cat("\n", "   RS", X, " ~~ 0*d", Z, i, sep="")
        cat("\n", "   RS", Y, " ~~ 0*d", Z, i, sep="")
        cat("\n", "   RS", Z, " ~~ 0*d", X, i, sep="")
        cat("\n", "   RS", Z, " ~~ 0*d", Y, i, sep="")
        cat("\n", "   RS", Z, " ~~ 0*d", Z, i, sep="")
      } # end (if Z != "NULL")
      if (W != "NULL") {
        cat("\n", "   RS", X, " ~~ 0*d", W, i, sep="")
        cat("\n", "   RS", Y, " ~~ 0*d", W, i, sep="")
        cat("\n", "   RS", Z, " ~~ 0*d", W, i, sep="")
        cat("\n", "   RS", W, " ~~ 0*d", X, i, sep="")
        cat("\n", "   RS", W, " ~~ 0*d", Y, i, sep="")
        cat("\n", "   RS", W, " ~~ 0*d", Z, i, sep="")
        cat("\n", "   RS", W, " ~~ 0*d", W, i, sep="")
      } # end (if W != "NULL")
    } # end (for i)

    # -- Estimate variance of impulses -- #
    cat(rep("\n",2), "  # -- Estimate variance of impulses -- #")
    for (i in 1:no.waves) {
      cat("\n", paste("  d", X, i, " ~~ iXX", i, "*d", X, i, sep=""))
      cat("\n", paste("  d", Y, i, " ~~ iYY", i, "*d", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  d", Z, i, " ~~ iZZ", i, "*d", Z, i, sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  d", W, i, " ~~ iWW", i, "*d", W, i, sep=""))
      } # end (if W)
    } # end (for i)

    # -- Constrain covariance of impulses to zero -- #
    cat(rep("\n",2), "  # -- Constrain covariance of impulses to zero -- #")
    for (j in 1:no.waves) {
      for (i in j:no.waves) {
        if (i != j) {
          cat("\n", paste("  d", X, j, " ~~ 0*d", X, i, sep=""))
          cat("\n", paste("  d", Y, j, " ~~ 0*d", Y, i, sep=""))
          cat("\n", paste("  d", X, j, " ~~ 0*d", Y, i, sep=""))
          cat("\n", paste("  d", Y, j, " ~~ 0*d", X, i, sep=""))
          if (Z != "NULL") {
            cat("\n", paste("  d", X, j, " ~~ 0*d", Z, i, sep=""))
            cat("\n", paste("  d", Y, j, " ~~ 0*d", Z, i, sep=""))
            cat("\n", paste("  d", Z, j, " ~~ 0*d", Z, i, sep=""))
          } # end (if Z)
          if (W != "NULL") {
            cat("\n", paste("  d", X, j, " ~~ 0*d", W, i, sep=""))
            cat("\n", paste("  d", Y, j, " ~~ 0*d", W, i, sep=""))
            cat("\n", paste("  d", Z, j, " ~~ 0*d", W, i, sep=""))
            cat("\n", paste("  d", W, j, " ~~ 0*d", W, i, sep=""))
          } # end (if W)
        } # end (if i != j)
      } # end (for i)
    } # end (for j)

    # -- Estimate Covariances (Co-movements) Between Impulses -- #
    cat(rep("\n",2), "  # -- Estimate covariances (co-movements) between impulses -- #")
    cat("\n", "  ###################################################################")
    cat("\n", "  # Remove the subscripts for eXY for invariant residual covariance #")
    cat("\n", "  ###################################################################")
    for (i in 1:no.waves) {
      cat("\n", paste("  d", X, i, " ~~ iXY", i, "*d", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  d", X, i, " ~~ iXZ", i, "*d", Z, i, sep=""))
        cat("\n", paste("  d", Y, i, " ~~ iYZ", i, "*d", Z, i, sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  d", X, i, " ~~ iXW", i, "*d", W, i, sep=""))
        cat("\n", paste("  d", Y, i, " ~~ iYW", i, "*d", W, i, sep=""))
        cat("\n", paste("  d", Z, i, " ~~ iZW", i, "*d", W, i, sep=""))
      } # end (if W)
    } # end ((for i)

    # -- Define New Parameters -- #
    cat(rep("\n",2), "  # -- Define New Parameters -- #")

    cat(rep("\n",2), "  '")

    # -- Run GCLMMLR -- #
    cat(rep("\n",2), "  GCLMMLR.fit <- lavaan::sem(GCLM, ")
    cat("\n", "   ", data.source, ",")
    cat("\n", "   missing = 'fiml',")
    cat("\n", "   meanstructure = TRUE,")
    cat("\n", "   information = 'observed',")
    cat("\n", "   estimator = 'MLR'")
    cat("\n", ")")

    # Request summary outputs
    cat(rep("\n",2), "  print(lavaan::summary(GCLMMLR.fit, fit.measure = TRUE, standardized = TRUE, rsq = TRUE))", rep("\n",3))

  sink() # Stop writing to file


  ## -- Execute GCLM.txt and request summary outputs-- ##
  source('GCLM.txt')
  if (lavaan::lavInspect(GCLMMLR.fit, what ="post.check") == FALSE) stop("The lavaan solution is non-admissible.")
  ## ------------------------------- ##


  ## -- Monte Carlo Simulation -- ##

  parEst <- lavaan::parameterEstimates(GCLMMLR.fit, remove.nonfree = TRUE)
  pest2 <- parEst[, "est"]  # Estimated Parameters
  pest3 <- lavaan::lavTech(GCLMMLR.fit, what = "vcov", add.labels = TRUE)  # Estimated Variance-Covariance of Estimated Parameters

  Invariance(parEst, pest2, pest3, no.path, no.waves, lag, Type1, Type1Adj, X, Y, Z, W)

  cat(rep("\n", 2))


}  # end (Function GCLM)

## ========================================================================================== ##





# ==================== Create Function "ML" ==================== #
#' Function ML (Multilevel Model)
#'
#' Multilevel Model allows for both MLM and MSEM and compares paths across levels
#'
#'
#' Define the model only once without specifying the level, and the function will compare the model across levels.
#'
#'
#' @param model User-specified measurement model.
#' @param data.source A data frame containing the observed variables used in the model.
#' @param Cluster Cluster variable for nested data. The Monte Carlo simulation method should be used for nested data.
#' @param missing missing = "listwise" or "fiml". lavaan uses listwise deletion for missing cases by default. missing = "fiml" in lavaan multilevel model is new and may have convergence problems.
#' @param L2 If the model should be replicated at L2 (default = TRUE). If L2=FALSE, covariance among variables will be estimated at L2.
#' @param mL2.variables L2 variables represented by group mean (gm_X, manifest L2 variables) and L1 variables represented by group-mean centered variables (gmc_X). By default, variables are partitioned into L1 and L2 by multilevel SEM.
#'
#' @return Estimated parameters and differences in estimated parameters across levels
#' @export
#' @examples
#'
#' ## == Example == ##
#' # Data file is "Data_A"; cluster variable is "ID"
#'
#' # Convert wide-format date file to long-format
#' # df_long <- wide2long(Data_A, no.waves=7, variables=c("EDU.", "AGE.", "BMI_S.",
#' # "EXPOSE.", "INTENS."), lag1=FALSE)
#'
#' ## Specify the path model - Model ##
#' Model <- '
#'  INTENS. ~ EDU. + AGE. + BMI_S.
#'  EXPOSE. ~ EDU. + AGE. + BMI_S.
#' '
#'
#' ## ===== Compare Estimated Parameters Across Levels ===== ##
#' ML(Model, df_long, Cluster = "id", L2=TRUE, mL2.variables=c("EDU.", "AGE.", "BMI_S."))
#'
ML <- function(model, data.source, Cluster="NULL", missing="listwise", L2=TRUE, mL2.variables=NULL) {

  arg1_char <- deparse(substitute(model))
  arg2_char <- deparse(substitute(data.source))
  arg3_char <- deparse(substitute(Cluster))
  missing <- tolower(missing)
  match.arg(missing, c("listwise","fiml"))

  # -- Run Single Level to get model information -- #
  Model.L1 <- lavaan::sem(model,
                          data.source,
                          cluster=Cluster,
                          estimator = 'MLR')

  temp <- lavaan::parameterEstimates(Model.L1)
  no.estimates <- nrow(temp)

  # -- Create manifest Level 2 variables and group-mean centered variables -- #
  if (is.null(mL2.variables) == FALSE) {
    mL2.data <<- mL2(data.source, Cluster, mL2.variables)
  } # end (if mL2.variables)

  # -- Create and Run Multilevel Model -- #
    ML.X <- paste0("Model.ML.X <- '", "\n", "level: 1", "\n")
    for (i in 1:no.estimates) {
      var.lhs <- temp[i, "lhs"]
      var.rhs <- temp[i, "rhs"]
      if (temp[i, "lhs"] %in% mL2.variables)  var.lhs <- paste0("gmc_", temp[i, "lhs"])
      if (temp[i, "rhs"] %in% mL2.variables)  var.rhs <- paste0("gmc_", temp[i, "rhs"])
      ML.X <- paste0(ML.X, var.lhs, temp[i, "op"], var.rhs, "\n")
    } # end (for i)
    if (L2 == TRUE) {
      ML.X <- paste0(ML.X, paste0("\n", "level: 2", "\n"))
      for (i in 1:no.estimates) {
        var.lhs <- temp[i, "lhs"]
        var.rhs <- temp[i, "rhs"]
        if (temp[i, "lhs"] %in% mL2.variables)  var.lhs <- paste0("gm_", temp[i, "lhs"])
        if (temp[i, "rhs"] %in% mL2.variables)  var.rhs <- paste0("gm_", temp[i, "rhs"])
        ML.X <- paste0(ML.X, var.lhs, temp[i, "op"], var.rhs, "\n")
      } # end (for i)
    } else {
      unique_variables <- unique(temp[,"lhs"])
      for (i in 1:length(unique_variables)) {
        if (unique_variables[i] %in% mL2.variables) {
          unique_variables[i] <- paste0("gm_", unique_variables[i])
        } # end (if unique_variables)
      } # end (for i)
      ML.X <- paste0(ML.X, paste0("\n", "level: 2", "\n"))
      for (j in 1:(length(unique_variables)-1)) {
        for (i in (j+1):(length(unique_variables))) {
          ML.X <- paste0(ML.X, unique_variables[j], "~~", unique_variables[i], "\n")
        } # end (for i)
      } # end (for j)
    } # end (if L2 = TRUE)
    ML.X <- paste0(ML.X, "'","\n")
    eval(parse(text = ML.X))

  sink('ML.txt')
    cat(ML.X)
    # -- Run Model.ML.X -- #
    cat(rep("\n",2), "  Model.L2.fit <- suppressWarnings(lavaan::sem(Model.ML.X,")
    if (is.null(mL2.variables) == FALSE) {
      cat("\n", "   data = mL2.data,")
    } else {
      cat("\n", paste0("    data = ", arg2_char,","))
    }
    if (missing == "fiml") { cat("\n", "   missing = 'fiml',") }
    cat("\n", "   information = 'observed',")
    cat("\n", "   estimator = 'MLR',")
    cat("\n", paste0("   cluster = ", arg3_char, ","))
    cat("\n", "   verbose = FALSE))")

    # Request summary outputs
    cat(rep("\n",2), "  print(lavaan::summary(Model.L2.fit, fit.measure = TRUE, standardized = TRUE, rsq = TRUE))", rep("\n",3))

  sink() # Stop writing to file

  ## -- Execute ML.txt and request summary outputs-- ##
  source('ML.txt')
  if (lavaan::lavInspect(Model.L2.fit, what ="post.check") == FALSE) stop("The lavaan solution is non-admissible.")
  ## ------------------------------- ##


  ## -- Estimate ICC -- ##
  est.var1 <- subset(lavaan::parameterEstimates(Model.L2.fit), op == "~~" & (lhs == rhs) & level == 1)
  est.var2 <- subset(lavaan::parameterEstimates(Model.L2.fit), op == "~~" & (lhs == rhs) & level == 2)
  unique_variables <- unique(temp[,"lhs"])
  cluster.size <- mean(lavaan::lavInspect(Model.L2.fit, what = "cluster.size"))
  no.variables <- length(unique_variables)
  ICC <- matrix(0, nrow=no.variables, ncol=7)
  colnames(ICC) <- c("Variable", "L1 label", "L2 label", "L1 variance", "L2 variance", "ICC(1)", "ICC(2)")
  rownames(ICC) <- unique_variables
  for (i in 1:no.variables) {
    ICC[i, 1] <- unique_variables[i]
    if (ICC[i, 1] %in% mL2.variables) {
      ICC[i, 2] <- paste0("gmc_", ICC[i, 1])
      ICC[i, 3] <- paste0("gm_", ICC[i, 1])
    } else {
      ICC[i, 2] <- unique_variables[i]
      ICC[i, 3] <- unique_variables[i]
    } # end (if ICC[i, 1])
    ICC[i, 4] <- round(est.var1[which(est.var1[, "lhs"] == ICC[i, 2]), "est"], digits=4)
    ICC[i, 5] <- round(est.var2[which(est.var1[, "lhs"] == ICC[i, 2]), "est"], digits=4)
    ICC[i, 6] <- round(as.numeric(ICC[i, 5]) / (as.numeric(ICC[i, 4]) + as.numeric(ICC[i, 5])), digits=4)
    ICC[i, 7] <- round(cluster.size*as.numeric(ICC[i, 6]) / (1 + (cluster.size-1)*as.numeric(ICC[i, 6])), digits=4)
  } # end (for i)
  cat("\n", "## ----- ICC(1) and ICC(2) ----- ##", "\n")
  print(as.data.frame(ICC),row.names=F)
  cat(rep("\n",2))

  ## -- Monte Carlo Simulation -- ##

  parEst.L1 <- subset(lavaan::parameterEstimates(Model.L2.fit, remove.nonfree = TRUE), level == 1)
  no.estimates.L1 <- nrow(parEst.L1)
  parEst <- lavaan::parameterEstimates(Model.L2.fit, remove.nonfree = TRUE)
  pest2 <- parEst[, "est"]  # Estimated Parameters
  pest3 <- lavaan::lavTech(Model.L2.fit, what = "vcov", add.labels = TRUE)  # Estimated Variance-Covariance of Estimated Parameters

  mcmc <- MASS::mvrnorm(n=1000000, mu=pest2, Sigma=pest3, tol = 1e-6)  # Run 1,000,000 simulations
  names(pest2) <-colnames(pest3)  # Save Parameter Names to Estimated Parameters
  b.no <- nrow(mcmc)  # No. of successful simulated samples

  ## -- Differences in estimated parameters -- ##
  for (j in 1:no.estimates.L1) {
    mcmcA <- (mcmc[, j] - mcmc[, j+no.estimates.L1])
    mcmc <- cbind(mcmc, mcmcA)
    colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste0("D", j)
    pest2A <- (pest2[j] - pest2[j+no.estimates.L1])
    names(pest2A) <- paste0("D", j)
    pest2 <- append(pest2, pest2A)
  } # end (for j)
  ## ----- end (Difference in estimated parameters ) ----- ##

  ## -- Print Results -- ##
  cat("\n", "## ----- Differences in Estimated Parameters Across Levels ----- ##")
  for (j in 1:no.estimates.L1) {
    cat("\n")
    print(lavaan::parameterEstimates(Model.L2.fit)[c(j, j+no.estimates),], row.names=FALSE)

    estM <- pest2[paste0("D", j)] # estimated parameter
    abM <- mcmc[, paste0("D", j)] # simulated parameter

    # Calculate Percentile Probability
    if (quantile(abM,probs=0.5)>0) {
      pPCI = 2*(sum(abM<0)/b.no)
    } else {
      pPCI = 2*(sum(abM>0)/b.no)
    } # end (if abM)
    cat(" Difference in estimated parameter = ", format(round(estM, 4), nsmall = 4), "p = ", format(round(pPCI, 4), nsmall = 4), "\n")
  } # end (for j)
  cat("\n")

#  list(lavaan.object = Model.L2.fit)

} ## End (Function ML)

# ==================== Finish Function "ML" ==================== #





# ==================== Creating Function "wide2long" ==================== #
#' Function wide2long (Convert Wide Data Format to Long Data Format)
#'
#' Convert wide data format to long data format with the capability to create lagged (t-1) variables. Variables should be coded with ".1, .2" to indicate the timelines. A new "id" variable will be created as the grouping variable, and a new "time" variable will be created as the timing variable.
#'
#' @param data.source name of dataframe.
#' @param no.waves number of waves of data.
#' @param variables variable names to create long data format.
#' @param lag1 whether to create lagged (t - 1) data of the variables. Default is FALSE.
#'
#' @return object as a dataframe.
#' @export
#' @examples
#'
#' ## -- Example -- ##
#'
#'
#' # Suppose df_long is the new dataframe name.
#' # Suppose data are "EXPOSE.1", "EXPOSE.2", "EXPOSE.3", ... "EXPOSE.7" in Data_A
#'
#' df_long <- wide2long(Data_A, 7, variables=c("EXPOSE.", "INTENS."), lag1=FALSE)
#' write.csv(df_long, "file_long.csv") # save dataframe to csv file
#'

## ===== Convert data file to long format ===== ##
wide2long <- function(data.source, no.waves, variables = c("X", "Y"), lag1=FALSE) {

  # -- Check inputs -- #
  if (is.logical(lag1) == FALSE) stop("lag1 (create time-lagged (t-1) variables can only be TRUE or FALSE")

  if (lag1 == FALSE) {
    df_long <- suppressWarnings(reshape(data.source,
      direction = "long",
      idvar = "id",
      varying = lapply(variables, function(x) paste0(x, 1:no.waves)),
      v.names = variables,
      timevar = "time",
      times = c(1: no.waves) # Optional: specify time values explicitly
    ))
  } else {
    df_long <- suppressWarnings(reshape(data.source,
      direction = "long",
      idvar = "id",
      varying = lapply(variables, function(x) paste0(x, 1:no.waves)),
      v.names = variables,
      timevar = "time",
      times = c(1: no.waves) # Optional: specify time values explicitly
    ))
    # -- Create Lag 1 Variables -- #
    df_long_lagged <- df_long %>%
    arrange(id, time) %>% # Ensure data is sorted by group and time
    group_by(id) %>%     # Group by the identifier
    mutate(across(all_of(variables), ~lag(.x, n = 1), .names = "{.col}_lag")) %>%
    ungroup()            # Ungroup the data when done
  } # end (if lag1)
} # end (function wide2long)

## ========================================================================================== ##





# ==================== Creating Function "long2wide" ==================== #
#' Function long2wide (Convert Long Data Format to Wide Data Format)
#'
#' Convert long data format to wide data format. New variables will be coded with ".1, .2" to indicate the timelines.
#'
#' @param data.source name of dataframe.
#' @param id variable indicates group.
#' @param time variable indicates time.
#'
#' @return object as a dataframe.
#' @export
#' @examples
#'
#' ## -- Example -- ##
#'
#' df_wide <- long2wide(df_long, "id", "time", variables=c("EXPOSE.", "INTENS."))
#' write.csv(df_wide, "file_wide.csv") # save dataframe to csv file
#'

## ===== Convert data file to wide format ===== ##
long2wide <- function(data.source, id="id", time="time", variables = c("X", "Y")) {
  df_wide <- suppressWarnings(reshape(
    data = data.source,
    idvar = id,
    timevar = time,
    v.names = variables,
    direction = "wide"
  ))
} # end (function long2wide)

## ========================================================================================== ##





