options("width"=210)


## ----- Sub-Function Invariance for testing invariance of parameters ----- ##

Invariance <- function(parEst, pest2, pest3, no.path, MIset, no.compare, no.waves, lag, varI.eq, p, X, Y, Z, W) {

  mcmc <- MASS::mvrnorm(n=1000000, mu=pest2, Sigma=pest3, tol = 1e-6)  # Run 1,000,000 simulations
  names(pest2) <-colnames(pest3)  # Save Parameter Names to Estimated Parameters
  b.no <- nrow(mcmc)  # No. of successful simulated samples

  ## -- Lag = 1 -- ##
  for (j in 3:no.waves) {
    for (i in j:no.waves) {
      mcmcA <- (mcmc[, paste("XX", i, i-1, sep="")] - mcmc[, paste("XX", i-j+2, i-j+1, sep="")])
      mcmc <- cbind(mcmc, mcmcA)
      colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("XX", i, i-1, "-XX", i-j+2, i-j+1, sep="")
      pest2A <- (pest2[paste("XX", i, i-1, sep="")] - pest2[paste("XX", i-j+2, i-j+1, sep="")])
      names(pest2A) <- paste("XX", i, i-1, "-XX", i-j+2, i-j+1, sep="")
      pest2 <- append(pest2, pest2A)

      mcmcA <- (mcmc[, paste("YY", i, i-1, sep="")] - mcmc[, paste("YY", i-j+2, i-j+1, sep="")])
      mcmc <- cbind(mcmc, mcmcA)
      colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("YY", i, i-1, "-YY", i-j+2, i-j+1, sep="")
      pest2A <- (pest2[paste("YY", i, i-1, sep="")] - pest2[paste("YY", i-j+2, i-j+1, sep="")])
      names(pest2A) <- paste("YY", i, i-1, "-YY", i-j+2, i-j+1, sep="")
      pest2 <- append(pest2, pest2A)

      mcmcA <- (mcmc[, paste("XY", i, i-1, sep="")] - mcmc[, paste("XY", i-j+2, i-j+1, sep="")])
      mcmc <- cbind(mcmc, mcmcA)
      colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("XX", i, i-1, "-XY", i-j+2, i-j+1, sep="")
      pest2A <- (pest2[paste("XY", i, i-1, sep="")] - pest2[paste("XY", i-j+2, i-j+1, sep="")])
      names(pest2A) <- paste("XY", i, i-1, "-XY", i-j+2, i-j+1, sep="")
      pest2 <- append(pest2, pest2A)

      mcmcA <- (mcmc[, paste("YX", i, i-1, sep="")] - mcmc[, paste("YX", i-j+2, i-j+1, sep="")])
      mcmc <- cbind(mcmc, mcmcA)
      colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("YX", i, i-1, "-YX", i-j+2, i-j+1, sep="")
      pest2A <- (pest2[paste("YX", i, i-1, sep="")] - pest2[paste("YX", i-j+2, i-j+1, sep="")])
      names(pest2A) <- paste("YX", i, i-1, "-YX", i-j+2, i-j+1, sep="")
      pest2 <- append(pest2, pest2A)

      if (Z != "NULL") {
        mcmcA <- (mcmc[, paste("ZZ", i, i-1, sep="")] - mcmc[, paste("ZZ", i-j+2, i-j+1, sep="")])
        mcmc <- cbind(mcmc, mcmcA)
        colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("ZZ", i, i-1, "-ZZ", i-j+2, i-j+1, sep="")
        pest2A <- (pest2[paste("ZZ", i, i-1, sep="")] - pest2[paste("ZZ", i-j+2, i-j+1, sep="")])
        names(pest2A) <- paste("ZZ", i, i-1, "-ZZ", i-j+2, i-j+1, sep="")
        pest2 <- append(pest2, pest2A)

        mcmcA <- (mcmc[, paste("XZ", i, i-1, sep="")] - mcmc[, paste("XZ", i-j+2, i-j+1, sep="")])
        mcmc <- cbind(mcmc, mcmcA)
        colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("XZ", i, i-1, "-XZ", i-j+2, i-j+1, sep="")
        pest2A <- (pest2[paste("XZ", i, i-1, sep="")] - pest2[paste("XZ", i-j+2, i-j+1, sep="")])
        names(pest2A) <- paste("XZ", i, i-1, "-XZ", i-j+2, i-j+1, sep="")
        pest2 <- append(pest2, pest2A)

        mcmcA <- (mcmc[, paste("YZ", i, i-1, sep="")] - mcmc[, paste("YZ", i-j+2, i-j+1, sep="")])
        mcmc <- cbind(mcmc, mcmcA)
        colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("YZ", i, i-1, "-YZ", i-j+2, i-j+1, sep="")
        pest2A <- (pest2[paste("YZ", i, i-1, sep="")] - pest2[paste("YZ", i-j+2, i-j+1, sep="")])
        names(pest2A) <- paste("YZ", i, i-1, "-YZ", i-j+2, i-j+1, sep="")
        pest2 <- append(pest2, pest2A)

        mcmcA <- (mcmc[, paste("ZX", i, i-1, sep="")] - mcmc[, paste("ZX", i-j+2, i-j+1, sep="")])
        mcmc <- cbind(mcmc, mcmcA)
        colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("ZX", i, i-1, "-ZX", i-j+2, i-j+1, sep="")
        pest2A <- (pest2[paste("ZX", i, i-1, sep="")] - pest2[paste("ZX", i-j+2, i-j+1, sep="")])
        names(pest2A) <- paste("ZX", i, i-1, "-ZX", i-j+2, i-j+1, sep="")
        pest2 <- append(pest2, pest2A)

        mcmcA <- (mcmc[, paste("ZY", i, i-1, sep="")] - mcmc[, paste("ZY", i-j+2, i-j+1, sep="")])
        mcmc <- cbind(mcmc, mcmcA)
        colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("ZY", i, i-1, "-ZY", i-j+2, i-j+1, sep="")
        pest2A <- (pest2[paste("ZY", i, i-1, sep="")] - pest2[paste("ZY", i-j+2, i-j+1, sep="")])
        names(pest2A) <- paste("ZY", i, i-1, "-ZY", i-j+2, i-j+1, sep="")
        pest2 <- append(pest2, pest2A)
      } # end (if Z != "NULL")

      if (W != "NULL") {
        mcmcA <- (mcmc[, paste("WW", i, i-1, sep="")] - mcmc[, paste("WW", i-j+2, i-j+1, sep="")])
        mcmc <- cbind(mcmc, mcmcA)
        colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("WW", i, i-1, "-WW", i-j+2, i-j+1, sep="")
        pest2A <- (pest2[paste("WW", i, i-1, sep="")] - pest2[paste("WW", i-j+2, i-j+1, sep="")])
        names(pest2A) <- paste("WW", i, i-1, "-WW", i-j+2, i-j+1, sep="")
        pest2 <- append(pest2, pest2A)

        mcmcA <- (mcmc[, paste("XW", i, i-1, sep="")] - mcmc[, paste("XW", i-j+2, i-j+1, sep="")])
        mcmc <- cbind(mcmc, mcmcA)
        colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("XW", i, i-1, "-XW", i-j+2, i-j+1, sep="")
        pest2A <- (pest2[paste("XW", i, i-1, sep="")] - pest2[paste("XW", i-j+2, i-j+1, sep="")])
        names(pest2A) <- paste("XW", i, i-1, "-XW", i-j+2, i-j+1, sep="")
        pest2 <- append(pest2, pest2A)

        mcmcA <- (mcmc[, paste("YW", i, i-1, sep="")] - mcmc[, paste("YW", i-j+2, i-j+1, sep="")])
        mcmc <- cbind(mcmc, mcmcA)
        colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("YW", i, i-1, "-YW", i-j+2, i-j+1, sep="")
        pest2A <- (pest2[paste("YW", i, i-1, sep="")] - pest2[paste("YW", i-j+2, i-j+1, sep="")])
        names(pest2A) <- paste("YW", i, i-1, "-YW", i-j+2, i-j+1, sep="")
        pest2 <- append(pest2, pest2A)

        mcmcA <- (mcmc[, paste("ZW", i, i-1, sep="")] - mcmc[, paste("ZW", i-j+2, i-j+1, sep="")])
        mcmc <- cbind(mcmc, mcmcA)
        colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("ZW", i, i-1, "-ZW", i-j+2, i-j+1, sep="")
        pest2A <- (pest2[paste("ZW", i, i-1, sep="")] - pest2[paste("ZW", i-j+2, i-j+1, sep="")])
        names(pest2A) <- paste("ZW", i, i-1, "-ZW", i-j+2, i-j+1, sep="")
        pest2 <- append(pest2, pest2A)

        mcmcA <- (mcmc[, paste("WX", i, i-1, sep="")] - mcmc[, paste("WX", i-j+2, i-j+1, sep="")])
        mcmc <- cbind(mcmc, mcmcA)
        colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("WX", i, i-1, "-WX", i-j+2, i-j+1, sep="")
        pest2A <- (pest2[paste("WX", i, i-1, sep="")] - pest2[paste("WX", i-j+2, i-j+1, sep="")])
        names(pest2A) <- paste("WX", i, i-1, "-WX", i-j+2, i-j+1, sep="")
        pest2 <- append(pest2, pest2A)

        mcmcA <- (mcmc[, paste("WY", i, i-1, sep="")] - mcmc[, paste("WY", i-j+2, i-j+1, sep="")])
        mcmc <- cbind(mcmc, mcmcA)
        colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("WY", i, i-1, "-WY", i-j+2, i-j+1, sep="")
        pest2A <- (pest2[paste("WY", i, i-1, sep="")] - pest2[paste("WY", i-j+2, i-j+1, sep="")])
        names(pest2A) <- paste("WY", i, i-1, "-WY", i-j+2, i-j+1, sep="")
        pest2 <- append(pest2, pest2A)

        mcmcA <- (mcmc[, paste("WZ", i, i-1, sep="")] - mcmc[, paste("WZ", i-j+2, i-j+1, sep="")])
        mcmc <- cbind(mcmc, mcmcA)
        colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("WZ", i, i-1, "-WZ", i-j+2, i-j+1, sep="")
        pest2A <- (pest2[paste("WZ", i, i-1, sep="")] - pest2[paste("WZ", i-j+2, i-j+1, sep="")])
        names(pest2A) <- paste("WZ", i, i-1, "-WZ", i-j+2, i-j+1, sep="")
        pest2 <- append(pest2, pest2A)
      }  # end (if W != "NULL")
    } # end (for i)
  } # end (for j)
  ## ----- end (Lag = 1) ----- ##


  ## -- Lag = 2 -- ##
  if (lag == 2 & no.waves > 3) {
    for (j in 4:no.waves) {
      for (i in j:no.waves) {
        mcmcA <- (mcmc[, paste("XX", i, i-2, sep="")] - mcmc[, paste("XX", i-j+3, i-j+1, sep="")])
        mcmc <- cbind(mcmc, mcmcA)
        colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("XX", i, i-2, "-XX", i-j+3, i-j+1, sep="")
        pest2A <- (pest2[paste("XX", i, i-2, sep="")] - pest2[paste("XX", i-j+3, i-j+1, sep="")])
        names(pest2A) <- paste("XX", i, i-2, "-XX", i-j+3, i-j+1, sep="")
        pest2 <- append(pest2, pest2A)

        mcmcA <- (mcmc[, paste("YY", i, i-2, sep="")] - mcmc[, paste("YY", i-j+3, i-j+1, sep="")])
        mcmc <- cbind(mcmc, mcmcA)
        colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("YY", i, i-2, "-YY", i-j+3, i-j+1, sep="")
        pest2A <- (pest2[paste("YY", i, i-2, sep="")] - pest2[paste("YY", i-j+3, i-j+1, sep="")])
        names(pest2A) <- paste("YY", i, i-2, "-YY", i-j+3, i-j+1, sep="")
        pest2 <- append(pest2, pest2A)

        mcmcA <- (mcmc[, paste("XY", i, i-2, sep="")] - mcmc[, paste("XY", i-j+3, i-j+1, sep="")])
        mcmc <- cbind(mcmc, mcmcA)
        colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("XY", i, i-2, "-XY", i-j+3, i-j+1, sep="")
        pest2A <- (pest2[paste("XY", i, i-2, sep="")] - pest2[paste("XY", i-j+3, i-j+1, sep="")])
        names(pest2A) <- paste("XY", i, i-2, "-XY", i-j+3, i-j+1, sep="")
        pest2 <- append(pest2, pest2A)

        mcmcA <- (mcmc[, paste("YX", i, i-2, sep="")] - mcmc[, paste("YX", i-j+3, i-j+1, sep="")])
        mcmc <- cbind(mcmc, mcmcA)
        colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("YX", i, i-2, "-YX", i-j+3, i-j+1, sep="")
        pest2A <- (pest2[paste("YX", i, i-2, sep="")] - pest2[paste("YX", i-j+3, i-j+1, sep="")])
        names(pest2A) <- paste("YX", i, i-2, "-YX", i-j+3, i-j+1, sep="")
        pest2 <- append(pest2, pest2A)

        if (Z != "NULL") {
          mcmcA <- (mcmc[, paste("ZZ", i, i-2, sep="")] - mcmc[, paste("ZZ", i-j+3, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("ZZ", i, i-2, "-ZZ", i-j+3, i-j+1, sep="")
          pest2A <- (pest2[paste("ZZ", i, i-2, sep="")] - pest2[paste("ZZ", i-j+3, i-j+1, sep="")])
          names(pest2A) <- paste("ZZ", i, i-2, "-ZZ", i-j+3, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("XZ", i, i-2, sep="")] - mcmc[, paste("XZ", i-j+3, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("XZ", i, i-2, "-XZ", i-j+3, i-j+1, sep="")
          pest2A <- (pest2[paste("XZ", i, i-2, sep="")] - pest2[paste("XZ", i-j+3, i-j+1, sep="")])
          names(pest2A) <- paste("XZ", i, i-2, "-XZ", i-j+3, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("YZ", i, i-2, sep="")] - mcmc[, paste("YZ", i-j+3, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("YZ", i, i-2, "-YZ", i-j+3, i-j+1, sep="")
          pest2A <- (pest2[paste("YZ", i, i-2, sep="")] - pest2[paste("YZ", i-j+3, i-j+1, sep="")])
          names(pest2A) <- paste("YZ", i, i-2, "-YZ", i-j+3, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("ZX", i, i-2, sep="")] - mcmc[, paste("ZX", i-j+3, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("ZX", i, i-2, "-ZX", i-j+3, i-j+1, sep="")
          pest2A <- (pest2[paste("ZX", i, i-2, sep="")] - pest2[paste("ZX", i-j+3, i-j+1, sep="")])
          names(pest2A) <- paste("ZX", i, i-2, "-ZX", i-j+3, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("ZY", i, i-2, sep="")] - mcmc[, paste("ZY", i-j+3, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("ZY", i, i-2, "-ZY", i-j+3, i-j+1, sep="")
          pest2A <- (pest2[paste("ZY", i, i-2, sep="")] - pest2[paste("ZY", i-j+3, i-j+1, sep="")])
          names(pest2A) <- paste("ZY", i, i-2, "-ZY", i-j+3, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)
        } # end (if Z != "NULL")

        if (W != "NULL") {
          mcmcA <- (mcmc[, paste("WW", i, i-2, sep="")] - mcmc[, paste("WW", i-j+3, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("WW", i, i-2, "-WW", i-j+3, i-j+1, sep="")
          pest2A <- (pest2[paste("WW", i, i-2, sep="")] - pest2[paste("WW", i-j+3, i-j+1, sep="")])
          names(pest2A) <- paste("WW", i, i-2, "-ZZ", i-j+3, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("XW", i, i-2, sep="")] - mcmc[, paste("XW", i-j+3, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("XW", i, i-2, "-XW", i-j+3, i-j+1, sep="")
          pest2A <- (pest2[paste("XW", i, i-2, sep="")] - pest2[paste("XW", i-j+3, i-j+1, sep="")])
          names(pest2A) <- paste("XW", i, i-2, "-XW", i-j+3, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("YW", i, i-2, sep="")] - mcmc[, paste("YW", i-j+3, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("YW", i, i-2, "-YW", i-j+3, i-j+1, sep="")
          pest2A <- (pest2[paste("YW", i, i-2, sep="")] - pest2[paste("YW", i-j+3, i-j+1, sep="")])
          names(pest2A) <- paste("YW", i, i-2, "-YW", i-j+3, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("ZW", i, i-2, sep="")] - mcmc[, paste("ZW", i-j+3, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("ZW", i, i-2, "-ZW", i-j+3, i-j+1, sep="")
          pest2A <- (pest2[paste("ZW", i, i-2, sep="")] - pest2[paste("ZW", i-j+3, i-j+1, sep="")])
          names(pest2A) <- paste("ZW", i, i-2, "-ZW", i-j+3, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("WX", i, i-2, sep="")] - mcmc[, paste("WX", i-j+3, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("WX", i, i-2, "-WX", i-j+3, i-j+1, sep="")
          pest2A <- (pest2[paste("WX", i, i-2, sep="")] - pest2[paste("WX", i-j+3, i-j+1, sep="")])
          names(pest2A) <- paste("WX", i, i-2, "-WX", i-j+3, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("WY", i, i-2, sep="")] - mcmc[, paste("WY", i-j+3, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("WY", i, i-2, "-WY", i-j+3, i-j+1, sep="")
          pest2A <- (pest2[paste("WY", i, i-2, sep="")] - pest2[paste("WY", i-j+3, i-j+1, sep="")])
          names(pest2A) <- paste("WY", i, i-2, "-WY", i-j+3, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("WZ", i, i-2, sep="")] - mcmc[, paste("WZ", i-j+3, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("WZ", i, i-2, "-WZ", i-j+3, i-j+1, sep="")
          pest2A <- (pest2[paste("WZ", i, i-2, sep="")] - pest2[paste("WZ", i-j+3, i-j+1, sep="")])
          names(pest2A) <- paste("WZ", i, i-2, "-WZ", i-j+3, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)
        }  # end (if W != "NULL")
      } # end (for i)
    } # end (for j)
  }  # end (lag == 2)


  ## -- Lag = 3 -- ##
  if (lag == 3 & no.waves > 4) {
    for (j in 5:no.waves) {
      for (i in j:no.waves) {
        mcmcA <- (mcmc[, paste("XX", i, i-3, sep="")] - mcmc[, paste("XX", i-j+4, i-j+1, sep="")])
        mcmc <- cbind(mcmc, mcmcA)
        colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("XX", i, i-3, "-XX", i-j+4, i-j+1, sep="")
        pest2A <- (pest2[paste("XX", i, i-3, sep="")] - pest2[paste("XX", i-j+4, i-j+1, sep="")])
        names(pest2A) <- paste("XX", i, i-3, "-XX", i-j+4, i-j+1, sep="")
        pest2 <- append(pest2, pest2A)

        mcmcA <- (mcmc[, paste("YY", i, i-3, sep="")] - mcmc[, paste("YY", i-j+4, i-j+1, sep="")])
        mcmc <- cbind(mcmc, mcmcA)
        colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("YY", i, i-3, "-YY", i-j+4, i-j+1, sep="")
        pest2A <- (pest2[paste("YY", i, i-3, sep="")] - pest2[paste("YY", i-j+4, i-j+1, sep="")])
        names(pest2A) <- paste("YY", i, i-3, "-YY", i-j+4, i-j+1, sep="")
        pest2 <- append(pest2, pest2A)

        mcmcA <- (mcmc[, paste("XY", i, i-3, sep="")] - mcmc[, paste("XY", i-j+4, i-j+1, sep="")])
        mcmc <- cbind(mcmc, mcmcA)
        colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("XY", i, i-3, "-XY", i-j+4, i-j+1, sep="")
        pest2A <- (pest2[paste("XY", i, i-3, sep="")] - pest2[paste("XY", i-j+4, i-j+1, sep="")])
        names(pest2A) <- paste("XY", i, i-3, "-XY", i-j+4, i-j+1, sep="")
        pest2 <- append(pest2, pest2A)

        mcmcA <- (mcmc[, paste("YX", i, i-3, sep="")] - mcmc[, paste("YX", i-j+4, i-j+1, sep="")])
        mcmc <- cbind(mcmc, mcmcA)
        colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("YX", i, i-3, "-YX", i-j+4, i-j+1, sep="")
        pest2A <- (pest2[paste("YX", i, i-3, sep="")] - pest2[paste("YX", i-j+4, i-j+1, sep="")])
        names(pest2A) <- paste("YX", i, i-3, "-YX", i-j+4, i-j+1, sep="")
        pest2 <- append(pest2, pest2A)

        if (Z != "NULL") {
          mcmcA <- (mcmc[, paste("ZZ", i, i-3, sep="")] - mcmc[, paste("ZZ", i-j+4, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("ZZ", i, i-3, "-ZZ", i-j+4, i-j+1, sep="")
          pest2A <- (pest2[paste("ZZ", i, i-3, sep="")] - pest2[paste("ZZ", i-j+4, i-j+1, sep="")])
          names(pest2A) <- paste("ZZ", i, i-3, "-ZZ", i-j+4, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("XZ", i, i-3, sep="")] - mcmc[, paste("XZ", i-j+4, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("XZ", i, i-3, "-XZ", i-j+4, i-j+1, sep="")
          pest2A <- (pest2[paste("XZ", i, i-3, sep="")] - pest2[paste("XZ", i-j+4, i-j+1, sep="")])
          names(pest2A) <- paste("XZ", i, i-3, "-XZ", i-j+4, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("YZ", i, i-3, sep="")] - mcmc[, paste("YZ", i-j+4, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("YZ", i, i-3, "-YZ", i-j+4, i-j+1, sep="")
          pest2A <- (pest2[paste("YZ", i, i-3, sep="")] - pest2[paste("YZ", i-j+4, i-j+1, sep="")])
          names(pest2A) <- paste("YZ", i, i-3, "-YZ", i-j+4, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("ZX", i, i-3, sep="")] - mcmc[, paste("ZX", i-j+4, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("ZX", i, i-3, "-ZX", i-j+4, i-j+1, sep="")
          pest2A <- (pest2[paste("ZX", i, i-3, sep="")] - pest2[paste("ZX", i-j+4, i-j+1, sep="")])
          names(pest2A) <- paste("ZX", i, i-3, "-ZX", i-j+4, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("ZY", i, i-3, sep="")] - mcmc[, paste("ZY", i-j+4, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("ZY", i, i-3, "-ZY", i-j+4, i-j+1, sep="")
          pest2A <- (pest2[paste("ZY", i, i-3, sep="")] - pest2[paste("ZY", i-j+4, i-j+1, sep="")])
          names(pest2A) <- paste("ZY", i, i-3, "-ZY", i-j+4, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)
        } # end (if Z != "NULL")

        if (W != "NULL") {
          mcmcA <- (mcmc[, paste("WW", i, i-3, sep="")] - mcmc[, paste("WW", i-j+4, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("WW", i, i-3, "-WW", i-j+4, i-j+1, sep="")
          pest2A <- (pest2[paste("WW", i, i-3, sep="")] - pest2[paste("WW", i-j+4, i-j+1, sep="")])
          names(pest2A) <- paste("WW", i, i-3, "-ZZ", i-j+4, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("XW", i, i-3, sep="")] - mcmc[, paste("XW", i-j+4, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("XW", i, i-3, "-XW", i-j+4, i-j+1, sep="")
          pest2A <- (pest2[paste("XW", i, i-3, sep="")] - pest2[paste("XW", i-j+4, i-j+1, sep="")])
          names(pest2A) <- paste("XW", i, i-3, "-XW", i-j+4, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("YW", i, i-3, sep="")] - mcmc[, paste("YW", i-j+4, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("YW", i, i-3, "-YW", i-j+4, i-j+1, sep="")
          pest2A <- (pest2[paste("YW", i, i-3, sep="")] - pest2[paste("YW", i-j+4, i-j+1, sep="")])
          names(pest2A) <- paste("YW", i, i-3, "-YW", i-j+4, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("ZW", i, i-3, sep="")] - mcmc[, paste("ZW", i-j+4, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("ZW", i, i-3, "-ZW", i-j+4, i-j+1, sep="")
          pest2A <- (pest2[paste("ZW", i, i-3, sep="")] - pest2[paste("ZW", i-j+4, i-j+1, sep="")])
          names(pest2A) <- paste("ZW", i, i-3, "-ZW", i-j+4, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("WX", i, i-3, sep="")] - mcmc[, paste("WX", i-j+4, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("WX", i, i-3, "-WX", i-j+4, i-j+1, sep="")
          pest2A <- (pest2[paste("WX", i, i-3, sep="")] - pest2[paste("WX", i-j+4, i-j+1, sep="")])
          names(pest2A) <- paste("WX", i, i-3, "-WX", i-j+4, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("WY", i, i-3, sep="")] - mcmc[, paste("WY", i-j+4, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("WY", i, i-3, "-WY", i-j+4, i-j+1, sep="")
          pest2A <- (pest2[paste("WY", i, i-3, sep="")] - pest2[paste("WY", i-j+4, i-j+1, sep="")])
          names(pest2A) <- paste("WY", i, i-3, "-WY", i-j+4, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("WZ", i, i-3, sep="")] - mcmc[, paste("WZ", i-j+4, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("WZ", i, i-3, "-WZ", i-j+4, i-j+1, sep="")
          pest2A <- (pest2[paste("WZ", i, i-3, sep="")] - pest2[paste("WZ", i-j+4, i-j+1, sep="")])
          names(pest2A) <- paste("WZ", i, i-3, "-WZ", i-j+4, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)
        }  # end (if W != "NULL")
      } # end (for i)
    } # end (for j)
  } # end (lag == 3)


  ## Lag = 4 ##
  if (lag == 4 & no.waves > 5) {
    for (j in 6:no.waves) {
      for (i in j:no.waves) {
        mcmcA <- (mcmc[, paste("XX", i, i-4, sep="")] - mcmc[, paste("XX", i-j+5, i-j+1, sep="")])
        mcmc <- cbind(mcmc, mcmcA)
        colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("XX", i, i-4, "-XX", i-j+5, i-j+1, sep="")
        pest2A <- (pest2[paste("XX", i, i-4, sep="")] - pest2[paste("XX", i-j+5, i-j+1, sep="")])
        names(pest2A) <- paste("XX", i, i-4, "-XX", i-j+5, i-j+1, sep="")
        pest2 <- append(pest2, pest2A)

        mcmcA <- (mcmc[, paste("YY", i, i-4, sep="")] - mcmc[, paste("YY", i-j+5, i-j+1, sep="")])
        mcmc <- cbind(mcmc, mcmcA)
        colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("YY", i, i-4, "-YY", i-j+5, i-j+1, sep="")
        pest2A <- (pest2[paste("YY", i, i-4, sep="")] - pest2[paste("YY", i-j+5, i-j+1, sep="")])
        names(pest2A) <- paste("YY", i, i-4, "-YY", i-j+5, i-j+1, sep="")
        pest2 <- append(pest2, pest2A)

        mcmcA <- (mcmc[, paste("XY", i, i-4, sep="")] - mcmc[, paste("XY", i-j+5, i-j+1, sep="")])
        mcmc <- cbind(mcmc, mcmcA)
        colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("XY", i, i-4, "-XY", i-j+5, i-j+1, sep="")
        pest2A <- (pest2[paste("XY", i, i-4, sep="")] - pest2[paste("XY", i-j+5, i-j+1, sep="")])
        names(pest2A) <- paste("XY", i, i-4, "-XY", i-j+5, i-j+1, sep="")
        pest2 <- append(pest2, pest2A)

        mcmcA <- (mcmc[, paste("YX", i, i-4, sep="")] - mcmc[, paste("YX", i-j+5, i-j+1, sep="")])
        mcmc <- cbind(mcmc, mcmcA)
        colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("YX", i, i-4, "-YX", i-j+5, i-j+1, sep="")
        pest2A <- (pest2[paste("YX", i, i-4, sep="")] - pest2[paste("YX", i-j+5, i-j+1, sep="")])
        names(pest2A) <- paste("YX", i, i-4, "-YX", i-j+5, i-j+1, sep="")
        pest2 <- append(pest2, pest2A)

        if (Z != "NULL") {
          mcmcA <- (mcmc[, paste("ZZ", i, i-4, sep="")] - mcmc[, paste("ZZ", i-j+5, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("ZZ", i, i-4, "-ZZ", i-j+5, i-j+1, sep="")
          pest2A <- (pest2[paste("ZZ", i, i-4, sep="")] - pest2[paste("ZZ", i-j+5, i-j+1, sep="")])
          names(pest2A) <- paste("ZZ", i, i-4, "-ZZ", i-j+5, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("XZ", i, i-4, sep="")] - mcmc[, paste("XZ", i-j+5, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("XZ", i, i-4, "-XZ", i-j+5, i-j+1, sep="")
          pest2A <- (pest2[paste("XZ", i, i-4, sep="")] - pest2[paste("XZ", i-j+5, i-j+1, sep="")])
          names(pest2A) <- paste("XZ", i, i-4, "-XZ", i-j+5, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("YZ", i, i-4, sep="")] - mcmc[, paste("YZ", i-j+5, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("YZ", i, i-4, "-YZ", i-j+5, i-j+1, sep="")
          pest2A <- (pest2[paste("YZ", i, i-4, sep="")] - pest2[paste("YZ", i-j+5, i-j+1, sep="")])
          names(pest2A) <- paste("YZ", i, i-4, "-YZ", i-j+5, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("ZX", i, i-4, sep="")] - mcmc[, paste("ZX", i-j+5, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("ZX", i, i-4, "-ZX", i-j+5, i-j+1, sep="")
          pest2A <- (pest2[paste("ZX", i, i-4, sep="")] - pest2[paste("ZX", i-j+5, i-j+1, sep="")])
          names(pest2A) <- paste("ZX", i, i-4, "-ZX", i-j+5, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("ZY", i, i-4, sep="")] - mcmc[, paste("ZY", i-j+5, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("ZY", i, i-4, "-ZY", i-j+5, i-j+1, sep="")
          pest2A <- (pest2[paste("ZY", i, i-4, sep="")] - pest2[paste("ZY", i-j+5, i-j+1, sep="")])
          names(pest2A) <- paste("ZY", i, i-4, "-ZY", i-j+5, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)
        } # end (if Z != "NULL")

        if (W != "NULL") {
          mcmcA <- (mcmc[, paste("WW", i, i-4, sep="")] - mcmc[, paste("WW", i-j+5, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("WW", i, i-4, "-WW", i-j+5, i-j+1, sep="")
          pest2A <- (pest2[paste("WW", i, i-4, sep="")] - pest2[paste("WW", i-j+5, i-j+1, sep="")])
          names(pest2A) <- paste("WW", i, i-4, "-ZZ", i-j+5, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("XW", i, i-4, sep="")] - mcmc[, paste("XW", i-j+5, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("XW", i, i-4, "-XW", i-j+5, i-j+1, sep="")
          pest2A <- (pest2[paste("XW", i, i-4, sep="")] - pest2[paste("XW", i-j+5, i-j+1, sep="")])
          names(pest2A) <- paste("XW", i, i-4, "-XW", i-j+5, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("YW", i, i-4, sep="")] - mcmc[, paste("YW", i-j+5, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("YW", i, i-4, "-YW", i-j+5, i-j+1, sep="")
          pest2A <- (pest2[paste("YW", i, i-4, sep="")] - pest2[paste("YW", i-j+5, i-j+1, sep="")])
          names(pest2A) <- paste("YW", i, i-4, "-YW", i-j+5, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("ZW", i, i-4, sep="")] - mcmc[, paste("ZW", i-j+5, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("ZW", i, i-4, "-ZW", i-j+5, i-j+1, sep="")
          pest2A <- (pest2[paste("ZW", i, i-4, sep="")] - pest2[paste("ZW", i-j+5, i-j+1, sep="")])
          names(pest2A) <- paste("ZW", i, i-4, "-ZW", i-j+5, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("WX", i, i-4, sep="")] - mcmc[, paste("WX", i-j+5, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("WX", i, i-4, "-WX", i-j+5, i-j+1, sep="")
          pest2A <- (pest2[paste("WX", i, i-4, sep="")] - pest2[paste("WX", i-j+5, i-j+1, sep="")])
          names(pest2A) <- paste("WX", i, i-4, "-WX", i-j+5, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("WY", i, i-4, sep="")] - mcmc[, paste("WY", i-j+5, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("WY", i, i-4, "-WY", i-j+5, i-j+1, sep="")
          pest2A <- (pest2[paste("WY", i, i-4, sep="")] - pest2[paste("WY", i-j+5, i-j+1, sep="")])
          names(pest2A) <- paste("WY", i, i-4, "-WY", i-j+5, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)

          mcmcA <- (mcmc[, paste("WZ", i, i-4, sep="")] - mcmc[, paste("WZ", i-j+5, i-j+1, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("WZ", i, i-4, "-WZ", i-j+5, i-j+1, sep="")
          pest2A <- (pest2[paste("WZ", i, i-4, sep="")] - pest2[paste("WZ", i-j+5, i-j+1, sep="")])
          names(pest2A) <- paste("WZ", i, i-4, "-WZ", i-j+5, i-j+1, sep="")
          pest2 <- append(pest2, pest2A)
        } # end (if W != "NULL")
      } # end (for i)
    } # end (for j)
  } # end (lag == 4)


  ## -- Differences in covariance -- ##
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
  ## ----- end (Difference in covariance) ----- ##

  ## -- Differences in Grand Means -- ##
  for (j in 2:no.waves) {
    for (i in j:no.waves) {
      mcmcA <- (mcmc[, paste("M", X, i, sep="")] - mcmc[, paste("M", X, i-1-j+2, sep="")])
      mcmc <- cbind(mcmc, mcmcA)
      colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("M", X, i, "-M", X, i-1-j+2, sep="")
      pest2A <- (pest2[paste("M", X, i, sep="")] - pest2[paste("M", X, i-1-j+2, sep="")])
      names(pest2A) <- paste("M", X, i, "-M", X, i-1-j+2, sep="")
      pest2 <- append(pest2, pest2A)

      mcmcA <- (mcmc[, paste("M", Y, i, sep="")] - mcmc[, paste("M", Y, i-1-j+2, sep="")])
      mcmc <- cbind(mcmc, mcmcA)
      colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("M", Y, i, "-M", Y, i-1-j+2, sep="")
      pest2A <- (pest2[paste("M", Y, i, sep="")] - pest2[paste("M", Y, i-1-j+2, sep="")])
      names(pest2A) <- paste("M", Y, i, "-M", Y, i-1-j+2, sep="")
      pest2 <- append(pest2, pest2A)

      if (Z != "NULL") {
        mcmcA <- (mcmc[, paste("M", Z, i, sep="")] - mcmc[, paste("M", Z, i-1-j+2, sep="")])
        mcmc <- cbind(mcmc, mcmcA)
        colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("M", Z, i, "-M", Z, i-1-j+2, sep="")
        pest2A <- (pest2[paste("M", Z, i, sep="")] - pest2[paste("M", Z, i-1-j+2, sep="")])
        names(pest2A) <- paste("M", Z, i, "-M", Z, i-1-j+2, sep="")
        pest2 <- append(pest2, pest2A)
      } # end (if Z != "NULL")

      if (W != "NULL") {
        mcmcA <- (mcmc[, paste("M", W, i, sep="")] - mcmc[, paste("M", W, i-1-j+2, sep="")])
        mcmc <- cbind(mcmc, mcmcA)
        colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("M", W, i, "-M", W, i-1-j+2, sep="")
        pest2A <- (pest2[paste("M", W, i, sep="")] - pest2[paste("M", W, i-1-j+2, sep="")])
        names(pest2A) <- paste("M", W, i, "-M", W, i-1-j+2, sep="")
        pest2 <- append(pest2, pest2A)
      } # End (if W != "NULL")
    } # end (for i)
  } # end (for j)
  ## ----- end (Differences in Grand Mean) ----- ##


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


  ## ---- List and Delete - Path Coefficients ---- ##

  no.path = no.waves - 1
  MIset <- no.waves - 3
  no.compare = (no.path - 1)*(no.path)/2
  no.compare.M = (no.waves - 1)*(no.waves)/2
  MIset.M <- no.waves - 2

  cat(rep("\n",3), "## ===== Identification of invariant paths (Lag = 1 wave) ===== ##")

  LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X", b="X", LDlag=1)  ## List and Delete - Path XX ##
  LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Y", b="Y", LDlag=1)  ## List and Delete - Path YY ##
  if (Z != "NULL") {
    LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Z", b="Z", LDlag=1)  ## List and Delete - Path ZZ ##
  } # end (if Z)
  if (W != "NULL") {
    LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="W", b="W", LDlag=1)  ## List and Delete - Path WW ##
  } # end (if W)

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

    LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X", b="X", LDlag=2)  ## List and Delete - Path XX ##
    LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Y", b="Y", LDlag=2)  ## List and Delete - Path YY ##
    if (Z != "NULL") {
      LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Z", b="Z", LDlag=2)  ## List and Delete - Path ZZ ##
    } # end (if Z != "NULL")
    if (W != "NULL") {
      LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="W", b="W", LDlag=2)  ## List and Delete - Path WW ##
    } # end (if W != "NULL")

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

    LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X", b="X", LDlag=3)  ## List and Delete - Path XX ##
    LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Y", b="Y", LDlag=3)  ## List and Delete - Path YY ##
    if (Z != "NULL") {
      LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Z", b="Z", LDlag=3)  ## List and Delete - Path ZZ ##
    } # end (if Z != "NULL")
    if (W != "NULL") {
      LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="W", b="W", LDlag=3)  ## List and Delete - Path WW ##
    } # end (if W != "NULL")

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

    LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X", b="X", LDlag=4)  ## List and Delete - Path XX ##
    LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Y", b="Y", LDlag=4)  ## List and Delete - Path YY ##
    if (Z != "NULL") {
      LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Z", b="Z", LDlag=4)  ## List and Delete - Path ZZ ##
    } # end (if Z != "NULL")
    if (W != "NULL") {
      LandD_Path(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="W", b="W", LDlag=4)  ## List and Delete - Path WW ##
    } # end (if W != "NULL")

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

  ## -------------------------------------------------------- ##


  ## ---- List and Delete - Residual covariance eXY ---- ##

  # -- Reset MISet and no.compare for residual covariance -- #
  no.path = no.waves - 1
  MIset <- no.waves - 3
  no.compare = (no.path - 1)*(no.path)/2

  cat(rep("\n",5), "## ===== Identification of invariant residual covariance ===== ##")

  LandD_COV(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X", b="Y")  ## List and Delete - eXY ##

  if (Z != "NULL") {
    LandD_COV(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X", b="Z")  ## List and Delete - eXZ ##
    LandD_COV(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Y", b="Z")  ## List and Delete - eYZ ##
  } # end (if Z != "NULL")

  if (W != "NULL") {
    LandD_COV(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X", b="W")  ## List and Delete - eXW ##
    LandD_COV(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Y", b="W")  ## List and Delete - eYW ##
    LandD_COV(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Z", b="W")  ## List and Delete - eZW ##
  } # end (if W != "NULL")

  ## -------------------------------------------------------- ##


  ## ---- List and Delete - Grand Mean ---- ##

  # -- Reset MISet and no.compare for grand mean  --#
  no.compare = (no.waves - 1)*(no.waves)/2
  MIset <- no.waves - 2

  cat(rep("\n",5), "## ===== Identification of invariant means ===== ##")

  LandD_MEAN(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X")  ## List and Delete - Grand mean of X ##
  LandD_MEAN(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Y")  ## List and Delete - Grand mean of Y ##

  if (Z != "NULL") {
    LandD_MEAN(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Z")  ## List and Delete - Grand mean of Z ##
  } # end (if Z != "NULL")

  if (W != "NULL") {
    LandD_MEAN(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="W")  ## List and Delete - Grand mean of W ##
  } # end (if W != "NULL")

  ## -------------------------------------------------------- ##


  ## -- Differences in Indicator Residual Variance -- ##
  if (isFALSE(varI.eq)) {
    for (j in 2:no.waves) {
      for (i in j:no.waves) {
        mcmcA <- (mcmc[, paste("varI", X, i, sep="")] - mcmc[, paste("varI", X, i-1-j+2, sep="")])
        mcmc <- cbind(mcmc, mcmcA)
        colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("varI", X, i, "-varI", X, i-1-j+2, sep="")
        pest2A <- (pest2[paste("varI", X, i, sep="")] - pest2[paste("varI", X, i-1-j+2, sep="")])
        names(pest2A) <- paste("varI", X, i, "-varI", X, i-1-j+2, sep="")
        pest2 <- append(pest2, pest2A)

        mcmcA <- (mcmc[, paste("varI", Y, i, sep="")] - mcmc[, paste("varI", Y, i-1-j+2, sep="")])
        mcmc <- cbind(mcmc, mcmcA)
        colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("varI", Y, i, "-varI", Y, i-1-j+2, sep="")
        pest2A <- (pest2[paste("varI", Y, i, sep="")] - pest2[paste("varI", Y, i-1-j+2, sep="")])
        names(pest2A) <- paste("varI", Y, i, "-varI", Y, i-1-j+2, sep="")
        pest2 <- append(pest2, pest2A)

        if (Z != "NULL") {
          mcmcA <- (mcmc[, paste("varI", Z, i, sep="")] - mcmc[, paste("varI", Z, i-1-j+2, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("varI", Z, i, "-varI", Z, i-1-j+2, sep="")
          pest2A <- (pest2[paste("varI", Z, i, sep="")] - pest2[paste("varI", Z, i-1-j+2, sep="")])
          names(pest2A) <- paste("varI", Z, i, "-varI", Z, i-1-j+2, sep="")
          pest2 <- append(pest2, pest2A)
        } # end (if Z != "NULL")

        if (W != "NULL") {
          mcmcA <- (mcmc[, paste("varI", W, i, sep="")] - mcmc[, paste("varI", W, i-1-j+2, sep="")])
          mcmc <- cbind(mcmc, mcmcA)
          colnames(mcmc)[colnames(mcmc) == "mcmcA"] = paste("varI", W, i, "-varI", W, i-1-j+2, sep="")
          pest2A <- (pest2[paste("varI", W, i, sep="")] - pest2[paste("varI", W, i-1-j+2, sep="")])
          names(pest2A) <- paste("varI", W, i, "-varI", W, i-1-j+2, sep="")
          pest2 <- append(pest2, pest2A)
        } # End (if W != "NULL")
      } # end (for i)
    } # end (for j)

    # -- Reset MISet and no.compare for grand mean  --#
    no.compare = (no.waves - 1)*(no.waves)/2
    MIset <- no.waves - 2

    cat(rep("\n",5), "## ===== Identification of invariant indicator variance ===== ##")

    LandD_varI(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X")  ## List and Delete - Indicator Variance of X ##
    LandD_varI(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Y")  ## List and Delete - Indicator Variance of Y ##

    if (Z != "NULL") {
      LandD_varI(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="Z")  ## List and Delete - Indicator Variance of Z ##
    } # end (if Z != "NULL")

    if (W != "NULL") {
      LandD_varI(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="W")  ## List and Delete - Indicator Variance of W ##
    } # end (if W != "NULL")
  } # end (if varI.eq == FALSE)

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
    p.TparEst <- paste("  ", a, b, i+LDlag, i, ":  ", Clhs, " ~ ", Crhs, " = ", format(round(TparEst["est"], digits=4), nsmall=4, scientific=FALSE),
                       ", p-value = ", format(round(TparEst["pvalue"], digits=4), nsmall=4, scientific=FALSE), sep="")
    cat("\n", p.TparEst)
    SumEst <- SumEst + TparEst["est"]
  }
  MeanEst <- SumEst/(no.path-LDlag+1)
  p.MeanEst <- paste("  Mean across paths = ", format(round(MeanEst, digits=4), nsmall=4, scientific=FALSE), sep="")
  cat("\n", p.MeanEst)

  # -- Save pairs of non-invariant paths -- #
  NI.path <- t(matrix(0,ncol=no.compare, nrow=2))
  k = 1
  for (j in 3:(no.waves-LDlag+1)) {
    for (i in j:(no.waves-LDlag+1)) {
      Clhs <- paste(a, b, (i+LDlag-1), i-1, "-", a, b, i-j+LDlag+1, i-j+1, sep="")
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
          mm.p[mm.p == i] <- paste(a, b, i+LDlag, i, sep="")
        }
        cat("\n", "    ", mm.p[ii,])
      } # end (if count4)
    } # end (for ii)
  } # end (for k)
}  # end (Function LandD_Path)
## --------------------------------------------------------- ##




## ----- Sub-Function List and Delete for Residual Covariance ----- ##
LandD_COV <- function(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X", b="Y") {

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
  p.MeanEst <- paste("   Mean across waves from T2 = ", format(round(MeanEst, digits=4), nsmall=4, scientific=FALSE), sep="")
  cat("\n", p.MeanEst)

  # -- Save pairs of non-invariant residual covariances -- #
  NI.path <- t(matrix(0,ncol=no.compare, nrow=2))
  k = 1
  for (j in 3:no.waves) {
    for (i in j:no.waves) {
      Clhs <- paste("e", a, b, i, "-e", a, b, i-j+2, sep="")
      if (pest2[Clhs,3] < p) {
        NI.path[k, 1] <- i-1
        NI.path[k, 2] <- i-j+1
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
        }
      }
    }
    for (ii in 1:NIset) {
      count4 <- sum(mm[ii,])
      if (count4 > 0) {
        mm.p <- mm
        for (i in 1:no.path) {
          mm.p[mm.p == i] <- paste("e", a, b, i+1, sep="")
        }
        cat("\n", "    ", mm.p[ii,])
      }
    }
  }
}  # end (Function LandD_COV)
## ---------------------------------------------------------------- ##





## ----- Sub-Function List and Delete for Grand Mean ----- ##
LandD_MEAN <- function(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X") {

  SumEst <- 0
  cat(rep("\n",3), paste("# -- Grand mean of ", a, " -- #", sep=""))
  aa <- get(a)
  for (i in 1: no.waves) {
    Clhs <- paste(aa, i, sep="")
    TparEst <- parEst[parEst["lhs"] == Clhs & parEst["op"] == "~1",]
    p.TparEst <- paste("  M", aa, i, ":  Mean of ", Clhs, " = ", format(round(TparEst["est"], digits=4), nsmall=4, scientific=FALSE),
                       ", p-value = ", format(round(TparEst["pvalue"], digits=4), nsmall=4, scientific=FALSE), sep="")
    cat("\n", p.TparEst)
    SumEst <- SumEst + TparEst["est"]
  } # end (for i)
  MeanEst <- SumEst/no.waves
  p.MeanEst <- paste("  Mean across waves = ", format(round(MeanEst, digits=4), nsmall=4, scientific=FALSE), sep="")
  cat("\n", p.MeanEst)

  # -- Save pairs of non-invariant Grand means -- #
  NI.path <- t(matrix(0,ncol=no.compare, nrow=2))
  k = 1
  for (j in 2:no.waves) {
    for (i in j:no.waves) {
      Clhs <- paste("M", aa, i, "-M", aa, i-j+1, sep="")
      if (pest2[Clhs,3] < p) {
        NI.path[k, 1] <- i
        NI.path[k, 2] <- i-j+1
        k <- k + 1
      } # end (if pest2)
    }
  }
  Count.NI.path = k - 1

  # -- Select sets of invariant grand means and print -- #
  cat(rep("\n",2), paste("# -- Sets of invariant Grand means ", aa, " -- #", sep=""))

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
        for (i in 0:no.waves) {
          mm.p[mm.p == i] <- paste("M", aa, i, sep="")
        } # end (for i)
        cat("\n", "    ", mm.p[ii,])
      } # end (if count4)
    } # end (for ii)
  } # end (for k)
}  # end (Function LandD_MEAN)
## ---------------------------------------------------------------- ##




## ----- Sub-Function List and Delete for Indicator Variance ----- ##
LandD_varI <- function(parEst, pest2, no.path, MIset, no.compare, no.waves, p, X, Y, Z, W, a="X") {

  SumEst <- 0
  cat(rep("\n",3), paste("# -- Indicator Variance of ", a, " -- #", sep=""))
  aa <- get(a)
  for (i in 1: no.waves) {
    Clhs <- paste(aa, i, sep="")
    TparEst <- parEst[parEst["lhs"] == Clhs & parEst["rhs"] == Clhs & parEst["op"] == "~~",]
    TparEst <- parEst[parEst["lhs"] == Clhs & parEst["rhs"] == Crhs & parEst["op"] == "~~",]
    p.TparEst <- paste("  varI", aa, i, ":  Indicator Variance of ", Clhs, " = ", format(round(TparEst["est"], digits=4), nsmall=4, scientific=FALSE),
                       ", p-value = ", format(round(TparEst["pvalue"], digits=4), nsmall=4, scientific=FALSE), sep="")
    cat("\n", p.TparEst)
    SumEst <- SumEst + TparEst["est"]
  } # end (for i)
  MeanEst <- SumEst/no.waves
  p.MeanEst <- paste("  Mean across waves = ", format(round(MeanEst, digits=4), nsmall=4, scientific=FALSE), sep="")
  cat("\n", p.MeanEst)

  # -- Save pairs of non-invariant Indicator Variance -- #
  NI.path <- t(matrix(0,ncol=no.compare, nrow=2))
  k = 1
  for (j in 2:no.waves) {
    for (i in j:no.waves) {
      Clhs <- paste("varI", aa, i, "-varI", aa, i-j+1, sep="")
      if (pest2[Clhs,3] < p) {
        NI.path[k, 1] <- i
        NI.path[k, 2] <- i-j+1
        k <- k + 1
      } # end (if pest2)
    }
  }
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
        for (i in 0:no.waves) {
          mm.p[mm.p == i] <- paste("varI", aa, i, sep="")
        } # end (for i)
        cat("\n", "    ", mm.p[ii,])
      } # end (if count4)
    } # end (for ii)
  } # end (for k)
}  # end (Sub-Function LandD_varI)
## ---------------------------------------------------------------- ##





# ==================== Creating Function "CLPM" ==================== #
#' Function Cross-Lagged Panel Model (CLPM)
#'
#' Cross-Lagged Panel Model
#'
#' @param data.source name of data.frame
#' @param no.waves number of waves (minimum = 3, must be grater than lag)
#' @param lag number of waves between two lags (minimum = 1, maximum = 4)
#' @param p critical p-value for pairwise comparisons (default is 0.001)
#' @param X name of variable X.
#' @param Y name of variable Y.
#' @param Z name of variable Z.
#' @param W name of variable W (Z must come before W).
#'
#' @return CLPM outputs.
#' @export
#' @examples
#'
#' ## -- Example -- ##
#'
#' CLPM(data.source="Data_A", 7, 2, X="EXPOSE", Y="INTENS")
#'

CLPM <- function(data.source, no.waves, lag=1, p = 0.001, X, Y, Z="NULL", W = "NULL") {

  ## -- Check inputs -- ##

  if (no.waves < 3) stop("Minimum number of waves is 3")

  if (lag < 1) stop("Minimum number of lag is 1")
  if (lag > 4) stop("Maximum number of lag is 4")

  if ((no.waves - lag) < 1) stop("Number of waves must be greater than lag plus 1")

  if (p > 0.05) stop("p > 0.05 is not recommended")
  if (p < 0.0001) stop("p < 0.0001 is not recommended")

  if (Z == "NULL" & W != "NULL") stop("Z must be defined before W")


  ## ----- Creating Model CLPM ----- ##
  sink('CLPM.txt')  # Start writing script to CLPM.txt

    cat("\n", "## ----- Specify the model (CLPM) ----- ##", "\n")
    cat("\n", "CLPM <- '")

    # -- Constrain Residual Variance of Observed Variables to Zero -- #
    cat(rep("\n",2), "  # -- Constrain residual variance of observed variables to zero -- #")
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

    # -- Create Latent Variables from Observed Variables -- #
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

    # -- Estimate Covariances Between Residuals of Latent Variables -- #
    cat(rep("\n",2), "  # -- Estimate covariances between residuals of latent variables -- #")
    for (i in 2:no.waves) {
      cat("\n", paste("  w", X, i, " ~~ eXY", i, "*w", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  w", X, i, " ~~ eXZ", i, "*w", Z, i, sep=""))
        cat("\n", paste("  w", Y, i, " ~~ eYZ", i, "*w", Z, i, sep=""))
      }  # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  w", X, i, " ~~ eXW", i, "*w", W, i, sep=""))
        cat("\n", paste("  w", Y, i, " ~~ eYW", i, "*w", W, i, sep=""))
        cat("\n", paste("  w", Z, i, " ~~ eZW", i, "*w", W, i, sep=""))
      }  # end (if W)
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

    # -- Estimate (Residual) Variance of Latent Variables -- #
    cat(rep("\n",2), "  # -- Estimate (residual) variance of latent variables -- #")
    for (i in 1:no.waves) {
      cat("\n", paste("  w", X, i, " ~~ ", "w", X, i, sep=""))
      cat("\n", paste("  w", Y, i, " ~~ ", "w", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  w", Z, i, " ~~ ", "w", Z, i, sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  w", W, i, " ~~ ", "w", W, i, sep=""))
      } # end (if W)
    } # end (for i)

    # -- Estimate Grand Means (Intercepts) of Observed Variables -- #
    cat(rep("\n",2), "  # -- Estimate grand means (intercepts) of observed variables -- #")
    for (i in 1:no.waves) {
      cat("\n", paste("  ", X, i, " ~ M", X, i, "*1", sep=""))
      cat("\n", paste("  ", Y, i, " ~ M", Y, i, "*1", sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  ", Z, i, " ~ M", Z, i, "*1", sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  ", W, i, " ~ M", W, i, "*1", sep=""))
      } # (if W)
    } # end (for i)

    # -- Estimate Lagged Effects Between Latent Variables -- #
    cat(rep("\n",2), "  # -- Estimate lagged effects between latent variables (Lag = 1 wave) -- #")
    for (i in 2:no.waves) {
      if (Z == "NULL") {
        cat("\n", paste("  w", X, i, " ~ XX", i,i-1, "*w", X,i-1, " + YX", i,i-1, "*w", Y,i-1, sep=""))
        cat("\n", paste("  w", Y, i, " ~ XY", i,i-1, "*w", X,i-1, " + YY", i,i-1, "*w", Y,i-1, sep=""))
      } else if (W != "NULL") {
        cat("\n", paste("  w", X, i, " ~ XX", i,i-1, "*w", X,i-1, " + YX", i,i-1, "*w", Y,i-1, " + ZX", i,i-1, "*w", Z,i-1, " + WX", i,i-1, "*w", W,i-1, sep=""))
        cat("\n", paste("  w", Y, i, " ~ XY", i,i-1, "*w", X,i-1, " + YY", i,i-1, "*w", Y,i-1, " + ZY", i,i-1, "*w", Z,i-1, " + WY", i,i-1, "*w", W,i-1, sep=""))
        cat("\n", paste("  w", Z, i, " ~ XZ", i,i-1, "*w", X,i-1, " + YZ", i,i-1, "*w", Y,i-1, " + ZZ", i,i-1, "*w", Z,i-1, " + WZ", i,i-1, "*w", W,i-1, sep=""))
        cat("\n", paste("  w", W, i, " ~ XW", i,i-1, "*w", X,i-1, " + YW", i,i-1, "*w", Y,i-1, " + ZW", i,i-1, "*w", Z,i-1, " + WW", i,i-1, "*w", W,i-1, sep=""))
      } else if (Z != "NULL") {
        cat("\n", paste("  w", X, i, " ~ XX", i,i-1, "*w", X,i-1, " + YX", i,i-1, "*w", Y,i-1, " + ZX", i,i-1, "*w", Z,i-1, sep=""))
        cat("\n", paste("  w", Y, i, " ~ XY", i,i-1, "*w", X,i-1, " + YY", i,i-1, "*w", Y,i-1, " + ZY", i,i-1, "*w", Z,i-1, sep=""))
        cat("\n", paste("  w", Z, i, " ~ XZ", i,i-1, "*w", X,i-1, " + YZ", i,i-1, "*w", Y,i-1, " + ZZ", i,i-1, "*w", Z,i-1, sep=""))
      } # end (if Z/W)
    } # end (for i)

    if (lag == 2) {
      cat(rep("\n",2), "  # -- Estimate lagged effects between within-person centered variables (Lag = 2 waves) -- #")
      for (i in 3:no.waves) {
        if (Z == "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-2, "*w", X,i-2, " + YX", i,i-2, "*w", Y,i-2, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-2, "*w", X,i-2, " + YY", i,i-2, "*w", Y,i-2, sep=""))
        } else if (W != "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-2, "*w", X,i-2, " + YX", i,i-2, "*w", Y,i-2, " + ZX", i,i-2, "*w", Z,i-2, " + WX", i,i-2, "*w", W,i-2, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-2, "*w", X,i-2, " + YY", i,i-2, "*w", Y,i-2, " + ZY", i,i-2, "*w", Z,i-2, " + WY", i,i-2, "*w", W,i-2, sep=""))
          cat("\n", paste("  w", Z,i, " ~ XZ", i,i-2, "*w", X,i-2, " + YZ", i,i-2, "*w", Y,i-2, " + ZZ", i,i-2, "*w", Z,i-2, " + WZ", i,i-2, "*w", W,i-2, sep=""))
          cat("\n", paste("  w", W,i, " ~ XW", i,i-2, "*w", X,i-2, " + YW", i,i-2, "*w", Y,i-2, " + ZW", i,i-2, "*w", Z,i-2, " + WW", i,i-2, "*w", W,i-2, sep=""))
        } else if (Z != "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-2, "*w", X,i-2, " + YX", i,i-2, "*w", Y,i-2, " + ZX", i,i-2, "*w", Z,i-2, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-2, "*w", X,i-2, " + YY", i,i-2, "*w", Y,i-2, " + ZY", i,i-2, "*w", Z,i-2, sep=""))
          cat("\n", paste("  w", Z,i, " ~ XZ", i,i-2, "*w", X,i-2, " + YZ", i,i-2, "*w", Y,i-2, " + ZZ", i,i-2, "*w", Z,i-2, sep=""))
        } # end (if Z/W)
      } # end (for i)
    } # end (if lag == 2)

    if (lag == 3) {
      cat(rep("\n",2), "  # -- Estimate lagged effects between within-person centered variables (Lag = 3 waves) -- #")
      for (i in 4:no.waves) {
        if (Z == "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-3, "*w", X,i-3, " + YX", i,i-3, "*w", Y,i-3, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-3, "*w", X,i-3, " + YY", i,i-3, "*w", Y,i-3, sep=""))
        } else if (W != "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-3, "*w", X,i-3, " + YX", i,i-3, "*w", Y,i-3, " + ZX", i,i-3, "*w", Z,i-3, " + WX", i,i-3, "*w", W,i-3, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-3, "*w", X,i-3, " + YY", i,i-3, "*w", Y,i-3, " + ZY", i,i-3, "*w", Z,i-3, " + WY", i,i-3, "*w", W,i-3, sep=""))
          cat("\n", paste("  w", Z,i, " ~ XZ", i,i-3, "*w", X,i-3, " + YZ", i,i-3, "*w", Y,i-3, " + ZZ", i,i-3, "*w", Z,i-3, " + WZ", i,i-3, "*w", W,i-3, sep=""))
          cat("\n", paste("  w", W,i, " ~ XW", i,i-3, "*w", X,i-3, " + YW", i,i-3, "*w", Y,i-3, " + ZW", i,i-3, "*w", Z,i-3, " + WW", i,i-3, "*w", W,i-3, sep=""))
        } else if (Z != "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-3, "*w", X,i-3, " + YX", i,i-3, "*w", Y,i-3, " + ZX", i,i-3, "*w", Z,i-3, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-3, "*w", X,i-3, " + YY", i,i-3, "*w", Y,i-3, " + ZY", i,i-3, "*w", Z,i-3, sep=""))
          cat("\n", paste("  w", Z,i, " ~ XZ", i,i-3, "*w", X,i-3, " + YZ", i,i-3, "*w", Y,i-3, " + ZZ", i,i-3, "*w", Z,i-3, sep=""))
        } # end (if Z/W)
      } # end (for i)
    } # end (lag == 3)

    if (lag == 4) {
      cat(rep("\n",2), "  # -- Estimate lagged effects between within-person centered variables (Lag = 4 waves) -- #")
      for (i in 5:no.waves) {
        if (Z == "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-4, "*w", X,i-4, " + YX", i,i-4, "*w", Y,i-4, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-4, "*w", X,i-4, " + YY", i,i-4, "*w", Y,i-4, sep=""))
        } else if (W != "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-4, "*w", X,i-4, " + YX", i,i-4, "*w", Y,i-4, " + ZX", i,i-4, "*w", Z,i-4, " + WX", i,i-4, "*w", W,i-4, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-4, "*w", X,i-4, " + YY", i,i-4, "*w", Y,i-4, " + ZY", i,i-4, "*w", Z,i-4, " + WY", i,i-4, "*w", W,i-4, sep=""))
          cat("\n", paste("  w", Z,i, " ~ XZ", i,i-4, "*w", X,i-4, " + YZ", i,i-4, "*w", Y,i-4, " + ZZ", i,i-4, "*w", Z,i-4, " + WZ", i,i-4, "*w", W,i-4, sep=""))
          cat("\n", paste("  w", W,i, " ~ XW", i,i-4, "*w", X,i-4, " + YW", i,i-4, "*w", Y,i-4, " + ZW", i,i-4, "*w", Z,i-4, " + WW", i,i-4, "*w", W,i-4, sep=""))
        } else if (Z != "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-4, "*w", X,i-4, " + YX", i,i-4, "*w", Y,i-4, " + ZX", i,i-4, "*w", Z,i-4, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-4, "*w", X,i-4, " + YY", i,i-4, "*w", Y,i-4, " + ZY", i,i-4, "*w", Z,i-4, sep=""))
          cat("\n", paste("  w", Z,i, " ~ XZ", i,i-4, "*w", X,i-4, " + YZ", i,i-4, "*w", Y,i-4, " + ZZ", i,i-4, "*w", Z,i-4, sep=""))
        } # end (if Z/W)
      } # end (for i)
    } # end (if lag == 4)

    cat(rep("\n",2), "  '")

    # -- Run Model CLPMMLR -- #
    cat(rep("\n",2), "# -- Run Model CLPMMLR -- #")
    cat(rep("\n",2), "  CLPMMLR.fit <- suppressWarnings(lavaan::sem(CLPM,")
    cat("\n", "   ", data.source, ",")
    cat("\n", "   missing = 'fiml',")
    cat("\n", "   meanstructure = TRUE,")
    cat("\n", "   information = 'observed',")
    cat("\n", "   estimator = 'MLR'")
    cat("\n", "))")

  sink()  # Stop writing to file CLPM.txt
  ## ------------------------------- ##


  ## -- Execute CLPM.txt and request summary outputs-- ##
  source('CLPM.txt')
  print(lavaan::summary(CLPMMLR.fit, fit.measure = TRUE, standardized = TRUE, rsq = TRUE))
  if (lavaan::lavInspect(CLPMMLR.fit, what ="post.check") == FALSE) stop("The lavaan solution is non-admissible.")
  ## ------------------------------- ##


  ## -- Monte Carlo Simulation -- ##

  parEst <- lavaan::parameterEstimates(CLPMMLR.fit, remove.nonfree = TRUE)
  pest2 <- parEst[,5]  # Estimated Parameters
  pest3 <- lavaan::lavTech(CLPMMLR.fit, what = "vcov", add.labels = TRUE)  # Estimated Variance-Covariance of Estimated Parameters

  varI.eq = TRUE
  Invariance(parEst, pest2, pest3, no.path, MIset, no.compare, no.waves, lag, varI.eq, p, X, Y, Z, W)

  cat(rep("\n", 2))


  ## ----- Start writing script to CLPM.txt ----- ##
  sink('CLPM.txt')
    cat("\n", "# == Specify the model (CLPM) == #", "\n")
    cat("\n", "CLPM <- '")
    cat(rep("\n",2), "  # -- Constrain residual variance of observed variables to zero -- #")
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

    cat(rep("\n",2), "  # -- Estimate lagged effects between latent variables (Lag = 1 wave) -- #")
    cat("\n", "  #############################################")
    cat("\n", "  # Remove the subscripts for invariant paths #")
    cat("\n", "  #############################################")
    for (i in 2:no.waves) {
      if (Z == "NULL") {
        cat("\n", paste("  w", X,i, " ~ XX", i,i-1, "*w", X,i-1, " + YX", i,i-1, "*w", Y,i-1, sep=""))
        cat("\n", paste("  w", Y,i, " ~ XY", i,i-1, "*w", X,i-1, " + YY", i,i-1, "*w", Y,i-1, sep=""))
      } else if (W != "NULL") {
        cat("\n", paste("  w", X,i, " ~ XX", i,i-1, "*w", X,i-1, " + YX", i,i-1, "*w", Y,i-1, " + ZX", i,i-1, "*w", Z,i-1, " + WX", i,i-1, "*w", W,i-1, sep=""))
        cat("\n", paste("  w", Y,i, " ~ XY", i,i-1, "*w", X,i-1, " + YY", i,i-1, "*w", Y,i-1, " + ZY", i,i-1, "*w", Z,i-1, " + WY", i,i-1, "*w", W,i-1, sep=""))
        cat("\n", paste("  w", Z,i, " ~ XZ", i,i-1, "*w", X,i-1, " + YZ", i,i-1, "*w", Y,i-1, " + ZZ", i,i-1, "*w", Z,i-1, " + WZ", i,i-1, "*w", W,i-1, sep=""))
        cat("\n", paste("  w", W,i, " ~ XW", i,i-1, "*w", X,i-1, " + YW", i,i-1, "*w", Y,i-1, " + ZW", i,i-1, "*w", Z,i-1, " + WW", i,i-1, "*w", W,i-1, sep=""))
      } else if (Z != "NULL") {
        cat("\n", paste("  w", X,i, " ~ XX", i,i-1, "*w", X,i-1, " + YX", i,i-1, "*w", Y,i-1, " + ZX", i,i-1, "*w", Z,i-1, sep=""))
        cat("\n", paste("  w", Y,i, " ~ XY", i,i-1, "*w", X,i-1, " + YY", i,i-1, "*w", Y,i-1, " + ZY", i,i-1, "*w", Z,i-1, sep=""))
        cat("\n", paste("  w", Z,i, " ~ XZ", i,i-1, "*w", X,i-1, " + YZ", i,i-1, "*w", Y,i-1, " + ZZ", i,i-1, "*w", Z,i-1, sep=""))
      } # end (if Z)
    } # end (for i)

    if (lag == 2) {
      cat(rep("\n",2), "  # -- Estimate lagged effects between latent variables (Lag = 2 waves) -- #")
      cat("\n", "  #############################################")
      cat("\n", "  # Remove the subscripts for invariant paths #")
      cat("\n", "  #############################################")
      for (i in 3:no.waves) {
        if (Z == "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-2, "*w", X,i-2, " + YX", i,i-2, "*w", Y,i-2, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-2, "*w", X,i-2, " + YY", i,i-2, "*w", Y,i-2, sep=""))
        } else if (W != "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-2, "*w", X,i-2, " + YX", i,i-2, "*w", Y,i-2, " + ZX", i,i-2, "*w", Z,i-2, " + WX", i,i-2, "*w", W,i-2, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-2, "*w", X,i-2, " + YY", i,i-2, "*w", Y,i-2, " + ZY", i,i-2, "*w", Z,i-2, " + WY", i,i-2, "*w", W,i-2, sep=""))
          cat("\n", paste("  w", Z,i, " ~ XZ", i,i-2, "*w", X,i-2, " + YZ", i,i-2, "*w", Y,i-2, " + ZZ", i,i-2, "*w", Z,i-2, " + WZ", i,i-2, "*w", W,i-2, sep=""))
          cat("\n", paste("  w", W,i, " ~ XW", i,i-2, "*w", X,i-2, " + YW", i,i-2, "*w", Y,i-2, " + ZW", i,i-2, "*w", Z,i-2, " + WW", i,i-2, "*w", W,i-2, sep=""))
        } else if (Z != "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-2, "*w", X,i-2, " + YX", i,i-2, "*w", Y,i-2, " + ZX", i,i-2, "*w", Z,i-2, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-2, "*w", X,i-2, " + YY", i,i-2, "*w", Y,i-2, " + ZY", i,i-2, "*w", Z,i-2, sep=""))
          cat("\n", paste("  w", Z,i, " ~ XZ", i,i-2, "*w", X,i-2, " + YZ", i,i-2, "*w", Y,i-2, " + ZZ", i,i-2, "*w", Z,i-2, sep=""))
        } # end (if Z/W)
      } # end (for i)
    } # end (if lag == 2)

    if (lag == 3) {
      cat(rep("\n",2), "  # -- Estimate lagged effects between latent variables (Lag = 3 waves) -- #")
      cat("\n", "  #############################################")
      cat("\n", "  # Remove the subscripts for invariant paths #")
      cat("\n", "  #############################################")
      for (i in 3:no.waves) {
        if (Z == "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-3, "*w", X,i-3, " + YX", i,i-3, "*w", Y,i-3, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-3, "*w", X,i-3, " + YY", i,i-3, "*w", Y,i-3, sep=""))
        } else if (W != "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-3, "*w", X,i-3, " + YX", i,i-3, "*w", Y,i-3, " + ZX", i,i-3, "*w", Z,i-3, " + WX", i,i-3, "*w", W,i-3, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-3, "*w", X,i-3, " + YY", i,i-3, "*w", Y,i-3, " + ZY", i,i-3, "*w", Z,i-3, " + WY", i,i-3, "*w", W,i-3, sep=""))
          cat("\n", paste("  w", Z,i, " ~ XZ", i,i-3, "*w", X,i-3, " + YZ", i,i-3, "*w", Y,i-3, " + ZZ", i,i-3, "*w", Z,i-3, " + WZ", i,i-3, "*w", W,i-3, sep=""))
          cat("\n", paste("  w", W,i, " ~ XW", i,i-3, "*w", X,i-3, " + YW", i,i-3, "*w", Y,i-3, " + ZW", i,i-3, "*w", Z,i-3, " + WW", i,i-3, "*w", W,i-3, sep=""))
        } else if (Z != "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-3, "*w", X,i-3, " + YX", i,i-3, "*w", Y,i-3, " + ZX", i,i-3, "*w", Z,i-3, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-3, "*w", X,i-3, " + YY", i,i-3, "*w", Y,i-3, " + ZY", i,i-3, "*w", Z,i-3, sep=""))
          cat("\n", paste("  w", Z,i, " ~ XZ", i,i-3, "*w", X,i-3, " + YZ", i,i-3, "*w", Y,i-3, " + ZZ", i,i-3, "*w", Z,i-3, sep=""))
        } # end (if Z/W)
      } # end (for i)
    } # end (lag == 3)


    if (lag == 4) {
      cat(rep("\n",2), "  # -- Estimate lagged effects between latent variables (Lag = 4 waves) -- #")
      cat("\n", "  #############################################")
      cat("\n", "  # Remove the subscripts for invariant paths #")
      cat("\n", "  #############################################")
      for (i in 3:no.waves) {
        if (Z == "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-4, "*w", X,i-4, " + YX", i,i-4, "*w", Y,i-4, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-4, "*w", X,i-4, " + YY", i,i-4, "*w", Y,i-4, sep=""))
        } else if (W != "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-4, "*w", X,i-4, " + YX", i,i-4, "*w", Y,i-4, " + ZX", i,i-4, "*w", Z,i-4, " + WX", i,i-4, "*w", W,i-4, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-4, "*w", X,i-4, " + YY", i,i-4, "*w", Y,i-4, " + ZY", i,i-4, "*w", Z,i-4, " + WY", i,i-4, "*w", W,i-4, sep=""))
          cat("\n", paste("  w", Z,i, " ~ XZ", i,i-4, "*w", X,i-4, " + YZ", i,i-4, "*w", Y,i-4, " + ZZ", i,i-4, "*w", Z,i-4, " + WZ", i,i-4, "*w", W,i-4, sep=""))
          cat("\n", paste("  w", W,i, " ~ XW", i,i-4, "*w", X,i-4, " + YW", i,i-4, "*w", Y,i-4, " + ZW", i,i-4, "*w", Z,i-4, " + WW", i,i-4, "*w", W,i-4, sep=""))
        } else if (Z != "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-4, "*w", X,i-4, " + YX", i,i-4, "*w", Y,i-4, " + ZX", i,i-4, "*w", Z,i-4, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-4, "*w", X,i-4, " + YY", i,i-4, "*w", Y,i-4, " + ZY", i,i-4, "*w", Z,i-4, sep=""))
          cat("\n", paste("  w", Z,i, " ~ XZ", i,i-4, "*w", X,i-4, " + YZ", i,i-4, "*w", Y,i-4, " + ZZ", i,i-4, "*w", Z,i-4, sep=""))
        } # end (if Z/W)
      } # end (for i)
    } # end (lag == 4)


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


    cat(rep("\n",2), "  # -- Estimate (residual) variance of latent variables -- #")
    for (i in 1:no.waves) {
      cat("\n", paste("  w", X, i, " ~~ ", "w", X, i, sep=""))
      cat("\n", paste("  w", Y, i, " ~~ ", "w", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  w", Z, i, " ~~ ", "w", Z, i, sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  w", W, i, " ~~ ", "w", W, i, sep=""))
      } # end (if W)
    } # end (for i)

    cat(rep("\n",2), "  # -- Estimate grand means (intercepts) of observed variables -- #")
    for (i in 1:no.waves) {
      cat("\n", paste("  ", X, i, " ~ M", X, i, "*1", sep=""))
      cat("\n", paste("  ", Y, i, " ~ M", Y, i, "*1", sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  ", Z, i, " ~ M", Z, i, "*1", sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  ", W, i, " ~ M", W, i, "*1", sep=""))
      } # end (if W)
    } # end (for i)


    cat(rep("\n",2), "  ##########################################")
    cat("\n", "  # Regression of observed variables on C1 #")
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
    cat(rep("\n",2), "  print(lavaan::summary(CLPMMLR.fit, fit.measure = TRUE, standardized = TRUE, rsq = TRUE))", "\n")

    cat(rep("\n",2))

  sink() # Stop writing to file

}  # end (Function CLPM)

## ========================================================================================== ##





# ==================== Creating Function "RICLPM" ==================== #
#' Function Random Intercept Cross-Lagged Panel Model (RICLPM)
#'
#' Random Intercept Cross-Lagged Panel Model (RI-CLPM)
#'
#' @param data.source name of data.frame
#' @param no.waves number of waves (minimum = 3, must be grater than lag)
#' @param lag number of waves between two lags (minimum = 1, maximum = 4)
#' @param p critical p-value for pairwise comparisons (default is 0.001)
#' @param X name of variable X.
#' @param Y name of variable Y.
#' @param Z name of variable Z.
#' @param W name of variable W (Z must come before W).
#'
#' @return RICLPM outputs.
#' @export
#' @examples
#'
#' ## -- Example -- ##
#'
#' RICLPM(data.source="Data_A", 7, 2, X="EXPOSE", Y="INTENS")
#'

RICLPM <- function(data.source, no.waves, lag=1, p = 0.001, X, Y, Z="NULL", W = "NULL") {

  ## -- Check inputs -- ##

  if (no.waves < 3) stop("Minimum number of waves is 3")

  if (lag < 1) stop("Minimum number of lag is 1")
  if (lag > 4) stop("Maximum number of lag is 4")

  if ((no.waves - lag) < 1) stop("Number of waves must be greater than lag plus 1")

  if (p > 0.05) stop("p > 0.05 is not recommended")
  if (p < 0.0001) stop("p < 0.0001 is not recommended")

  if (Z == "NULL" & W != "NULL") stop("Z must be defined before W")


  ## ----- Creating Model RICLPM ----- ###
  sink('RICLPM.txt') # Start writing script to RICLPM.txt

    cat("\n", "## ----- Specify the model (RICLPM) ----- ##", "\n")
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


    # -- Constrain Residual Variance of Observed Variables to Zero -- #
    cat(rep("\n",2), "  # -- Constrain residual variance of observed variables to zero -- #")
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

    # -- Create Latent Variables from Observed Variables -- #
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

    # -- Estimate Covariances Between Residuals of Latent Variables -- #
    cat(rep("\n",2), "  # -- Estimate covariances between residuals of latent variables -- #")
    for (i in 2:no.waves) {
      cat("\n", paste("  w", X, i, " ~~ eXY", i, "*w", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  w", X, i, " ~~ eXZ", i, "*w", Z, i, sep=""))
        cat("\n", paste("  w", Y, i, " ~~ eYZ", i, "*w", Z, i, sep=""))
      }  # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  w", X, i, " ~~ eXW", i, "*w", W, i, sep=""))
        cat("\n", paste("  w", Y, i, " ~~ eYW", i, "*w", W, i, sep=""))
        cat("\n", paste("  w", Z, i, " ~~ eZW", i, "*w", W, i, sep=""))
      }  # end (if W)
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


    # -- Estimate Variance and Covariance of Random Intercepts -- #
    cat(rep("\n",2), "  # -- Estimate variance and covariance of random intercepts -- #")
    cat("\n", "   RI", X, " ~~ RI", X, sep="")
    cat("\n", "   RI", Y, " ~~ RI", Y, sep="")
    if (Z != "NULL") {
      cat("\n", "   RI", Z, " ~~ RI", Z, sep="")
    } # end (if Z != "NULL")
    if (W != "NULL") {
      cat("\n", "   RI", W, " ~~ RI", W, sep="")
    } # end (if W != "NULL")

    cat("\n", "   RI", X, " ~~ RI", Y, sep="")
    if (Z != "NULL") {
      cat("\n", "   RI", X, " ~~ RI", Z, sep="")
      cat("\n", "   RI", Y, " ~~ RI", Z, sep="")
    } # end (if Z != "NULL")
    if (W != "NULL") {
      cat("\n", "   RI", X, " ~~ RI", W, sep="")
      cat("\n", "   RI", Y, " ~~ RI", W, sep="")
      cat("\n", "   RI", Z, " ~~ RI", W, sep="")
    } # end (if W != "NULL")


    # -- Estimate (Residual) Variance of Latent Variables -- #
    cat(rep("\n",2), "  # -- Estimate (residual) variance of latent variables -- #")
    for (i in 1:no.waves) {
      cat("\n", paste("  w", X, i, " ~~ ", "w", X, i, sep=""))
      cat("\n", paste("  w", Y, i, " ~~ ", "w", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  w", Z, i, " ~~ ", "w", Z, i, sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  w", W, i, " ~~ ", "w", W, i, sep=""))
      } # end (if W)
    } # end (for i)

    # -- Estimate Grand Means (Intercepts) of Observed Variables -- #
    cat(rep("\n",2), "  # -- Estimate grand means (intercepts) of observed variables -- #")
    for (i in 1:no.waves) {
      cat("\n", paste("  ", X, i, " ~ M", X, i, "*1", sep=""))
      cat("\n", paste("  ", Y, i, " ~ M", Y, i, "*1", sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  ", Z, i, " ~ M", Z, i, "*1", sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  ", W, i, " ~ M", W, i, "*1", sep=""))
      } # (if W)
    } # end (for i)

    # -- Constrain covariance among RIX RIY, RIZ, RIW and wx1, wy1, wz1, ww1 -- #
    cat(rep("\n",2), "  # -- Constrain covariance among RIx, RIy, RIz, RIw and wx1, wy1, wz1, ww1 -- #")
    cat("\n", "   RI", X, " ~~ 0*w", X, "1", sep="")
    cat("\n", "   RI", X, " ~~ 0*w", Y, "1", sep="")
    if (Z != "NULL") {
      cat("\n", "   RI", X, " ~~ 0*w", Z, "1", sep="")
    } # end (if Z != "NULL")
    if (W != "NULL") {
      cat("\n", "   RI", X, " ~~ 0*w", W, "1", sep="")
    } # end (if W != "NULL")

    cat("\n", "   RI", Y, " ~~ 0*w", X, "1", sep="")
    cat("\n", "   RI", Y, " ~~ 0*w", Y, "1", sep="")
    if (Z != "NULL") {
      cat("\n", "   RI", Y, " ~~ 0*w", Z, "1", sep="")
    } # end (if Z != "NULL")
    if (W != "NULL") {
      cat("\n", "   RI", Y, " ~~ 0*w", W, "1", sep="")
    } # end (if W != "NULL")
    if (Z != "NULL") {
      cat("\n", "   RI", Z, " ~~ 0*w", X, "1", sep="")
      cat("\n", "   RI", Z, " ~~ 0*w", Y, "1", sep="")
      cat("\n", "   RI", Z, " ~~ 0*w", Z, "1", sep="")
    } # end (if Z != "NULL")
    if (W != "NULL") {
      cat("\n", "   RI", Z, " ~~ 0*w", W, "1", sep="")
      cat("\n", "   RI", W, " ~~ 0*w", X, "1", sep="")
      cat("\n", "   RI", W, " ~~ 0*w", Y, "1", sep="")
      cat("\n", "   RI", W, " ~~ 0*w", Z, "1", sep="")
      cat("\n", "   RI", W, " ~~ 0*w", W, "1", sep="")
    } # end (if W != "NULL")

    # -- Estimate Lagged Effects Between Latent Variables -- #
    cat(rep("\n",2), "  # -- Estimate lagged effects between latent variables (Lag = 1 wave) -- #")
    for (i in 2:no.waves) {
      if (Z == "NULL") {
        cat("\n", paste("  w", X, i, " ~ XX", i,i-1, "*w", X,i-1, " + YX", i,i-1, "*w", Y,i-1, sep=""))
        cat("\n", paste("  w", Y, i, " ~ XY", i,i-1, "*w", X,i-1, " + YY", i,i-1, "*w", Y,i-1, sep=""))
      } else if (W != "NULL") {
        cat("\n", paste("  w", X, i, " ~ XX", i,i-1, "*w", X,i-1, " + YX", i,i-1, "*w", Y,i-1, " + ZX", i,i-1, "*w", Z,i-1, " + WX", i,i-1, "*w", W,i-1, sep=""))
        cat("\n", paste("  w", Y, i, " ~ XY", i,i-1, "*w", X,i-1, " + YY", i,i-1, "*w", Y,i-1, " + ZY", i,i-1, "*w", Z,i-1, " + WY", i,i-1, "*w", W,i-1, sep=""))
        cat("\n", paste("  w", Z, i, " ~ XZ", i,i-1, "*w", X,i-1, " + YZ", i,i-1, "*w", Y,i-1, " + ZZ", i,i-1, "*w", Z,i-1, " + WZ", i,i-1, "*w", W,i-1, sep=""))
        cat("\n", paste("  w", W, i, " ~ XW", i,i-1, "*w", X,i-1, " + YW", i,i-1, "*w", Y,i-1, " + ZW", i,i-1, "*w", Z,i-1, " + WW", i,i-1, "*w", W,i-1, sep=""))
      } else if (Z != "NULL") {
        cat("\n", paste("  w", X, i, " ~ XX", i,i-1, "*w", X,i-1, " + YX", i,i-1, "*w", Y,i-1, " + ZX", i,i-1, "*w", Z,i-1, sep=""))
        cat("\n", paste("  w", Y, i, " ~ XY", i,i-1, "*w", X,i-1, " + YY", i,i-1, "*w", Y,i-1, " + ZY", i,i-1, "*w", Z,i-1, sep=""))
        cat("\n", paste("  w", Z, i, " ~ XZ", i,i-1, "*w", X,i-1, " + YZ", i,i-1, "*w", Y,i-1, " + ZZ", i,i-1, "*w", Z,i-1, sep=""))
      } # end (if Z/W)
    } # end (for i)

    if (lag == 2) {
      cat(rep("\n",2), "  # -- Estimate lagged effects between within-person centered variables (Lag = 2 waves) -- #")
      for (i in 3:no.waves) {
        if (Z == "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-2, "*w", X,i-2, " + YX", i,i-2, "*w", Y,i-2, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-2, "*w", X,i-2, " + YY", i,i-2, "*w", Y,i-2, sep=""))
        } else if (W != "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-2, "*w", X,i-2, " + YX", i,i-2, "*w", Y,i-2, " + ZX", i,i-2, "*w", Z,i-2, " + WX", i,i-2, "*w", W,i-2, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-2, "*w", X,i-2, " + YY", i,i-2, "*w", Y,i-2, " + ZY", i,i-2, "*w", Z,i-2, " + WY", i,i-2, "*w", W,i-2, sep=""))
          cat("\n", paste("  w", Z,i, " ~ XZ", i,i-2, "*w", X,i-2, " + YZ", i,i-2, "*w", Y,i-2, " + ZZ", i,i-2, "*w", Z,i-2, " + WZ", i,i-2, "*w", W,i-2, sep=""))
          cat("\n", paste("  w", W,i, " ~ XW", i,i-2, "*w", X,i-2, " + YW", i,i-2, "*w", Y,i-2, " + ZW", i,i-2, "*w", Z,i-2, " + WW", i,i-2, "*w", W,i-2, sep=""))
        } else if (Z != "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-2, "*w", X,i-2, " + YX", i,i-2, "*w", Y,i-2, " + ZX", i,i-2, "*w", Z,i-2, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-2, "*w", X,i-2, " + YY", i,i-2, "*w", Y,i-2, " + ZY", i,i-2, "*w", Z,i-2, sep=""))
          cat("\n", paste("  w", Z,i, " ~ XZ", i,i-2, "*w", X,i-2, " + YZ", i,i-2, "*w", Y,i-2, " + ZZ", i,i-2, "*w", Z,i-2, sep=""))
        } # end (if Z/W)
      } # end (for i)
    } # end (if lag == 2)

    if (lag == 3) {
      cat(rep("\n",2), "  # -- Estimate lagged effects between within-person centered variables (Lag = 3 waves) -- #")
      for (i in 4:no.waves) {
        if (Z == "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-3, "*w", X,i-3, " + YX", i,i-3, "*w", Y,i-3, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-3, "*w", X,i-3, " + YY", i,i-3, "*w", Y,i-3, sep=""))
        } else if (W != "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-3, "*w", X,i-3, " + YX", i,i-3, "*w", Y,i-3, " + ZX", i,i-3, "*w", Z,i-3, " + WX", i,i-3, "*w", W,i-3, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-3, "*w", X,i-3, " + YY", i,i-3, "*w", Y,i-3, " + ZY", i,i-3, "*w", Z,i-3, " + WY", i,i-3, "*w", W,i-3, sep=""))
          cat("\n", paste("  w", Z,i, " ~ XZ", i,i-3, "*w", X,i-3, " + YZ", i,i-3, "*w", Y,i-3, " + ZZ", i,i-3, "*w", Z,i-3, " + WZ", i,i-3, "*w", W,i-3, sep=""))
          cat("\n", paste("  w", W,i, " ~ XW", i,i-3, "*w", X,i-3, " + YW", i,i-3, "*w", Y,i-3, " + ZW", i,i-3, "*w", Z,i-3, " + WW", i,i-3, "*w", W,i-3, sep=""))
        } else if (Z != "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-3, "*w", X,i-3, " + YX", i,i-3, "*w", Y,i-3, " + ZX", i,i-3, "*w", Z,i-3, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-3, "*w", X,i-3, " + YY", i,i-3, "*w", Y,i-3, " + ZY", i,i-3, "*w", Z,i-3, sep=""))
          cat("\n", paste("  w", Z,i, " ~ XZ", i,i-3, "*w", X,i-3, " + YZ", i,i-3, "*w", Y,i-3, " + ZZ", i,i-3, "*w", Z,i-3, sep=""))
        } # end (if Z/W)
      } # end (for i)
    } # end (lag == 3)

    if (lag == 4) {
      cat(rep("\n",2), "  # -- Estimate lagged effects between within-person centered variables (Lag = 4 waves) -- #")
      for (i in 5:no.waves) {
        if (Z == "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-4, "*w", X,i-4, " + YX", i,i-4, "*w", Y,i-4, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-4, "*w", X,i-4, " + YY", i,i-4, "*w", Y,i-4, sep=""))
        } else if (W != "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-4, "*w", X,i-4, " + YX", i,i-4, "*w", Y,i-4, " + ZX", i,i-4, "*w", Z,i-4, " + WX", i,i-4, "*w", W,i-4, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-4, "*w", X,i-4, " + YY", i,i-4, "*w", Y,i-4, " + ZY", i,i-4, "*w", Z,i-4, " + WY", i,i-4, "*w", W,i-4, sep=""))
          cat("\n", paste("  w", Z,i, " ~ XZ", i,i-4, "*w", X,i-4, " + YZ", i,i-4, "*w", Y,i-4, " + ZZ", i,i-4, "*w", Z,i-4, " + WZ", i,i-4, "*w", W,i-4, sep=""))
          cat("\n", paste("  w", W,i, " ~ XW", i,i-4, "*w", X,i-4, " + YW", i,i-4, "*w", Y,i-4, " + ZW", i,i-4, "*w", Z,i-4, " + WW", i,i-4, "*w", W,i-4, sep=""))
        } else if (Z != "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-4, "*w", X,i-4, " + YX", i,i-4, "*w", Y,i-4, " + ZX", i,i-4, "*w", Z,i-4, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-4, "*w", X,i-4, " + YY", i,i-4, "*w", Y,i-4, " + ZY", i,i-4, "*w", Z,i-4, sep=""))
          cat("\n", paste("  w", Z,i, " ~ XZ", i,i-4, "*w", X,i-4, " + YZ", i,i-4, "*w", Y,i-4, " + ZZ", i,i-4, "*w", Z,i-4, sep=""))
        } # end (if Z/W)
      } # end (for i)
    } # end (if lag == 4)

    cat(rep("\n",2), "  '")

    # -- Run Model RICLPMMLR -- #
    cat(rep("\n",2), "# -- Run Model RICLPMMLR -- #")
    cat(rep("\n",2), "  RICLPMMLR.fit <- suppressWarnings(lavaan::sem(RICLPM,")
    cat("\n", "   ", data.source, ",")
    cat("\n", "   missing = 'fiml',")
    cat("\n", "   meanstructure = TRUE,")
    cat("\n", "   information = 'observed',")
    cat("\n", "   estimator = 'MLR'")
    cat("\n", "))")

  sink()  # Stop writing to file RICLPM.txt
  ## ------------------------------- ##


  ## -- Execute RICLPM.txt and request summary outputs-- ##
  source('RICLPM.txt')
  print(lavaan::summary(RICLPMMLR.fit, fit.measure = TRUE, standardized = TRUE, rsq = TRUE))
  if (lavaan::lavInspect(RICLPMMLR.fit, what ="post.check") == FALSE) stop("The lavaan solution is non-admissible.")
  ## ------------------------------- ##


  ## -- Monte Carlo Simulation -- ##

  parEst <- lavaan::parameterEstimates(RICLPMMLR.fit, remove.nonfree = TRUE)
  pest2 <- parEst[,5]  # Estimated Parameters
  pest3 <- lavaan::lavTech(RICLPMMLR.fit, what = "vcov", add.labels = TRUE)  # Estimated Variance-Covariance of Estimated Parameters

  varI.eq = TRUE
  Invariance(parEst, pest2, pest3, no.path, MIset, no.compare, no.waves, lag, varI.eq, p, X, Y, Z, W)

  cat(rep("\n", 2))


  ## ----- Start writing script to RICLPM.txt ----- ##
  sink('RICLPM.txt')
    cat("\n", "# Specify the model (RICLPM)", "\n")
    cat("\n", "RICLPM <- '")

    cat(rep("\n",2), "  # -- Create between components (random intercepts) -- #")
    BX <- paste("  RI", X, " =~ 1*", X, "1", sep="")
    BY <- paste("  RI", Y, " =~ 1*", Y, "1", sep="")
    if (Z != "NULL") {
      BY <- paste("  RI", Z, " =~ 1*", Z, "1", sep="")
    } # end (if Z != "NULL")
    if (W != "NULL") {
      BY <- paste("  RI", W, " =~ 1*", W, "1", sep="")
    } # end (if W != "NULL")

    for (i in 2:no.waves) {
      BX <- paste(BX, " +1*", X, i, sep="")
      BY <- paste(BY, " +1*", Y, i, sep="")
      if (Z != "NULL") {
        BZ <- paste(BZ, " +1*", Z, i, sep="")
      } # end (if Z != "NULL")
      if (W != "NULL") {
        BW <- paste(BW, " +1*", W, i, sep="")
      } # end (if W != "NULL")
    } # end (for i)

    cat("\n", BX)
    cat("\n", BY)
    if (Z != "NULL") {
      cat("\n", BZ)
    } # end (if Z != "NULL")
    if (W != "NULL") {
      cat("\n", BW)
    } # end (if W != "NULL")

    # -- Constrain Residual Variance of Observed Variables to Zero -- #
    cat(rep("\n",2), "  # -- Constrain residual variance of observed variables to zero -- #")
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

    # -- Create Latent Variables from Observed Variables -- #
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

    cat(rep("\n",2), "  # -- Estimate lagged effects between latent variables (Lag = 1 wave) -- #")
    cat("\n", "  #############################################")
    cat("\n", "  # Remove the subscripts for invariant paths #")
    cat("\n", "  #############################################")
    for (i in 2:no.waves) {
      if (Z == "NULL") {
        cat("\n", paste("  w", X,i, " ~ XX", i,i-1, "*w", X,i-1, " + YX", i,i-1, "*w", Y,i-1, sep=""))
        cat("\n", paste("  w", Y,i, " ~ XY", i,i-1, "*w", X,i-1, " + YY", i,i-1, "*w", Y,i-1, sep=""))
      } else if (W != "NULL") {
        cat("\n", paste("  w", X,i, " ~ XX", i,i-1, "*w", X,i-1, " + YX", i,i-1, "*w", Y,i-1, " + ZX", i,i-1, "*w", Z,i-1, " + WX", i,i-1, "*w", W,i-1, sep=""))
        cat("\n", paste("  w", Y,i, " ~ XY", i,i-1, "*w", X,i-1, " + YY", i,i-1, "*w", Y,i-1, " + ZY", i,i-1, "*w", Z,i-1, " + WY", i,i-1, "*w", W,i-1, sep=""))
        cat("\n", paste("  w", Z,i, " ~ XZ", i,i-1, "*w", X,i-1, " + YZ", i,i-1, "*w", Y,i-1, " + ZZ", i,i-1, "*w", Z,i-1, " + WZ", i,i-1, "*w", W,i-1, sep=""))
        cat("\n", paste("  w", W,i, " ~ XW", i,i-1, "*w", X,i-1, " + YW", i,i-1, "*w", Y,i-1, " + ZW", i,i-1, "*w", Z,i-1, " + WW", i,i-1, "*w", W,i-1, sep=""))
      } else if (Z != "NULL") {
        cat("\n", paste("  w", X,i, " ~ XX", i,i-1, "*w", X,i-1, " + YX", i,i-1, "*w", Y,i-1, " + ZX", i,i-1, "*w", Z,i-1, sep=""))
        cat("\n", paste("  w", Y,i, " ~ XY", i,i-1, "*w", X,i-1, " + YY", i,i-1, "*w", Y,i-1, " + ZY", i,i-1, "*w", Z,i-1, sep=""))
        cat("\n", paste("  w", Z,i, " ~ XZ", i,i-1, "*w", X,i-1, " + YZ", i,i-1, "*w", Y,i-1, " + ZZ", i,i-1, "*w", Z,i-1, sep=""))
      } # end (if Z)
    } # end (for i)

    if (lag == 2) {
      cat(rep("\n",2), "  # -- Estimate lagged effects between latent variables (Lag = 2 waves) -- #")
      cat("\n", "  #############################################")
      cat("\n", "  # Remove the subscripts for invariant paths #")
      cat("\n", "  #############################################")
      for (i in 3:no.waves) {
        if (Z == "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-2, "*w", X,i-2, " + YX", i,i-2, "*w", Y,i-2, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-2, "*w", X,i-2, " + YY", i,i-2, "*w", Y,i-2, sep=""))
        } else if (W != "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-2, "*w", X,i-2, " + YX", i,i-2, "*w", Y,i-2, " + ZX", i,i-2, "*w", Z,i-2, " + WX", i,i-2, "*w", W,i-2, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-2, "*w", X,i-2, " + YY", i,i-2, "*w", Y,i-2, " + ZY", i,i-2, "*w", Z,i-2, " + WY", i,i-2, "*w", W,i-2, sep=""))
          cat("\n", paste("  w", Z,i, " ~ XZ", i,i-2, "*w", X,i-2, " + YZ", i,i-2, "*w", Y,i-2, " + ZZ", i,i-2, "*w", Z,i-2, " + WZ", i,i-2, "*w", W,i-2, sep=""))
          cat("\n", paste("  w", W,i, " ~ XW", i,i-2, "*w", X,i-2, " + YW", i,i-2, "*w", Y,i-2, " + ZW", i,i-2, "*w", Z,i-2, " + WW", i,i-2, "*w", W,i-2, sep=""))
        } else if (Z != "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-2, "*w", X,i-2, " + YX", i,i-2, "*w", Y,i-2, " + ZX", i,i-2, "*w", Z,i-2, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-2, "*w", X,i-2, " + YY", i,i-2, "*w", Y,i-2, " + ZY", i,i-2, "*w", Z,i-2, sep=""))
          cat("\n", paste("  w", Z,i, " ~ XZ", i,i-2, "*w", X,i-2, " + YZ", i,i-2, "*w", Y,i-2, " + ZZ", i,i-2, "*w", Z,i-2, sep=""))
        } # end (if Z/W)
      } # end (for i)
    } # end (if lag == 2)

    if (lag == 3) {
      cat(rep("\n",2), "  # -- Estimate lagged effects between latent variables (Lag = 3 waves) -- #")
      cat("\n", "  #############################################")
      cat("\n", "  # Remove the subscripts for invariant paths #")
      cat("\n", "  #############################################")
      for (i in 3:no.waves) {
        if (Z == "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-3, "*w", X,i-3, " + YX", i,i-3, "*w", Y,i-3, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-3, "*w", X,i-3, " + YY", i,i-3, "*w", Y,i-3, sep=""))
        } else if (W != "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-3, "*w", X,i-3, " + YX", i,i-3, "*w", Y,i-3, " + ZX", i,i-3, "*w", Z,i-3, " + WX", i,i-3, "*w", W,i-3, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-3, "*w", X,i-3, " + YY", i,i-3, "*w", Y,i-3, " + ZY", i,i-3, "*w", Z,i-3, " + WY", i,i-3, "*w", W,i-3, sep=""))
          cat("\n", paste("  w", Z,i, " ~ XZ", i,i-3, "*w", X,i-3, " + YZ", i,i-3, "*w", Y,i-3, " + ZZ", i,i-3, "*w", Z,i-3, " + WZ", i,i-3, "*w", W,i-3, sep=""))
          cat("\n", paste("  w", W,i, " ~ XW", i,i-3, "*w", X,i-3, " + YW", i,i-3, "*w", Y,i-3, " + ZW", i,i-3, "*w", Z,i-3, " + WW", i,i-3, "*w", W,i-3, sep=""))
        } else if (Z != "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-3, "*w", X,i-3, " + YX", i,i-3, "*w", Y,i-3, " + ZX", i,i-3, "*w", Z,i-3, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-3, "*w", X,i-3, " + YY", i,i-3, "*w", Y,i-3, " + ZY", i,i-3, "*w", Z,i-3, sep=""))
          cat("\n", paste("  w", Z,i, " ~ XZ", i,i-3, "*w", X,i-3, " + YZ", i,i-3, "*w", Y,i-3, " + ZZ", i,i-3, "*w", Z,i-3, sep=""))
        } # end (if Z/W)
      } # end (for i)
    } # end (lag == 3)

   if (lag == 4) {
      cat(rep("\n",2), "  # -- Estimate lagged effects between latent variables (Lag = 4 waves) -- #")
      cat("\n", "  #############################################")
      cat("\n", "  # Remove the subscripts for invariant paths #")
      cat("\n", "  #############################################")
      for (i in 3:no.waves) {
        if (Z == "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-4, "*w", X,i-4, " + YX", i,i-4, "*w", Y,i-4, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-4, "*w", X,i-4, " + YY", i,i-4, "*w", Y,i-4, sep=""))
        } else if (W != "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-4, "*w", X,i-4, " + YX", i,i-4, "*w", Y,i-4, " + ZX", i,i-4, "*w", Z,i-4, " + WX", i,i-4, "*w", W,i-4, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-4, "*w", X,i-4, " + YY", i,i-4, "*w", Y,i-4, " + ZY", i,i-4, "*w", Z,i-4, " + WY", i,i-4, "*w", W,i-4, sep=""))
          cat("\n", paste("  w", Z,i, " ~ XZ", i,i-4, "*w", X,i-4, " + YZ", i,i-4, "*w", Y,i-4, " + ZZ", i,i-4, "*w", Z,i-4, " + WZ", i,i-4, "*w", W,i-4, sep=""))
          cat("\n", paste("  w", W,i, " ~ XW", i,i-4, "*w", X,i-4, " + YW", i,i-4, "*w", Y,i-4, " + ZW", i,i-4, "*w", Z,i-4, " + WW", i,i-4, "*w", W,i-4, sep=""))
        } else if (Z != "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-4, "*w", X,i-4, " + YX", i,i-4, "*w", Y,i-4, " + ZX", i,i-4, "*w", Z,i-4, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-4, "*w", X,i-4, " + YY", i,i-4, "*w", Y,i-4, " + ZY", i,i-4, "*w", Z,i-4, sep=""))
          cat("\n", paste("  w", Z,i, " ~ XZ", i,i-4, "*w", X,i-4, " + YZ", i,i-4, "*w", Y,i-4, " + ZZ", i,i-4, "*w", Z,i-4, sep=""))
        } # end (if Z/W)
      } # end (for i)
    } # end (lag == 4)


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


    cat(rep("\n",2), "  # -- Estimate variance and covariance of random intercepts -- #")
    cat("\n", "   RI", X, " ~~ RI", X, sep="")
    cat("\n", "   RI", Y, " ~~ RI", Y, sep="")
    if (Z != "NULL") {
      cat("\n", "   RI", Z, " ~~ RI", Z, sep="")
    } # end (if Z != "NULL")
    if (W != "NULL") {
      cat("\n", "   RI", W, " ~~ RI", W, sep="")
    } # end (if W != "NULL")
    cat("\n", "   RI", X, " ~~ RI", Y, sep="")
    if (Z != "NULL") {
      cat("\n", "   RI", X, " ~~ RI", Z, sep="")
      cat("\n", "   RI", Y, " ~~ RI", Z, sep="")
    } # end (if Z != "NULL")
    if (W != "NULL") {
      cat("\n", "   RI", X, " ~~ RI", W, sep="")
      cat("\n", "   RI", Y, " ~~ RI", W, sep="")
      cat("\n", "   RI", Z, " ~~ RI", W, sep="")
    } # end (if W != "NULL")


    cat(rep("\n",2), "  # -- Estimate (residual) variance of latent variables -- #")
    for (i in 1:no.waves) {
      cat("\n", paste("  w", X, i, " ~~ ", "w", X, i, sep=""))
      cat("\n", paste("  w", Y, i, " ~~ ", "w", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  w", Z, i, " ~~ ", "w", Z, i, sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  w", W, i, " ~~ ", "w", W, i, sep=""))
      } # end (if W)
    } # end (for i)


    cat(rep("\n",2), "  # -- Estimate grand means (intercepts) of observed variables -- #")
    for (i in 1:no.waves) {
      cat("\n", paste("  ", X, i, " ~ M", X, i, "*1", sep=""))
      cat("\n", paste("  ", Y, i, " ~ M", Y, i, "*1", sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  ", Z, i, " ~ M", Z, i, "*1", sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  ", W, i, " ~ M", W, i, "*1", sep=""))
      } # end (if W)
    } # end (for i)


    # -- Constrain covariance among RIX RIY, RIZ, RIW and wx1, wy1, wz1, ww1 -- #
    cat(rep("\n",2), "  # -- Constrain covariance among RIx, RIy, RIz, RIw and wx1, wy1, wz1, ww1 -- #")
    cat("\n", "   RI", X, " ~~ 0*w", X, "1", sep="")
    cat("\n", "   RI", X, " ~~ 0*w", Y, "1", sep="")
    if (Z != "NULL") {
      cat("\n", "   RI", X, " ~~ 0*w", Z, "1", sep="")
    } # end (if Z)
    if (W != "NULL") {
      cat("\n", "   RI", X, " ~~ 0*w", W, "1", sep="")
    } # end (if W)
    cat("\n", "   RI", Y, " ~~ 0*w", X, "1", sep="")
    cat("\n", "   RI", Y, " ~~ 0*w", Y, "1", sep="")
    if (Z != "NULL") {
      cat("\n", "   RI", Y, " ~~ 0*w", Z, "1", sep="")
    } # end (if Z)
    if (W != "NULL") {
      cat("\n", "   RI", Y, " ~~ 0*w", W, "1", sep="")
    } # end (if W)
    if (Z != "NULL") {
      cat("\n", "   RI", Z, " ~~ 0*w", X, "1", sep="")
      cat("\n", "   RI", Z, " ~~ 0*w", Y, "1", sep="")
      cat("\n", "   RI", Z, " ~~ 0*w", Z, "1", sep="")
    } # end (if Z)
    if (W != "NULL") {
      cat("\n", "   RI", Z, " ~~ 0*w", W, "1", sep="")
      cat("\n", "   RI", W, " ~~ 0*w", X, "1", sep="")
      cat("\n", "   RI", W, " ~~ 0*w", Y, "1", sep="")
      cat("\n", "   RI", W, " ~~ 0*w", Z, "1", sep="")
      cat("\n", "   RI", W, " ~~ 0*w", W, "1", sep="")
    } # end (if W)


    cat(rep("\n",2), "  ##########################################")
    cat("\n", "  # Regression of observed variables on C1 #")
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
    cat(rep("\n",2), "  print(lavaan::summary(RICLPMMLR.fit, fit.measure = TRUE, standardized = TRUE, rsq = TRUE))", "\n")

    cat(rep("\n",2))

  sink() # Stop writing to file

}  # end (Function RICLPM)

## ========================================================================================== ##





# ==================== Creating Function "LGCMSR" ==================== #
#' Function Latent Growth Curve Model with Structural Residuals (LGCMSR)
#'
#' Latent Growth Curve Model with Structural Residuals (LGCM-SR)
#'
#' @param data.source name of data.frame
#' @param no.waves number of waves (minimum = 3, must be grater than lag)
#' @param lag number of waves between two lags (minimum = 1, maximum = 4)
#' @param p critical p-value for pairwise comparisons (default is 0.001)
#' @param X name of variable X.
#' @param Y name of variable Y.
#' @param Z name of variable Z.
#' @param W name of variable W (Z must come before W).
#'
#' @return LGCM-SR outputs.
#' @export
#' @examples
#'
#' ## -- Example -- ##
#'
#' LGCMSR(data.source="Data_A", 7, 2, X="EXPOSE", Y="INTENS")
#'

LGCMSR <- function(data.source, no.waves, lag=1, p = 0.001, X, Y, Z="NULL", W = "NULL") {

  ## -- Check inputs -- ##

  if (no.waves < 3) stop("Minimum number of waves is 3")

  if (lag < 1) stop("Minimum number of lag is 1")
  if (lag > 4) stop("Maximum number of lag is 4")

  if ((no.waves - lag) < 1) stop("Number of waves must be greater than lag plus 1")

  if (p > 0.05) stop("p > 0.05 is not recommended")
  if (p < 0.0001) stop("p < 0.0001 is not recommended")

  if (Z == "NULL" & W != "NULL") stop("Z must be defined before W")



  ## ----- Creating Model LGCMSR ----- ###
  sink('LGCMSR.txt') # Start writing script to LGCMSR.txt

    cat("\n", "## ----- Specify the model (LGCMSR) ----- ##", "\n")
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


    # -- Constrain Residual Variance of Observed Variables to Zero -- #
    cat(rep("\n",2), "  # -- Constrain residual variance of observed variables to zero -- #")
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

    # -- Create Latent Variables from Observed Variables -- #
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

    # -- Estimate Covariances Between Residuals of Latent Variables -- #
    cat(rep("\n",2), "  # -- Estimate covariances between residuals of latent variables -- #")
    for (i in 2:no.waves) {
      cat("\n", paste("  w", X, i, " ~~ eXY", i, "*w", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  w", X, i, " ~~ eXZ", i, "*w", Z, i, sep=""))
        cat("\n", paste("  w", Y, i, " ~~ eYZ", i, "*w", Z, i, sep=""))
      }  # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  w", X, i, " ~~ eXW", i, "*w", W, i, sep=""))
        cat("\n", paste("  w", Y, i, " ~~ eYW", i, "*w", W, i, sep=""))
        cat("\n", paste("  w", Z, i, " ~~ eZW", i, "*w", W, i, sep=""))
      }  # end (if W)
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


    # -- Estimate (Residual) Variance of Latent Variables -- #
    cat(rep("\n",2), "  # -- Estimate (residual) variance of latent variables -- #")
    for (i in 1:no.waves) {
      cat("\n", paste("  w", X, i, " ~~ ", "w", X, i, sep=""))
      cat("\n", paste("  w", Y, i, " ~~ ", "w", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  w", Z, i, " ~~ ", "w", Z, i, sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  w", W, i, " ~~ ", "w", W, i, sep=""))
      } # end (if W)
    } # end (for i)

    # -- Estimate Grand Means (Intercepts) of Observed Variables -- #
    cat(rep("\n",2), "  # -- Estimate grand means (intercepts) of observed variables -- #")
    for (i in 1:no.waves) {
      cat("\n", paste("  ", X, i, " ~ M", X, i, "*1", sep=""))
      cat("\n", paste("  ", Y, i, " ~ M", Y, i, "*1", sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  ", Z, i, " ~ M", Z, i, "*1", sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  ", W, i, " ~ M", W, i, "*1", sep=""))
      } # (if W)
    } # end (for i)

    # -- Constrain covariance among RIX RIY, RIZ, RIW and wx1, wy1, wz1, ww1 -- #
    cat(rep("\n",2), "  # -- Constrain covariance among RIx, RIy, RIz, RIw and wx1, wy1, wz1, ww1 -- #")
    cat("\n", "   RI", X, " ~~ 0*w", X, "1", sep="")
    cat("\n", "   RI", X, " ~~ 0*w", Y, "1", sep="")
    if (Z != "NULL") {
      cat("\n", "   RI", X, " ~~ 0*w", Z, "1", sep="")
    } # end (if Z != "NULL")
    if (W != "NULL") {
      cat("\n", "   RI", X, " ~~ 0*w", W, "1", sep="")
    } # end (if W != "NULL")

    cat("\n", "   RI", Y, " ~~ 0*w", X, "1", sep="")
    cat("\n", "   RI", Y, " ~~ 0*w", Y, "1", sep="")
    if (Z != "NULL") {
      cat("\n", "   RI", Y, " ~~ 0*w", Z, "1", sep="")
    } # end (if Z != "NULL")
    if (W != "NULL") {
      cat("\n", "   RI", Y, " ~~ 0*w", W, "1", sep="")
    } # end (if W != "NULL")
    if (Z != "NULL") {
      cat("\n", "   RI", Z, " ~~ 0*w", X, "1", sep="")
      cat("\n", "   RI", Z, " ~~ 0*w", Y, "1", sep="")
      cat("\n", "   RI", Z, " ~~ 0*w", Z, "1", sep="")
    } # end (if Z != "NULL")
    if (W != "NULL") {
      cat("\n", "   RI", Z, " ~~ 0*w", W, "1", sep="")
      cat("\n", "   RI", W, " ~~ 0*w", X, "1", sep="")
      cat("\n", "   RI", W, " ~~ 0*w", Y, "1", sep="")
      cat("\n", "   RI", W, " ~~ 0*w", Z, "1", sep="")
      cat("\n", "   RI", W, " ~~ 0*w", W, "1", sep="")
    } # end (if W != "NULL")


    # -- Constrain covariance among RSX RSY, RSZ, RSW and wx1, wy1, wz1, ww1 -- #
    cat(rep("\n",2), "  # -- Constrain covariance among RSx, RSy, RSz, RSw and wx1, wy1, wz1, ww1 -- #")
    cat("\n", "   RS", X, " ~~ 0*w", X, "1", sep="")
    cat("\n", "   RS", X, " ~~ 0*w", Y, "1", sep="")
    if (Z != "NULL") {
      cat("\n", "   RS", X, " ~~ 0*w", Z, "1", sep="")
    } # end (if Z != "NULL")
    if (W != "NULL") {
      cat("\n", "   RS", X, " ~~ 0*w", W, "1", sep="")
    } # end (if W != "NULL")

    cat("\n", "   RS", Y, " ~~ 0*w", X, "1", sep="")
    cat("\n", "   RS", Y, " ~~ 0*w", Y, "1", sep="")
    if (Z != "NULL") {
      cat("\n", "   RS", Y, " ~~ 0*w", Z, "1", sep="")
    } # end (if Z != "NULL")
    if (W != "NULL") {
      cat("\n", "   RS", Y, " ~~ 0*w", W, "1", sep="")
    } # end (if W != "NULL")
    if (Z != "NULL") {
      cat("\n", "   RS", Z, " ~~ 0*w", X, "1", sep="")
      cat("\n", "   RS", Z, " ~~ 0*w", Y, "1", sep="")
      cat("\n", "   RS", Z, " ~~ 0*w", Z, "1", sep="")
    } # end (if Z != "NULL")
    if (W != "NULL") {
      cat("\n", "   RS", Z, " ~~ 0*w", W, "1", sep="")
      cat("\n", "   RS", W, " ~~ 0*w", X, "1", sep="")
      cat("\n", "   RS", W, " ~~ 0*w", Y, "1", sep="")
      cat("\n", "   RS", W, " ~~ 0*w", Z, "1", sep="")
      cat("\n", "   RS", W, " ~~ 0*w", W, "1", sep="")
    } # end (if W != "NULL")


    # -- Estimate Lagged Effects Between Latent Variables -- #
    cat(rep("\n",2), "  # -- Estimate lagged effects between latent variables (Lag = 1 wave) -- #")
    for (i in 2:no.waves) {
      if (Z == "NULL") {
        cat("\n", paste("  w", X, i, " ~ XX", i,i-1, "*w", X,i-1, " + YX", i,i-1, "*w", Y,i-1, sep=""))
        cat("\n", paste("  w", Y, i, " ~ XY", i,i-1, "*w", X,i-1, " + YY", i,i-1, "*w", Y,i-1, sep=""))
      } else if (W != "NULL") {
        cat("\n", paste("  w", X, i, " ~ XX", i,i-1, "*w", X,i-1, " + YX", i,i-1, "*w", Y,i-1, " + ZX", i,i-1, "*w", Z,i-1, " + WX", i,i-1, "*w", W,i-1, sep=""))
        cat("\n", paste("  w", Y, i, " ~ XY", i,i-1, "*w", X,i-1, " + YY", i,i-1, "*w", Y,i-1, " + ZY", i,i-1, "*w", Z,i-1, " + WY", i,i-1, "*w", W,i-1, sep=""))
        cat("\n", paste("  w", Z, i, " ~ XZ", i,i-1, "*w", X,i-1, " + YZ", i,i-1, "*w", Y,i-1, " + ZZ", i,i-1, "*w", Z,i-1, " + WZ", i,i-1, "*w", W,i-1, sep=""))
        cat("\n", paste("  w", W, i, " ~ XW", i,i-1, "*w", X,i-1, " + YW", i,i-1, "*w", Y,i-1, " + ZW", i,i-1, "*w", Z,i-1, " + WW", i,i-1, "*w", W,i-1, sep=""))
      } else if (Z != "NULL") {
        cat("\n", paste("  w", X, i, " ~ XX", i,i-1, "*w", X,i-1, " + YX", i,i-1, "*w", Y,i-1, " + ZX", i,i-1, "*w", Z,i-1, sep=""))
        cat("\n", paste("  w", Y, i, " ~ XY", i,i-1, "*w", X,i-1, " + YY", i,i-1, "*w", Y,i-1, " + ZY", i,i-1, "*w", Z,i-1, sep=""))
        cat("\n", paste("  w", Z, i, " ~ XZ", i,i-1, "*w", X,i-1, " + YZ", i,i-1, "*w", Y,i-1, " + ZZ", i,i-1, "*w", Z,i-1, sep=""))
      } # end (if Z/W)
    } # end (for i)

    if (lag == 2) {
      cat(rep("\n",2), "  # -- Estimate lagged effects between within-person centered variables (Lag = 2 waves) -- #")
      for (i in 3:no.waves) {
        if (Z == "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-2, "*w", X,i-2, " + YX", i,i-2, "*w", Y,i-2, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-2, "*w", X,i-2, " + YY", i,i-2, "*w", Y,i-2, sep=""))
        } else if (W != "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-2, "*w", X,i-2, " + YX", i,i-2, "*w", Y,i-2, " + ZX", i,i-2, "*w", Z,i-2, " + WX", i,i-2, "*w", W,i-2, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-2, "*w", X,i-2, " + YY", i,i-2, "*w", Y,i-2, " + ZY", i,i-2, "*w", Z,i-2, " + WY", i,i-2, "*w", W,i-2, sep=""))
          cat("\n", paste("  w", Z,i, " ~ XZ", i,i-2, "*w", X,i-2, " + YZ", i,i-2, "*w", Y,i-2, " + ZZ", i,i-2, "*w", Z,i-2, " + WZ", i,i-2, "*w", W,i-2, sep=""))
          cat("\n", paste("  w", W,i, " ~ XW", i,i-2, "*w", X,i-2, " + YW", i,i-2, "*w", Y,i-2, " + ZW", i,i-2, "*w", Z,i-2, " + WW", i,i-2, "*w", W,i-2, sep=""))
        } else if (Z != "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-2, "*w", X,i-2, " + YX", i,i-2, "*w", Y,i-2, " + ZX", i,i-2, "*w", Z,i-2, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-2, "*w", X,i-2, " + YY", i,i-2, "*w", Y,i-2, " + ZY", i,i-2, "*w", Z,i-2, sep=""))
          cat("\n", paste("  w", Z,i, " ~ XZ", i,i-2, "*w", X,i-2, " + YZ", i,i-2, "*w", Y,i-2, " + ZZ", i,i-2, "*w", Z,i-2, sep=""))
        } # end (if Z/W)
      } # end (for i)
    } # end (if lag == 2)

    if (lag == 3) {
      cat(rep("\n",2), "  # -- Estimate lagged effects between within-person centered variables (Lag = 3 waves) -- #")
      for (i in 4:no.waves) {
        if (Z == "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-3, "*w", X,i-3, " + YX", i,i-3, "*w", Y,i-3, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-3, "*w", X,i-3, " + YY", i,i-3, "*w", Y,i-3, sep=""))
        } else if (W != "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-3, "*w", X,i-3, " + YX", i,i-3, "*w", Y,i-3, " + ZX", i,i-3, "*w", Z,i-3, " + WX", i,i-3, "*w", W,i-3, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-3, "*w", X,i-3, " + YY", i,i-3, "*w", Y,i-3, " + ZY", i,i-3, "*w", Z,i-3, " + WY", i,i-3, "*w", W,i-3, sep=""))
          cat("\n", paste("  w", Z,i, " ~ XZ", i,i-3, "*w", X,i-3, " + YZ", i,i-3, "*w", Y,i-3, " + ZZ", i,i-3, "*w", Z,i-3, " + WZ", i,i-3, "*w", W,i-3, sep=""))
          cat("\n", paste("  w", W,i, " ~ XW", i,i-3, "*w", X,i-3, " + YW", i,i-3, "*w", Y,i-3, " + ZW", i,i-3, "*w", Z,i-3, " + WW", i,i-3, "*w", W,i-3, sep=""))
        } else if (Z != "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-3, "*w", X,i-3, " + YX", i,i-3, "*w", Y,i-3, " + ZX", i,i-3, "*w", Z,i-3, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-3, "*w", X,i-3, " + YY", i,i-3, "*w", Y,i-3, " + ZY", i,i-3, "*w", Z,i-3, sep=""))
          cat("\n", paste("  w", Z,i, " ~ XZ", i,i-3, "*w", X,i-3, " + YZ", i,i-3, "*w", Y,i-3, " + ZZ", i,i-3, "*w", Z,i-3, sep=""))
        } # end (if Z/W)
      } # end (for i)
    } # end (lag == 3)

    if (lag == 4) {
      cat(rep("\n",2), "  # -- Estimate lagged effects between within-person centered variables (Lag = 4 waves) -- #")
      for (i in 5:no.waves) {
        if (Z == "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-4, "*w", X,i-4, " + YX", i,i-4, "*w", Y,i-4, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-4, "*w", X,i-4, " + YY", i,i-4, "*w", Y,i-4, sep=""))
        } else if (W != "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-4, "*w", X,i-4, " + YX", i,i-4, "*w", Y,i-4, " + ZX", i,i-4, "*w", Z,i-4, " + WX", i,i-4, "*w", W,i-4, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-4, "*w", X,i-4, " + YY", i,i-4, "*w", Y,i-4, " + ZY", i,i-4, "*w", Z,i-4, " + WY", i,i-4, "*w", W,i-4, sep=""))
          cat("\n", paste("  w", Z,i, " ~ XZ", i,i-4, "*w", X,i-4, " + YZ", i,i-4, "*w", Y,i-4, " + ZZ", i,i-4, "*w", Z,i-4, " + WZ", i,i-4, "*w", W,i-4, sep=""))
          cat("\n", paste("  w", W,i, " ~ XW", i,i-4, "*w", X,i-4, " + YW", i,i-4, "*w", Y,i-4, " + ZW", i,i-4, "*w", Z,i-4, " + WW", i,i-4, "*w", W,i-4, sep=""))
        } else if (Z != "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-4, "*w", X,i-4, " + YX", i,i-4, "*w", Y,i-4, " + ZX", i,i-4, "*w", Z,i-4, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-4, "*w", X,i-4, " + YY", i,i-4, "*w", Y,i-4, " + ZY", i,i-4, "*w", Z,i-4, sep=""))
          cat("\n", paste("  w", Z,i, " ~ XZ", i,i-4, "*w", X,i-4, " + YZ", i,i-4, "*w", Y,i-4, " + ZZ", i,i-4, "*w", Z,i-4, sep=""))
        } # end (if Z/W)
      } # end (for i)
    } # end (if lag == 4)

    cat(rep("\n",2), "  '")

    # -- Run Model LGCMSRMLR -- #
    cat(rep("\n",2), "# -- Run Model LGCMSRMLR -- #")
    cat(rep("\n",2), "  LGCMSRMLR.fit <- suppressWarnings(lavaan::sem(LGCMSR,")
    cat("\n", "   ", data.source, ",")
    cat("\n", "   missing = 'fiml',")
    cat("\n", "   meanstructure = TRUE,")
    cat("\n", "   information = 'observed',")
    cat("\n", "   estimator = 'MLR'")
    cat("\n", "))")

  sink()  # Stop writing to file LGCMSR.txt
  ## ------------------------------- ##


  ## -- Execute LGCMSR.txt and request summary outputs-- ##
  source('LGCMSR.txt')
  print(lavaan::summary(LGCMSRMLR.fit, fit.measure = TRUE, standardized = TRUE, rsq = TRUE))
  if (lavaan::lavInspect(LGCMSRMLR.fit, what ="post.check") == FALSE) stop("The lavaan solution is non-admissible.")
  ## ------------------------------- ##



  ## -- Monte Carlo Simulation -- ##

  parEst <- lavaan::parameterEstimates(LGCMSRMLR.fit, remove.nonfree = TRUE)
  pest2 <- parEst[,5]  # Estimated Parameters
  pest3 <- lavaan::lavTech(LGCMSRMLR.fit, what = "vcov", add.labels = TRUE)  # Estimated Variance-Covariance of Estimated Parameters
 
  varI.eq = TRUE
  Invariance(parEst, pest2, pest3, no.path, MIset, no.compare, no.waves, lag, varI.eq, p, X, Y, Z, W)

  cat(rep("\n", 2))

  ## ----- Start writing script to LGCMSR.txt ----- ##
  sink('LGCMSR.txt')
    cat("\n", "# Specify the model (LGCMSR)", "\n")
    cat("\n", "LGCMSR <- '")

    cat(rep("\n",2), "  # -- Create between components (random intercepts) -- #")
    BX <- paste("  RI", X, " =~ 1*", X, "1", sep="")
    BY <- paste("  RI", Y, " =~ 1*", Y, "1", sep="")
    if (Z != "NULL") {
      BY <- paste("  RI", Z, " =~ 1*", Z, "1", sep="")
    } # end (if Z != "NULL")
    if (W != "NULL") {
      BY <- paste("  RI", W, " =~ 1*", W, "1", sep="")
    } # end (if W != "NULL")

    for (i in 2:no.waves) {
      BX <- paste(BX, " +1*", X, i, sep="")
      BY <- paste(BY, " +1*", Y, i, sep="")
      if (Z != "NULL") {
        BZ <- paste(BZ, " +1*", Z, i, sep="")
      } # end (if Z != "NULL")
      if (W != "NULL") {
        BW <- paste(BW, " +1*", W, i, sep="")
      } # end (if W != "NULL")
    } # end (for i)

    cat("\n", BX)
    cat("\n", BY)
    if (Z != "NULL") {
      cat("\n", BZ)
    } # end (if Z != "NULL")
    if (W != "NULL") {
      cat("\n", BW)
    } # end (if W != "NULL")


    cat(rep("\n",2), "  # -- Create between components (random slopes) -- #")
    BX <- paste("  RS", X, " =~ 0*", X, "1", sep="")
    BY <- paste("  RS", Y, " =~ 0*", Y, "1", sep="")
    if (Z != "NULL") {
      BY <- paste("  RS", Z, " =~ 0*", Z, "1", sep="")
    } # end (if Z != "NULL")
    if (W != "NULL") {
      BY <- paste("  RS", W, " =~ 0*", W, "1", sep="")
    } # end (if W != "NULL")

    for (i in 2:no.waves) {
      BX <- paste(BX, " +", (i-1), "*", X, i, sep="")
      BY <- paste(BY, " +", (i-1), "*", Y, i, sep="")
      if (Z != "NULL") {
        BZ <- paste(BZ, " +", (i-1), "*", Z, i, sep="")
      } # end (if Z != "NULL")
      if (W != "NULL") {
        BW <- paste(BW, " +", (i-1), "*", W, i, sep="")
      } # end (if W != "NULL")
    } # end (for i)

    cat("\n", BX)
    cat("\n", BY)
    if (Z != "NULL") {
      cat("\n", BZ)
    } # end (if Z != "NULL")
    if (W != "NULL") {
      cat("\n", BW)
    } # end (if W != "NULL")

    # -- Constrain Residual Variance of Observed Variables to Zero -- #
    cat(rep("\n",2), "  # -- Constrain residual variance of observed variables to zero -- #")
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

    # -- Create Latent Variables from Observed Variables -- #
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

    cat(rep("\n",2), "  # -- Estimate lagged effects between latent variables (Lag = 1 wave) -- #")
    cat("\n", "  #############################################")
    cat("\n", "  # Remove the subscripts for invariant paths #")
    cat("\n", "  #############################################")
    for (i in 2:no.waves) {
      if (Z == "NULL") {
        cat("\n", paste("  w", X,i, " ~ XX", i,i-1, "*w", X,i-1, " + YX", i,i-1, "*w", Y,i-1, sep=""))
        cat("\n", paste("  w", Y,i, " ~ XY", i,i-1, "*w", X,i-1, " + YY", i,i-1, "*w", Y,i-1, sep=""))
      } else if (W != "NULL") {
        cat("\n", paste("  w", X,i, " ~ XX", i,i-1, "*w", X,i-1, " + YX", i,i-1, "*w", Y,i-1, " + ZX", i,i-1, "*w", Z,i-1, " + WX", i,i-1, "*w", W,i-1, sep=""))
        cat("\n", paste("  w", Y,i, " ~ XY", i,i-1, "*w", X,i-1, " + YY", i,i-1, "*w", Y,i-1, " + ZY", i,i-1, "*w", Z,i-1, " + WY", i,i-1, "*w", W,i-1, sep=""))
        cat("\n", paste("  w", Z,i, " ~ XZ", i,i-1, "*w", X,i-1, " + YZ", i,i-1, "*w", Y,i-1, " + ZZ", i,i-1, "*w", Z,i-1, " + WZ", i,i-1, "*w", W,i-1, sep=""))
        cat("\n", paste("  w", W,i, " ~ XW", i,i-1, "*w", X,i-1, " + YW", i,i-1, "*w", Y,i-1, " + ZW", i,i-1, "*w", Z,i-1, " + WW", i,i-1, "*w", W,i-1, sep=""))
      } else if (Z != "NULL") {
        cat("\n", paste("  w", X,i, " ~ XX", i,i-1, "*w", X,i-1, " + YX", i,i-1, "*w", Y,i-1, " + ZX", i,i-1, "*w", Z,i-1, sep=""))
        cat("\n", paste("  w", Y,i, " ~ XY", i,i-1, "*w", X,i-1, " + YY", i,i-1, "*w", Y,i-1, " + ZY", i,i-1, "*w", Z,i-1, sep=""))
        cat("\n", paste("  w", Z,i, " ~ XZ", i,i-1, "*w", X,i-1, " + YZ", i,i-1, "*w", Y,i-1, " + ZZ", i,i-1, "*w", Z,i-1, sep=""))
      } # end (if Z)
    } # end (for i)

    if (lag == 2) {
      cat(rep("\n",2), "  # -- Estimate lagged effects between latent variables (Lag = 2 waves) -- #")
      cat("\n", "  #############################################")
      cat("\n", "  # Remove the subscripts for invariant paths #")
      cat("\n", "  #############################################")
      for (i in 3:no.waves) {
        if (Z == "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-2, "*w", X,i-2, " + YX", i,i-2, "*w", Y,i-2, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-2, "*w", X,i-2, " + YY", i,i-2, "*w", Y,i-2, sep=""))
        } else if (W != "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-2, "*w", X,i-2, " + YX", i,i-2, "*w", Y,i-2, " + ZX", i,i-2, "*w", Z,i-2, " + WX", i,i-2, "*w", W,i-2, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-2, "*w", X,i-2, " + YY", i,i-2, "*w", Y,i-2, " + ZY", i,i-2, "*w", Z,i-2, " + WY", i,i-2, "*w", W,i-2, sep=""))
          cat("\n", paste("  w", Z,i, " ~ XZ", i,i-2, "*w", X,i-2, " + YZ", i,i-2, "*w", Y,i-2, " + ZZ", i,i-2, "*w", Z,i-2, " + WZ", i,i-2, "*w", W,i-2, sep=""))
          cat("\n", paste("  w", W,i, " ~ XW", i,i-2, "*w", X,i-2, " + YW", i,i-2, "*w", Y,i-2, " + ZW", i,i-2, "*w", Z,i-2, " + WW", i,i-2, "*w", W,i-2, sep=""))
        } else if (Z != "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-2, "*w", X,i-2, " + YX", i,i-2, "*w", Y,i-2, " + ZX", i,i-2, "*w", Z,i-2, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-2, "*w", X,i-2, " + YY", i,i-2, "*w", Y,i-2, " + ZY", i,i-2, "*w", Z,i-2, sep=""))
          cat("\n", paste("  w", Z,i, " ~ XZ", i,i-2, "*w", X,i-2, " + YZ", i,i-2, "*w", Y,i-2, " + ZZ", i,i-2, "*w", Z,i-2, sep=""))
        } # end (if Z/W)
      } # end (for i)
    } # end (if lag == 2)

    if (lag == 3) {
      cat(rep("\n",2), "  # -- Estimate lagged effects between latent variables (Lag = 3 waves) -- #")
      cat("\n", "  #############################################")
      cat("\n", "  # Remove the subscripts for invariant paths #")
      cat("\n", "  #############################################")
      for (i in 3:no.waves) {
        if (Z == "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-3, "*w", X,i-3, " + YX", i,i-3, "*w", Y,i-3, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-3, "*w", X,i-3, " + YY", i,i-3, "*w", Y,i-3, sep=""))
        } else if (W != "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-3, "*w", X,i-3, " + YX", i,i-3, "*w", Y,i-3, " + ZX", i,i-3, "*w", Z,i-3, " + WX", i,i-3, "*w", W,i-3, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-3, "*w", X,i-3, " + YY", i,i-3, "*w", Y,i-3, " + ZY", i,i-3, "*w", Z,i-3, " + WY", i,i-3, "*w", W,i-3, sep=""))
          cat("\n", paste("  w", Z,i, " ~ XZ", i,i-3, "*w", X,i-3, " + YZ", i,i-3, "*w", Y,i-3, " + ZZ", i,i-3, "*w", Z,i-3, " + WZ", i,i-3, "*w", W,i-3, sep=""))
          cat("\n", paste("  w", W,i, " ~ XW", i,i-3, "*w", X,i-3, " + YW", i,i-3, "*w", Y,i-3, " + ZW", i,i-3, "*w", Z,i-3, " + WW", i,i-3, "*w", W,i-3, sep=""))
        } else if (Z != "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-3, "*w", X,i-3, " + YX", i,i-3, "*w", Y,i-3, " + ZX", i,i-3, "*w", Z,i-3, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-3, "*w", X,i-3, " + YY", i,i-3, "*w", Y,i-3, " + ZY", i,i-3, "*w", Z,i-3, sep=""))
          cat("\n", paste("  w", Z,i, " ~ XZ", i,i-3, "*w", X,i-3, " + YZ", i,i-3, "*w", Y,i-3, " + ZZ", i,i-3, "*w", Z,i-3, sep=""))
        } # end (if Z/W)
      } # end (for i)
    } # end (lag == 3)

   if (lag == 4) {
      cat(rep("\n",2), "  # -- Estimate lagged effects between latent variables (Lag = 4 waves) -- #")
      cat("\n", "  #############################################")
      cat("\n", "  # Remove the subscripts for invariant paths #")
      cat("\n", "  #############################################")
      for (i in 3:no.waves) {
        if (Z == "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-4, "*w", X,i-4, " + YX", i,i-4, "*w", Y,i-4, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-4, "*w", X,i-4, " + YY", i,i-4, "*w", Y,i-4, sep=""))
        } else if (W != "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-4, "*w", X,i-4, " + YX", i,i-4, "*w", Y,i-4, " + ZX", i,i-4, "*w", Z,i-4, " + WX", i,i-4, "*w", W,i-4, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-4, "*w", X,i-4, " + YY", i,i-4, "*w", Y,i-4, " + ZY", i,i-4, "*w", Z,i-4, " + WY", i,i-4, "*w", W,i-4, sep=""))
          cat("\n", paste("  w", Z,i, " ~ XZ", i,i-4, "*w", X,i-4, " + YZ", i,i-4, "*w", Y,i-4, " + ZZ", i,i-4, "*w", Z,i-4, " + WZ", i,i-4, "*w", W,i-4, sep=""))
          cat("\n", paste("  w", W,i, " ~ XW", i,i-4, "*w", X,i-4, " + YW", i,i-4, "*w", Y,i-4, " + ZW", i,i-4, "*w", Z,i-4, " + WW", i,i-4, "*w", W,i-4, sep=""))
        } else if (Z != "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-4, "*w", X,i-4, " + YX", i,i-4, "*w", Y,i-4, " + ZX", i,i-4, "*w", Z,i-4, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-4, "*w", X,i-4, " + YY", i,i-4, "*w", Y,i-4, " + ZY", i,i-4, "*w", Z,i-4, sep=""))
          cat("\n", paste("  w", Z,i, " ~ XZ", i,i-4, "*w", X,i-4, " + YZ", i,i-4, "*w", Y,i-4, " + ZZ", i,i-4, "*w", Z,i-4, sep=""))
        } # end (if Z/W)
      } # end (for i)
    } # end (lag == 4)


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


    cat(rep("\n",2), "  # -- Estimate (residual) variance of latent variables -- #")
    for (i in 1:no.waves) {
      cat("\n", paste("  w", X, i, " ~~ ", "w", X, i, sep=""))
      cat("\n", paste("  w", Y, i, " ~~ ", "w", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  w", Z, i, " ~~ ", "w", Z, i, sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  w", W, i, " ~~ ", "w", W, i, sep=""))
      } # end (if W)
    } # end (for i)


    cat(rep("\n",2), "  # -- Estimate grand means (intercepts) of observed variables -- #")
    for (i in 1:no.waves) {
      cat("\n", paste("  ", X, i, " ~ M", X, i, "*1", sep=""))
      cat("\n", paste("  ", Y, i, " ~ M", Y, i, "*1", sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  ", Z, i, " ~ M", Z, i, "*1", sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  ", W, i, " ~ M", W, i, "*1", sep=""))
      } # end (if W)
    } # end (for i)


    # -- Constrain covariance among RIX RIY, RIZ, RIW and wx1, wy1, wz1, ww1 -- #
    cat(rep("\n",2), "  # -- Constrain covariance among RIx, RIy, RIz, RIw and wx1, wy1, wz1, ww1 -- #")
    cat("\n", "   RI", X, " ~~ 0*w", X, "1", sep="")
    cat("\n", "   RI", X, " ~~ 0*w", Y, "1", sep="")
    if (Z != "NULL") {
      cat("\n", "   RI", X, " ~~ 0*w", Z, "1", sep="")
    } # end (if Z)
    if (W != "NULL") {
      cat("\n", "   RI", X, " ~~ 0*w", W, "1", sep="")
    } # end (if W)
    cat("\n", "   RI", Y, " ~~ 0*w", X, "1", sep="")
    cat("\n", "   RI", Y, " ~~ 0*w", Y, "1", sep="")
    if (Z != "NULL") {
      cat("\n", "   RI", Y, " ~~ 0*w", Z, "1", sep="")
    } # end (if Z)
    if (W != "NULL") {
      cat("\n", "   RI", Y, " ~~ 0*w", W, "1", sep="")
    } # end (if W)
    if (Z != "NULL") {
      cat("\n", "   RI", Z, " ~~ 0*w", X, "1", sep="")
      cat("\n", "   RI", Z, " ~~ 0*w", Y, "1", sep="")
      cat("\n", "   RI", Z, " ~~ 0*w", Z, "1", sep="")
    } # end (if Z)
    if (W != "NULL") {
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
    if (Z != "NULL") {
      cat("\n", "   RS", X, " ~~ 0*w", Z, "1", sep="")
    } # end (if Z != "NULL")
    if (W != "NULL") {
      cat("\n", "   RS", X, " ~~ 0*w", W, "1", sep="")
    } # end (if W != "NULL")

    cat("\n", "   RS", Y, " ~~ 0*w", X, "1", sep="")
    cat("\n", "   RS", Y, " ~~ 0*w", Y, "1", sep="")
    if (Z != "NULL") {
      cat("\n", "   RS", Y, " ~~ 0*w", Z, "1", sep="")
    } # end (if Z != "NULL")
    if (W != "NULL") {
      cat("\n", "   RS", Y, " ~~ 0*w", W, "1", sep="")
    } # end (if W != "NULL")
    if (Z != "NULL") {
      cat("\n", "   RS", Z, " ~~ 0*w", X, "1", sep="")
      cat("\n", "   RS", Z, " ~~ 0*w", Y, "1", sep="")
      cat("\n", "   RS", Z, " ~~ 0*w", Z, "1", sep="")
    } # end (if Z != "NULL")
    if (W != "NULL") {
      cat("\n", "   RS", Z, " ~~ 0*w", W, "1", sep="")
      cat("\n", "   RS", W, " ~~ 0*w", X, "1", sep="")
      cat("\n", "   RS", W, " ~~ 0*w", Y, "1", sep="")
      cat("\n", "   RS", W, " ~~ 0*w", Z, "1", sep="")
      cat("\n", "   RS", W, " ~~ 0*w", W, "1", sep="")
    } # end (if W != "NULL")


    cat(rep("\n",2), "  ##########################################")
    cat("\n", "  # Regression of observed variables on C1 #")
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
    cat(rep("\n",2), "  print(lavaan::summary(LGCMSRMLR.fit, fit.measure = TRUE, standardized = TRUE, rsq = TRUE))", "\n")

    cat(rep("\n",2))

  sink() # Stop writing to file

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
#' @param p critical p-value for pairwise comparisons (default is 0.001).
#' @param X name of variable X.
#' @param Y name of variable Y.
#' @param Z name of variable Z.
#' @param W name of variable W (Z must come before W).
#'
#' @return STARTS outputs.
#' @export
#' @examples
#'
#' ## -- Example -- ##
#'
#' STARTS(data.source="Data_A", 7, 2, X="EXPOSE", Y="INTENS")
#'

STARTS <- function(data.source, no.waves, lag=1, varI.eq = FALSE, p = 0.001, X, Y, Z="NULL", W = "NULL") {

  ## -- Check inputs -- ##

  if (no.waves < 3) stop("Minimum number of waves is 3")

  if (lag < 1) stop("Minimum number of lag is 1")
  if (lag > 4) stop("Maximum number of lag is 4")

  if ((no.waves - lag) < 1) stop("Number of waves must be greater than lag plus 1")

  if (p > 0.05) stop("p > 0.05 is not recommended")
  if (p < 0.0001) stop("p < 0.0001 is not recommended")

  if (Z == "NULL" & W != "NULL") stop("Z must be defined before W")

  if (is.logical(varI.eq) == FALSE) stop("varI.eq (equivalence of indicator residual variance) can only be TRUE or FALSE")


  ## ----- Creating Model STARTS ----- ###
  sink('STARTS.txt') # Start writing script to STARTS.txt

    cat("\n", "## ----- Specify the model (STARTS) ----- ##", "\n")
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


    # -- Create Residual Variance of Observed Variables -- #
    if (isFALSE(varI.eq)) {
      cat(rep("\n",2), "  # -- Create residual variance of observed variables -- #")
      for (i in 1:no.waves) {
        cat("\n", paste("  ", X, i, " ~~ varI", X, i, "*", X, i, sep=""))
        cat("\n", paste("  ", Y, i, " ~~ varI", Y, i, "*", Y, i, sep=""))
        if (Z != "NULL") {
          cat("\n", paste("  ", Z, i, " ~~ varI", Z, i, "*", Z, i, sep=""))
        } # end (if Z)
        if (W != "NULL") {
          cat("\n", paste("  ", W, i, " ~~ varI", W, i, "*", W, i, sep=""))
        } # end (if W)
      } # end (for i)###   ###
    } else {
       cat(rep("\n",2), "  # -- Create residual variance of observed variables -- #")
      for (i in 1:no.waves) {
        cat("\n", paste("  ", X, i, " ~~ varI", X, "*", X, i, sep=""))
        cat("\n", paste("  ", Y, i, " ~~ varI", Y, "*", Y, i, sep=""))
        if (Z != "NULL") {
          cat("\n", paste("  ", Z, i, " ~~ varI", Z, "*", Z, i, sep=""))
        } # end (if Z)
        if (W != "NULL") {
          cat("\n", paste("  ", W, i, " ~~ varI", W, "*", W, i, sep=""))
        } # end (if W)
      } # end (for i)###   ###
    } # end (if varI)

    # -- Create Latent Variables from Observed Variables -- #
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

    # -- Estimate Covariances Between Residuals of Latent Variables -- #
    cat(rep("\n",2), "  # -- Estimate covariances between residuals of latent variables -- #")
    for (i in 2:no.waves) {
      cat("\n", paste("  w", X, i, " ~~ eXY", i, "*w", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  w", X, i, " ~~ eXZ", i, "*w", Z, i, sep=""))
        cat("\n", paste("  w", Y, i, " ~~ eYZ", i, "*w", Z, i, sep=""))
      }  # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  w", X, i, " ~~ eXW", i, "*w", W, i, sep=""))
        cat("\n", paste("  w", Y, i, " ~~ eYW", i, "*w", W, i, sep=""))
        cat("\n", paste("  w", Z, i, " ~~ eZW", i, "*w", W, i, sep=""))
      }  # end (if W)
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


    # -- Estimate Variance and Covariance of Random Intercepts -- #
    cat(rep("\n",2), "  # -- Estimate variance and covariance of random intercepts -- #")
    cat("\n", "   RI", X, " ~~ RI", X, sep="")
    cat("\n", "   RI", Y, " ~~ RI", Y, sep="")
    if (Z != "NULL") {
      cat("\n", "   RI", Z, " ~~ RI", Z, sep="")
    } # end (if Z != "NULL")
    if (W != "NULL") {
      cat("\n", "   RI", W, " ~~ RI", W, sep="")
    } # end (if W != "NULL")

    cat("\n", "   RI", X, " ~~ RI", Y, sep="")
    if (Z != "NULL") {
      cat("\n", "   RI", X, " ~~ RI", Z, sep="")
      cat("\n", "   RI", Y, " ~~ RI", Z, sep="")
    } # end (if Z != "NULL")
    if (W != "NULL") {
      cat("\n", "   RI", X, " ~~ RI", W, sep="")
      cat("\n", "   RI", Y, " ~~ RI", W, sep="")
      cat("\n", "   RI", Z, " ~~ RI", W, sep="")
    } # end (if W != "NULL")


    # -- Estimate (Residual) Variance of Latent Variables -- #
    cat(rep("\n",2), "  # -- Estimate (residual) variance of latent variables -- #")
    for (i in 1:no.waves) {
      cat("\n", paste("  w", X, i, " ~~ ", "w", X, i, sep=""))
      cat("\n", paste("  w", Y, i, " ~~ ", "w", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  w", Z, i, " ~~ ", "w", Z, i, sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  w", W, i, " ~~ ", "w", W, i, sep=""))
      } # end (if W)
    } # end (for i)

    # -- Estimate Grand Means (Intercepts) of Observed Variables -- #
    cat(rep("\n",2), "  # -- Estimate grand means (intercepts) of observed variables -- #")
    for (i in 1:no.waves) {
      cat("\n", paste("  ", X, i, " ~ M", X, i, "*1", sep=""))
      cat("\n", paste("  ", Y, i, " ~ M", Y, i, "*1", sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  ", Z, i, " ~ M", Z, i, "*1", sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  ", W, i, " ~ M", W, i, "*1", sep=""))
      } # (if W)
    } # end (for i)

    # -- Constrain covariance among RIX RIY, RIZ, RIW and wx1, wy1, wz1, ww1 -- #
    cat(rep("\n",2), "  # -- Constrain covariance among RIx, RIy, RIz, RIw and wx1, wy1, wz1, ww1 -- #")
    cat("\n", "   RI", X, " ~~ 0*w", X, "1", sep="")
    cat("\n", "   RI", X, " ~~ 0*w", Y, "1", sep="")
    if (Z != "NULL") {
      cat("\n", "   RI", X, " ~~ 0*w", Z, "1", sep="")
    } # end (if Z != "NULL")
    if (W != "NULL") {
      cat("\n", "   RI", X, " ~~ 0*w", W, "1", sep="")
    } # end (if W != "NULL")

    cat("\n", "   RI", Y, " ~~ 0*w", X, "1", sep="")
    cat("\n", "   RI", Y, " ~~ 0*w", Y, "1", sep="")
    if (Z != "NULL") {
      cat("\n", "   RI", Y, " ~~ 0*w", Z, "1", sep="")
    } # end (if Z != "NULL")
    if (W != "NULL") {
      cat("\n", "   RI", Y, " ~~ 0*w", W, "1", sep="")
    } # end (if W != "NULL")
    if (Z != "NULL") {
      cat("\n", "   RI", Z, " ~~ 0*w", X, "1", sep="")
      cat("\n", "   RI", Z, " ~~ 0*w", Y, "1", sep="")
      cat("\n", "   RI", Z, " ~~ 0*w", Z, "1", sep="")
    } # end (if Z != "NULL")
    if (W != "NULL") {
      cat("\n", "   RI", Z, " ~~ 0*w", W, "1", sep="")
      cat("\n", "   RI", W, " ~~ 0*w", X, "1", sep="")
      cat("\n", "   RI", W, " ~~ 0*w", Y, "1", sep="")
      cat("\n", "   RI", W, " ~~ 0*w", Z, "1", sep="")
      cat("\n", "   RI", W, " ~~ 0*w", W, "1", sep="")
    } # end (if W != "NULL")

    # -- Estimate Lagged Effects Between Latent Variables -- #
    cat(rep("\n",2), "  # -- Estimate lagged effects between latent variables (Lag = 1 wave) -- #")
    for (i in 2:no.waves) {
      if (Z == "NULL") {
        cat("\n", paste("  w", X, i, " ~ XX", i,i-1, "*w", X,i-1, " + YX", i,i-1, "*w", Y,i-1, sep=""))
        cat("\n", paste("  w", Y, i, " ~ XY", i,i-1, "*w", X,i-1, " + YY", i,i-1, "*w", Y,i-1, sep=""))
      } else if (W != "NULL") {
        cat("\n", paste("  w", X, i, " ~ XX", i,i-1, "*w", X,i-1, " + YX", i,i-1, "*w", Y,i-1, " + ZX", i,i-1, "*w", Z,i-1, " + WX", i,i-1, "*w", W,i-1, sep=""))
        cat("\n", paste("  w", Y, i, " ~ XY", i,i-1, "*w", X,i-1, " + YY", i,i-1, "*w", Y,i-1, " + ZY", i,i-1, "*w", Z,i-1, " + WY", i,i-1, "*w", W,i-1, sep=""))
        cat("\n", paste("  w", Z, i, " ~ XZ", i,i-1, "*w", X,i-1, " + YZ", i,i-1, "*w", Y,i-1, " + ZZ", i,i-1, "*w", Z,i-1, " + WZ", i,i-1, "*w", W,i-1, sep=""))
        cat("\n", paste("  w", W, i, " ~ XW", i,i-1, "*w", X,i-1, " + YW", i,i-1, "*w", Y,i-1, " + ZW", i,i-1, "*w", Z,i-1, " + WW", i,i-1, "*w", W,i-1, sep=""))
      } else if (Z != "NULL") {
        cat("\n", paste("  w", X, i, " ~ XX", i,i-1, "*w", X,i-1, " + YX", i,i-1, "*w", Y,i-1, " + ZX", i,i-1, "*w", Z,i-1, sep=""))
        cat("\n", paste("  w", Y, i, " ~ XY", i,i-1, "*w", X,i-1, " + YY", i,i-1, "*w", Y,i-1, " + ZY", i,i-1, "*w", Z,i-1, sep=""))
        cat("\n", paste("  w", Z, i, " ~ XZ", i,i-1, "*w", X,i-1, " + YZ", i,i-1, "*w", Y,i-1, " + ZZ", i,i-1, "*w", Z,i-1, sep=""))
      } # end (if Z/W)
    } # end (for i)

    if (lag == 2) {
      cat(rep("\n",2), "  # -- Estimate lagged effects between within-person centered variables (Lag = 2 waves) -- #")
      for (i in 3:no.waves) {
        if (Z == "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-2, "*w", X,i-2, " + YX", i,i-2, "*w", Y,i-2, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-2, "*w", X,i-2, " + YY", i,i-2, "*w", Y,i-2, sep=""))
        } else if (W != "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-2, "*w", X,i-2, " + YX", i,i-2, "*w", Y,i-2, " + ZX", i,i-2, "*w", Z,i-2, " + WX", i,i-2, "*w", W,i-2, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-2, "*w", X,i-2, " + YY", i,i-2, "*w", Y,i-2, " + ZY", i,i-2, "*w", Z,i-2, " + WY", i,i-2, "*w", W,i-2, sep=""))
          cat("\n", paste("  w", Z,i, " ~ XZ", i,i-2, "*w", X,i-2, " + YZ", i,i-2, "*w", Y,i-2, " + ZZ", i,i-2, "*w", Z,i-2, " + WZ", i,i-2, "*w", W,i-2, sep=""))
          cat("\n", paste("  w", W,i, " ~ XW", i,i-2, "*w", X,i-2, " + YW", i,i-2, "*w", Y,i-2, " + ZW", i,i-2, "*w", Z,i-2, " + WW", i,i-2, "*w", W,i-2, sep=""))
        } else if (Z != "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-2, "*w", X,i-2, " + YX", i,i-2, "*w", Y,i-2, " + ZX", i,i-2, "*w", Z,i-2, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-2, "*w", X,i-2, " + YY", i,i-2, "*w", Y,i-2, " + ZY", i,i-2, "*w", Z,i-2, sep=""))
          cat("\n", paste("  w", Z,i, " ~ XZ", i,i-2, "*w", X,i-2, " + YZ", i,i-2, "*w", Y,i-2, " + ZZ", i,i-2, "*w", Z,i-2, sep=""))
        } # end (if Z/W)
      } # end (for i)
    } # end (if lag == 2)

    if (lag == 3) {
      cat(rep("\n",2), "  # -- Estimate lagged effects between within-person centered variables (Lag = 3 waves) -- #")
      for (i in 4:no.waves) {
        if (Z == "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-3, "*w", X,i-3, " + YX", i,i-3, "*w", Y,i-3, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-3, "*w", X,i-3, " + YY", i,i-3, "*w", Y,i-3, sep=""))
        } else if (W != "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-3, "*w", X,i-3, " + YX", i,i-3, "*w", Y,i-3, " + ZX", i,i-3, "*w", Z,i-3, " + WX", i,i-3, "*w", W,i-3, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-3, "*w", X,i-3, " + YY", i,i-3, "*w", Y,i-3, " + ZY", i,i-3, "*w", Z,i-3, " + WY", i,i-3, "*w", W,i-3, sep=""))
          cat("\n", paste("  w", Z,i, " ~ XZ", i,i-3, "*w", X,i-3, " + YZ", i,i-3, "*w", Y,i-3, " + ZZ", i,i-3, "*w", Z,i-3, " + WZ", i,i-3, "*w", W,i-3, sep=""))
          cat("\n", paste("  w", W,i, " ~ XW", i,i-3, "*w", X,i-3, " + YW", i,i-3, "*w", Y,i-3, " + ZW", i,i-3, "*w", Z,i-3, " + WW", i,i-3, "*w", W,i-3, sep=""))
        } else if (Z != "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-3, "*w", X,i-3, " + YX", i,i-3, "*w", Y,i-3, " + ZX", i,i-3, "*w", Z,i-3, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-3, "*w", X,i-3, " + YY", i,i-3, "*w", Y,i-3, " + ZY", i,i-3, "*w", Z,i-3, sep=""))
          cat("\n", paste("  w", Z,i, " ~ XZ", i,i-3, "*w", X,i-3, " + YZ", i,i-3, "*w", Y,i-3, " + ZZ", i,i-3, "*w", Z,i-3, sep=""))
        } # end (if Z/W)
      } # end (for i)
    } # end (lag == 3)

    if (lag == 4) {
      cat(rep("\n",2), "  # -- Estimate lagged effects between within-person centered variables (Lag = 4 waves) -- #")
      for (i in 5:no.waves) {
        if (Z == "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-4, "*w", X,i-4, " + YX", i,i-4, "*w", Y,i-4, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-4, "*w", X,i-4, " + YY", i,i-4, "*w", Y,i-4, sep=""))
        } else if (W != "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-4, "*w", X,i-4, " + YX", i,i-4, "*w", Y,i-4, " + ZX", i,i-4, "*w", Z,i-4, " + WX", i,i-4, "*w", W,i-4, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-4, "*w", X,i-4, " + YY", i,i-4, "*w", Y,i-4, " + ZY", i,i-4, "*w", Z,i-4, " + WY", i,i-4, "*w", W,i-4, sep=""))
          cat("\n", paste("  w", Z,i, " ~ XZ", i,i-4, "*w", X,i-4, " + YZ", i,i-4, "*w", Y,i-4, " + ZZ", i,i-4, "*w", Z,i-4, " + WZ", i,i-4, "*w", W,i-4, sep=""))
          cat("\n", paste("  w", W,i, " ~ XW", i,i-4, "*w", X,i-4, " + YW", i,i-4, "*w", Y,i-4, " + ZW", i,i-4, "*w", Z,i-4, " + WW", i,i-4, "*w", W,i-4, sep=""))
        } else if (Z != "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-4, "*w", X,i-4, " + YX", i,i-4, "*w", Y,i-4, " + ZX", i,i-4, "*w", Z,i-4, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-4, "*w", X,i-4, " + YY", i,i-4, "*w", Y,i-4, " + ZY", i,i-4, "*w", Z,i-4, sep=""))
          cat("\n", paste("  w", Z,i, " ~ XZ", i,i-4, "*w", X,i-4, " + YZ", i,i-4, "*w", Y,i-4, " + ZZ", i,i-4, "*w", Z,i-4, sep=""))
        } # end (if Z/W)
      } # end (for i)
    } # end (if lag == 4)

    cat(rep("\n",2), "  '")

    # -- Run Model STARTSMLR -- #
    cat(rep("\n",2), "# -- Run Model STARTSMLR -- #")
    cat(rep("\n",2), "  STARTSMLR.fit <- suppressWarnings(lavaan::sem(STARTS,")
    cat("\n", "   ", data.source, ",")
    cat("\n", "   missing = 'fiml',")
    cat("\n", "   meanstructure = TRUE,")
    cat("\n", "   information = 'observed',")
    cat("\n", "   estimator = 'MLR'")
    cat("\n", "))")

  sink()  # Stop writing to file STARTS.txt
  ## ------------------------------- ##


  ## -- Execute STARTS.txt and request summary outputs-- ##
  source('STARTS.txt')
  print(lavaan::summary(STARTSMLR.fit, fit.measure = TRUE, standardized = TRUE, rsq = TRUE))
  if (lavaan::lavInspect(STARTSMLR.fit, what ="post.check") == FALSE) stop("The lavaan solution is non-admissible.")
  ## ------------------------------- ##



  ## -- Monte Carlo Simulation -- ##

  parEst <- lavaan::parameterEstimates(STARTSMLR.fit, remove.nonfree = TRUE)
  pest2 <- parEst[,5]  # Estimated Parameters
  pest3 <- lavaan::lavTech(STARTSMLR.fit, what = "vcov", add.labels = TRUE)  # Estimated Variance-Covariance of Estimated Parameters

    Invariance(parEst, pest2, pest3, no.path, MIset, no.compare, no.waves, lag, varI.eq, p, X, Y, Z, W)

    cat(rep("\n", 2))

 
  ## ----- Start writing script to STARTS.txt ----- ##
  sink('STARTS.txt')
    cat("\n", "# Specify the model (STARTS)", "\n")
    cat("\n", "STARTS <- '")

    cat(rep("\n",2), "  # -- Create between components (random intercepts) -- #")
    BX <- paste("  RI", X, " =~ 1*", X, "1", sep="")
    BY <- paste("  RI", Y, " =~ 1*", Y, "1", sep="")
    if (Z != "NULL") {
      BY <- paste("  RI", Z, " =~ 1*", Z, "1", sep="")
    } # end (if Z != "NULL")
    if (W != "NULL") {
      BY <- paste("  RI", W, " =~ 1*", W, "1", sep="")
    } # end (if W != "NULL")

    for (i in 2:no.waves) {
      BX <- paste(BX, " +1*", X, i, sep="")
      BY <- paste(BY, " +1*", Y, i, sep="")
      if (Z != "NULL") {
        BZ <- paste(BZ, " +1*", Z, i, sep="")
      } # end (if Z != "NULL")
      if (W != "NULL") {
        BW <- paste(BW, " +1*", W, i, sep="")
      } # end (if W != "NULL")
    } # end (for i)

    cat("\n", BX)
    cat("\n", BY)
    if (Z != "NULL") {
      cat("\n", BZ)
    } # end (if Z != "NULL")
    if (W != "NULL") {
      cat("\n", BW)
    } # end (if W != "NULL")

    # -- Create Residual Variance of Observed Variables -- #
    if (isFALSE(varI.eq)) {
      cat(rep("\n",2), "  # -- Create residual variance of observed variables -- #")
      cat("\n", "  ##########################################################")
      cat("\n", "  # Remove the subscripts for invariant indicator variance #")
      cat("\n", "  ##########################################################")
      for (i in 1:no.waves) {
        cat("\n", paste("  ", X, i, " ~~ varI", X, i, "*", X, i, sep=""))
        cat("\n", paste("  ", Y, i, " ~~ varI", Y, i, "*", Y, i, sep=""))
        if (Z != "NULL") {
          cat("\n", paste("  ", Z, i, " ~~ varI", Z, i, "*", Z, i, sep=""))
        } # end (if Z)
        if (W != "NULL") {
          cat("\n", paste("  ", W, i, " ~~ varI", W, i, "*", W, i, sep=""))
        } # end (if W)
      } # end (for i)###   ###
    } else {
      cat(rep("\n",2), "  # -- Create residual variance of observed variables -- #")
      cat("\n", "  ##########################################################")
      cat("\n", "  # Remove the subscripts for invariant indicator variance #")
      cat("\n", "  ##########################################################")
      for (i in 1:no.waves) {
        cat("\n", paste("  ", X, i, " ~~ varI", X, "*", X, i, sep=""))
        cat("\n", paste("  ", Y, i, " ~~ varI", Y, "*", Y, i, sep=""))
        if (Z != "NULL") {
          cat("\n", paste("  ", Z, i, " ~~ varI", Z, "*", Z, i, sep=""))
        } # end (if Z)
        if (W != "NULL") {
          cat("\n", paste("  ", W, i, " ~~ varI", W, "*", W, i, sep=""))
        } # end (if W)
      } # end (for i)###   ###
    } # end (if varI.eq == FALSE)

    # -- Create Latent Variables from Observed Variables -- #
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

    cat(rep("\n",2), "  # -- Estimate lagged effects between latent variables (Lag = 1 wave) -- #")
    cat("\n", "  #############################################")
    cat("\n", "  # Remove the subscripts for invariant paths #")
    cat("\n", "  #############################################")
    for (i in 2:no.waves) {
      if (Z == "NULL") {
        cat("\n", paste("  w", X,i, " ~ XX", i,i-1, "*w", X,i-1, " + YX", i,i-1, "*w", Y,i-1, sep=""))
        cat("\n", paste("  w", Y,i, " ~ XY", i,i-1, "*w", X,i-1, " + YY", i,i-1, "*w", Y,i-1, sep=""))
      } else if (W != "NULL") {
        cat("\n", paste("  w", X,i, " ~ XX", i,i-1, "*w", X,i-1, " + YX", i,i-1, "*w", Y,i-1, " + ZX", i,i-1, "*w", Z,i-1, " + WX", i,i-1, "*w", W,i-1, sep=""))
        cat("\n", paste("  w", Y,i, " ~ XY", i,i-1, "*w", X,i-1, " + YY", i,i-1, "*w", Y,i-1, " + ZY", i,i-1, "*w", Z,i-1, " + WY", i,i-1, "*w", W,i-1, sep=""))
        cat("\n", paste("  w", Z,i, " ~ XZ", i,i-1, "*w", X,i-1, " + YZ", i,i-1, "*w", Y,i-1, " + ZZ", i,i-1, "*w", Z,i-1, " + WZ", i,i-1, "*w", W,i-1, sep=""))
        cat("\n", paste("  w", W,i, " ~ XW", i,i-1, "*w", X,i-1, " + YW", i,i-1, "*w", Y,i-1, " + ZW", i,i-1, "*w", Z,i-1, " + WW", i,i-1, "*w", W,i-1, sep=""))
      } else if (Z != "NULL") {
        cat("\n", paste("  w", X,i, " ~ XX", i,i-1, "*w", X,i-1, " + YX", i,i-1, "*w", Y,i-1, " + ZX", i,i-1, "*w", Z,i-1, sep=""))
        cat("\n", paste("  w", Y,i, " ~ XY", i,i-1, "*w", X,i-1, " + YY", i,i-1, "*w", Y,i-1, " + ZY", i,i-1, "*w", Z,i-1, sep=""))
        cat("\n", paste("  w", Z,i, " ~ XZ", i,i-1, "*w", X,i-1, " + YZ", i,i-1, "*w", Y,i-1, " + ZZ", i,i-1, "*w", Z,i-1, sep=""))
      } # end (if Z)
    } # end (for i)

    if (lag == 2) {
      cat(rep("\n",2), "  # -- Estimate lagged effects between latent variables (Lag = 2 waves) -- #")
      cat("\n", "  #############################################")
      cat("\n", "  # Remove the subscripts for invariant paths #")
      cat("\n", "  #############################################")
      for (i in 3:no.waves) {
        if (Z == "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-2, "*w", X,i-2, " + YX", i,i-2, "*w", Y,i-2, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-2, "*w", X,i-2, " + YY", i,i-2, "*w", Y,i-2, sep=""))
        } else if (W != "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-2, "*w", X,i-2, " + YX", i,i-2, "*w", Y,i-2, " + ZX", i,i-2, "*w", Z,i-2, " + WX", i,i-2, "*w", W,i-2, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-2, "*w", X,i-2, " + YY", i,i-2, "*w", Y,i-2, " + ZY", i,i-2, "*w", Z,i-2, " + WY", i,i-2, "*w", W,i-2, sep=""))
          cat("\n", paste("  w", Z,i, " ~ XZ", i,i-2, "*w", X,i-2, " + YZ", i,i-2, "*w", Y,i-2, " + ZZ", i,i-2, "*w", Z,i-2, " + WZ", i,i-2, "*w", W,i-2, sep=""))
          cat("\n", paste("  w", W,i, " ~ XW", i,i-2, "*w", X,i-2, " + YW", i,i-2, "*w", Y,i-2, " + ZW", i,i-2, "*w", Z,i-2, " + WW", i,i-2, "*w", W,i-2, sep=""))
        } else if (Z != "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-2, "*w", X,i-2, " + YX", i,i-2, "*w", Y,i-2, " + ZX", i,i-2, "*w", Z,i-2, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-2, "*w", X,i-2, " + YY", i,i-2, "*w", Y,i-2, " + ZY", i,i-2, "*w", Z,i-2, sep=""))
          cat("\n", paste("  w", Z,i, " ~ XZ", i,i-2, "*w", X,i-2, " + YZ", i,i-2, "*w", Y,i-2, " + ZZ", i,i-2, "*w", Z,i-2, sep=""))
        } # end (if Z/W)
      } # end (for i)
    } # end (if lag == 2)

    if (lag == 3) {
      cat(rep("\n",2), "  # -- Estimate lagged effects between latent variables (Lag = 3 waves) -- #")
      cat("\n", "  #############################################")
      cat("\n", "  # Remove the subscripts for invariant paths #")
      cat("\n", "  #############################################")
      for (i in 3:no.waves) {
        if (Z == "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-3, "*w", X,i-3, " + YX", i,i-3, "*w", Y,i-3, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-3, "*w", X,i-3, " + YY", i,i-3, "*w", Y,i-3, sep=""))
        } else if (W != "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-3, "*w", X,i-3, " + YX", i,i-3, "*w", Y,i-3, " + ZX", i,i-3, "*w", Z,i-3, " + WX", i,i-3, "*w", W,i-3, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-3, "*w", X,i-3, " + YY", i,i-3, "*w", Y,i-3, " + ZY", i,i-3, "*w", Z,i-3, " + WY", i,i-3, "*w", W,i-3, sep=""))
          cat("\n", paste("  w", Z,i, " ~ XZ", i,i-3, "*w", X,i-3, " + YZ", i,i-3, "*w", Y,i-3, " + ZZ", i,i-3, "*w", Z,i-3, " + WZ", i,i-3, "*w", W,i-3, sep=""))
          cat("\n", paste("  w", W,i, " ~ XW", i,i-3, "*w", X,i-3, " + YW", i,i-3, "*w", Y,i-3, " + ZW", i,i-3, "*w", Z,i-3, " + WW", i,i-3, "*w", W,i-3, sep=""))
        } else if (Z != "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-3, "*w", X,i-3, " + YX", i,i-3, "*w", Y,i-3, " + ZX", i,i-3, "*w", Z,i-3, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-3, "*w", X,i-3, " + YY", i,i-3, "*w", Y,i-3, " + ZY", i,i-3, "*w", Z,i-3, sep=""))
          cat("\n", paste("  w", Z,i, " ~ XZ", i,i-3, "*w", X,i-3, " + YZ", i,i-3, "*w", Y,i-3, " + ZZ", i,i-3, "*w", Z,i-3, sep=""))
        } # end (if Z/W)
      } # end (for i)
    } # end (lag == 3)

   if (lag == 4) {
      cat(rep("\n",2), "  # -- Estimate lagged effects between latent variables (Lag = 4 waves) -- #")
      cat("\n", "  #############################################")
      cat("\n", "  # Remove the subscripts for invariant paths #")
      cat("\n", "  #############################################")
      for (i in 3:no.waves) {
        if (Z == "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-4, "*w", X,i-4, " + YX", i,i-4, "*w", Y,i-4, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-4, "*w", X,i-4, " + YY", i,i-4, "*w", Y,i-4, sep=""))
        } else if (W != "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-4, "*w", X,i-4, " + YX", i,i-4, "*w", Y,i-4, " + ZX", i,i-4, "*w", Z,i-4, " + WX", i,i-4, "*w", W,i-4, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-4, "*w", X,i-4, " + YY", i,i-4, "*w", Y,i-4, " + ZY", i,i-4, "*w", Z,i-4, " + WY", i,i-4, "*w", W,i-4, sep=""))
          cat("\n", paste("  w", Z,i, " ~ XZ", i,i-4, "*w", X,i-4, " + YZ", i,i-4, "*w", Y,i-4, " + ZZ", i,i-4, "*w", Z,i-4, " + WZ", i,i-4, "*w", W,i-4, sep=""))
          cat("\n", paste("  w", W,i, " ~ XW", i,i-4, "*w", X,i-4, " + YW", i,i-4, "*w", Y,i-4, " + ZW", i,i-4, "*w", Z,i-4, " + WW", i,i-4, "*w", W,i-4, sep=""))
        } else if (Z != "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-4, "*w", X,i-4, " + YX", i,i-4, "*w", Y,i-4, " + ZX", i,i-4, "*w", Z,i-4, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-4, "*w", X,i-4, " + YY", i,i-4, "*w", Y,i-4, " + ZY", i,i-4, "*w", Z,i-4, sep=""))
          cat("\n", paste("  w", Z,i, " ~ XZ", i,i-4, "*w", X,i-4, " + YZ", i,i-4, "*w", Y,i-4, " + ZZ", i,i-4, "*w", Z,i-4, sep=""))
        } # end (if Z/W)
      } # end (for i)
    } # end (lag == 4)


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


    cat(rep("\n",2), "  # -- Estimate variance and covariance of random intercepts -- #")
    cat("\n", "   RI", X, " ~~ RI", X, sep="")
    cat("\n", "   RI", Y, " ~~ RI", Y, sep="")
    if (Z != "NULL") {
      cat("\n", "   RI", Z, " ~~ RI", Z, sep="")
    } # end (if Z != "NULL")
    if (W != "NULL") {
      cat("\n", "   RI", W, " ~~ RI", W, sep="")
    } # end (if W != "NULL")
    cat("\n", "   RI", X, " ~~ RI", Y, sep="")
    if (Z != "NULL") {
      cat("\n", "   RI", X, " ~~ RI", Z, sep="")
      cat("\n", "   RI", Y, " ~~ RI", Z, sep="")
    } # end (if Z != "NULL")
    if (W != "NULL") {
      cat("\n", "   RI", X, " ~~ RI", W, sep="")
      cat("\n", "   RI", Y, " ~~ RI", W, sep="")
      cat("\n", "   RI", Z, " ~~ RI", W, sep="")
    } # end (if W != "NULL")


    cat(rep("\n",2), "  # -- Estimate (residual) variance of latent variables -- #")
    for (i in 1:no.waves) {
      cat("\n", paste("  w", X, i, " ~~ ", "w", X, i, sep=""))
      cat("\n", paste("  w", Y, i, " ~~ ", "w", Y, i, sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  w", Z, i, " ~~ ", "w", Z, i, sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  w", W, i, " ~~ ", "w", W, i, sep=""))
      } # end (if W)
    } # end (for i)


    cat(rep("\n",2), "  # -- Estimate grand means (intercepts) of observed variables -- #")
    for (i in 1:no.waves) {
      cat("\n", paste("  ", X, i, " ~ M", X, i, "*1", sep=""))
      cat("\n", paste("  ", Y, i, " ~ M", Y, i, "*1", sep=""))
      if (Z != "NULL") {
        cat("\n", paste("  ", Z, i, " ~ M", Z, i, "*1", sep=""))
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", paste("  ", W, i, " ~ M", W, i, "*1", sep=""))
      } # end (if W)
    } # end (for i)


    # -- Constrain covariance among RIX RIY, RIZ, RIW and wx1, wy1, wz1, ww1 -- #
    cat(rep("\n",2), "  # -- Constrain covariance among RIx, RIy, RIz, RIw and wx1, wy1, wz1, ww1 -- #")
    cat("\n", "   RI", X, " ~~ 0*w", X, "1", sep="")
    cat("\n", "   RI", X, " ~~ 0*w", Y, "1", sep="")
    if (Z != "NULL") {
      cat("\n", "   RI", X, " ~~ 0*w", Z, "1", sep="")
    } # end (if Z)
    if (W != "NULL") {
      cat("\n", "   RI", X, " ~~ 0*w", W, "1", sep="")
    } # end (if W)
    cat("\n", "   RI", Y, " ~~ 0*w", X, "1", sep="")
    cat("\n", "   RI", Y, " ~~ 0*w", Y, "1", sep="")
    if (Z != "NULL") {
      cat("\n", "   RI", Y, " ~~ 0*w", Z, "1", sep="")
    } # end (if Z)
    if (W != "NULL") {
      cat("\n", "   RI", Y, " ~~ 0*w", W, "1", sep="")
    } # end (if W)
    if (Z != "NULL") {
      cat("\n", "   RI", Z, " ~~ 0*w", X, "1", sep="")
      cat("\n", "   RI", Z, " ~~ 0*w", Y, "1", sep="")
      cat("\n", "   RI", Z, " ~~ 0*w", Z, "1", sep="")
    } # end (if Z)
    if (W != "NULL") {
      cat("\n", "   RI", Z, " ~~ 0*w", W, "1", sep="")
      cat("\n", "   RI", W, " ~~ 0*w", X, "1", sep="")
      cat("\n", "   RI", W, " ~~ 0*w", Y, "1", sep="")
      cat("\n", "   RI", W, " ~~ 0*w", Z, "1", sep="")
      cat("\n", "   RI", W, " ~~ 0*w", W, "1", sep="")
    } # end (if W)


    cat(rep("\n",2), "  ##########################################")
    cat("\n", "  # Regression of observed variables on C1 #")
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

    # -- Run STARTSMLR -- #
    cat(rep("\n",2), "  STARTSMLR.fit <- lavaan::sem(STARTS, ")
    cat("\n", "   ", data.source, ",")
    cat("\n", "   missing = 'fiml',")
    cat("\n", "   meanstructure = TRUE,")
    cat("\n", "   information = 'observed',")
    cat("\n", "   estimator = 'MLR'")
    cat("\n", ")")

    # Request summary outputs
    cat(rep("\n",2), "  print(lavaan::summary(STARTSMLR.fit, fit.measure = TRUE, standardized = TRUE, rsq = TRUE))", "\n")

    cat(rep("\n",2))

  sink() # Stop writing to file

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
#' @param p critical p-value for pairwise comparisons (default is 0.001)
#' @param X name of variable X.
#' @param Y name of variable Y.
#' @param Z name of variable Z.
#' @param W name of variable W (Z must come before W).
#'
#' @return ALT outputs.
#' @export
#' @examples
#'
#' ## -- Example -- ##
#'
#' ALT(data.source="Data_A", 7, 2, X="EXPOSE", Y="INTENS")
#'

ALT <- function(data.source, no.waves, lag=1, p = 0.001, X, Y, Z="NULL", W = "NULL") {
 
  ## -- Check inputs -- ##
 
  if (no.waves < 3) stop("Minimum number of waves is 3")
 
  if (lag < 1) stop("Minimum number of lag is 1")
  if (lag > 4) stop("Maximum number of lag is 4")
 
  if ((no.waves - lag) < 1) stop("Number of waves must be greater than lag plus 1")
 
  if (p > 0.05) stop("p > 0.05 is not recommended")
  if (p < 0.0001) stop("p < 0.0001 is not recommended")
 
  if (Z == "NULL" & W != "NULL") stop("Z must be defined before W")
 
 
 
  ## ----- Creating Model ALT ----- ###
  sink('ALT.txt') # Start writing script to ALT.txt
 
    cat("\n", "## ----- Specify the model (ALT) ----- ##", "\n")
    cat("\n", "ALT <- '")
 
    # -- Constrain Residual Variance of Observed Variables to Zero -- #
    cat(rep("\n",2), "  # -- Constrain residual variance of observed variables to zero -- #")
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
 
    # -- Create Latent Variables from Observed Variables -- #
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

      # -- Create Between Components (Random Intercepts) from Latent Variables -- #
      cat(rep("\n",2), "  # -- Create between components (random intercepts) -- #")
      BX <- paste("  RI", X, " =~ 0*w", X, "1", sep="")
      BY <- paste("  RI", Y, " =~ 0*w", Y, "1", sep="")
      for (i in 2:no.waves) {
        BX <- paste(BX, " +1*w", X, i, sep="")
        BY <- paste(BY, " +1*w", Y, i, sep="")
      } # end (for i)
      cat("\n", BX)
      cat("\n", BY)
 
      if (Z != "NULL") {
        BZ <- paste("  RI", Z, " =~ 0*w", Z, "1", sep="")
        for (i in 2:no.waves) {
          BZ <- paste(BZ, " +1*w", Z, i, sep="")
        } # end (for i)
        cat("\n", BZ)
      } # end (if Z)
      if (W != "NULL") {
        BW <- paste("  RI", W, " =~ 0*w", W, "1", sep="")
        for (i in 2:no.waves) {
          BW <- paste(BW, " +1*w", W, i, sep="")
        } # end (for i)
        cat("\n", BW)
      } # end (if W)
 
      # -- Create Between Components (Random Slopes) from Latent Variables -- #
      cat(rep("\n",2), "  # -- Create between components (random slopes) -- #")
      BX <- paste("  RS", X, " =~ 0*w", X, "1", sep="")
      BY <- paste("  RS", Y, " =~ 0*w", Y, "1", sep="")
      for (i in 2:no.waves) {
        BX <- paste(BX, " +", (i-1), "*w", X, i, sep="")
        BY <- paste(BY, " +", (i-1), "*w", Y, i, sep="")
      } # end (for i)
      cat("\n", BX)
      cat("\n", BY)
 
      if (Z != "NULL") {
        BZ <- paste("  RS", Z, " =~ 0*w", Z, "1", sep="")
        for (i in 2:no.waves) {
          BZ <- paste(BZ, " +", (i-1), "*w", Z, i, sep="")
        } # end (for i)
        cat("\n", BZ)
      } # end (if Z)
      if (W != "NULL") {
        BW <- paste("  RS", W, " =~ 0*w", W, "1", sep="")
        for (i in 2:no.waves) {
          BW <- paste(BW, " +", (i-1), "*w", W, i, sep="")
        } # end (for i)
        cat("\n", BW)
      } # end (if W)
 
      # -- Estimate Covariances Between Residuals of Latent Variables -- #
      cat(rep("\n",2), "  # -- Estimate covariances between residuals of latent variables -- #")
      for (i in 2:no.waves) {
        cat("\n", paste("  w", X, i, " ~~ eXY", i, "*w", Y, i, sep=""))
        if (Z != "NULL") {
          cat("\n", paste("  w", X, i, " ~~ eXZ", i, "*w", Z, i, sep=""))
          cat("\n", paste("  w", Y, i, " ~~ eYZ", i, "*w", Z, i, sep=""))
        }  # end (if Z)
        if (W != "NULL") {
          cat("\n", paste("  w", X, i, " ~~ eXW", i, "*w", W, i, sep=""))
          cat("\n", paste("  w", Y, i, " ~~ eYW", i, "*w", W, i, sep=""))
          cat("\n", paste("  w", Z, i, " ~~ eZW", i, "*w", W, i, sep=""))
        }  # end (if W)
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
 
 
      # -- Estimate (Residual) Variance of Latent Variables -- #
      cat(rep("\n",2), "  # -- Estimate (residual) variance of latent variables -- #")
      for (i in 1:no.waves) {
        cat("\n", paste("  w", X, i, " ~~ ", "w", X, i, sep=""))
        cat("\n", paste("  w", Y, i, " ~~ ", "w", Y, i, sep=""))
        if (Z != "NULL") {
          cat("\n", paste("  w", Z, i, " ~~ ", "w", Z, i, sep=""))
        } # end (if Z)
        if (W != "NULL") {
          cat("\n", paste("  w", W, i, " ~~ ", "w", W, i, sep=""))
        } # end (if W)
      } # end (for i)
 
      # -- Estimate Grand Means (Intercepts) of Observed Variables -- #
      cat(rep("\n",2), "  # -- Estimate grand means (intercepts) of observed variables -- #")
      for (i in 1:no.waves) {
        cat("\n", paste("  ", X, i, " ~ M", X, i, "*1", sep=""))
        cat("\n", paste("  ", Y, i, " ~ M", Y, i, "*1", sep=""))
        if (Z != "NULL") {
          cat("\n", paste("  ", Z, i, " ~ M", Z, i, "*1", sep=""))
        } # end (if Z)
        if (W != "NULL") {
          cat("\n", paste("  ", W, i, " ~ M", W, i, "*1", sep=""))
        } # (if W)
      } # end (for i)
 
      # -- Estimate covariance among RIX RIY, RIZ, RIW and wx1, wy1, wz1, ww1 -- #
      cat(rep("\n",2), "  # -- Estimate covariance among RIx, RIy, RIz, RIw and wx1, wy1, wz1, ww1 -- #")
      cat("\n", "   RI", X, " ~~ w", X, "1", sep="")
      cat("\n", "   RI", X, " ~~ w", Y, "1", sep="")
      if (Z != "NULL") {
        cat("\n", "   RI", X, " ~~ w", Z, "1", sep="")
      } # end (if Z != "NULL")
      if (W != "NULL") {
        cat("\n", "   RI", X, " ~~ w", W, "1", sep="")
      } # end (if W != "NULL")
 
      cat("\n", "   RI", Y, " ~~ w", X, "1", sep="")
      cat("\n", "   RI", Y, " ~~ w", Y, "1", sep="")
      if (Z != "NULL") {
        cat("\n", "   RI", Y, " ~~ w", Z, "1", sep="")
      } # end (if Z != "NULL")
      if (W != "NULL") {
        cat("\n", "   RI", Y, " ~~ w", W, "1", sep="")
      } # end (if W != "NULL")
      if (Z != "NULL") {
        cat("\n", "   RI", Z, " ~~ w", X, "1", sep="")
        cat("\n", "   RI", Z, " ~~ w", Y, "1", sep="")
        cat("\n", "   RI", Z, " ~~ w", Z, "1", sep="")
      } # end (if Z != "NULL")
      if (W != "NULL") {
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
      if (Z != "NULL") {
        cat("\n", "   RS", X, " ~~ w", Z, "1", sep="")
      } # end (if Z != "NULL")
      if (W != "NULL") {
        cat("\n", "   RS", X, " ~~ w", W, "1", sep="")
      } # end (if W != "NULL")
 
      cat("\n", "   RS", Y, " ~~ w", X, "1", sep="")
      cat("\n", "   RS", Y, " ~~ w", Y, "1", sep="")
      if (Z != "NULL") {
        cat("\n", "   RS", Y, " ~~ w", Z, "1", sep="")
      } # end (if Z != "NULL")
      if (W != "NULL") {
        cat("\n", "   RS", Y, " ~~ w", W, "1", sep="")
      } # end (if W != "NULL")
      if (Z != "NULL") {
        cat("\n", "   RS", Z, " ~~ w", X, "1", sep="")
        cat("\n", "   RS", Z, " ~~ w", Y, "1", sep="")
        cat("\n", "   RS", Z, " ~~ w", Z, "1", sep="")
      } # end (if Z != "NULL")
      if (W != "NULL") {
        cat("\n", "   RS", Z, " ~~ w", W, "1", sep="")
        cat("\n", "   RS", W, " ~~ w", X, "1", sep="")
        cat("\n", "   RS", W, " ~~ w", Y, "1", sep="")
        cat("\n", "   RS", W, " ~~ w", Z, "1", sep="")
        cat("\n", "   RS", W, " ~~ w", W, "1", sep="")
      } # end (if W != "NULL")
 
 
      # -- Estimate Lagged Effects Between Latent Variables -- #
      cat(rep("\n",2), "  # -- Estimate lagged effects between latent variables (Lag = 1 wave) -- #")
      for (i in 2:no.waves) {
        if (Z == "NULL") {
          cat("\n", paste("  w", X, i, " ~ XX", i,i-1, "*w", X,i-1, " + YX", i,i-1, "*w", Y,i-1, sep=""))
          cat("\n", paste("  w", Y, i, " ~ XY", i,i-1, "*w", X,i-1, " + YY", i,i-1, "*w", Y,i-1, sep=""))
        } else if (W != "NULL") {
          cat("\n", paste("  w", X, i, " ~ XX", i,i-1, "*w", X,i-1, " + YX", i,i-1, "*w", Y,i-1, " + ZX", i,i-1, "*w", Z,i-1, " + WX", i,i-1, "*w", W,i-1, sep=""))
          cat("\n", paste("  w", Y, i, " ~ XY", i,i-1, "*w", X,i-1, " + YY", i,i-1, "*w", Y,i-1, " + ZY", i,i-1, "*w", Z,i-1, " + WY", i,i-1, "*w", W,i-1, sep=""))
          cat("\n", paste("  w", Z, i, " ~ XZ", i,i-1, "*w", X,i-1, " + YZ", i,i-1, "*w", Y,i-1, " + ZZ", i,i-1, "*w", Z,i-1, " + WZ", i,i-1, "*w", W,i-1, sep=""))
          cat("\n", paste("  w", W, i, " ~ XW", i,i-1, "*w", X,i-1, " + YW", i,i-1, "*w", Y,i-1, " + ZW", i,i-1, "*w", Z,i-1, " + WW", i,i-1, "*w", W,i-1, sep=""))
        } else if (Z != "NULL") {
          cat("\n", paste("  w", X, i, " ~ XX", i,i-1, "*w", X,i-1, " + YX", i,i-1, "*w", Y,i-1, " + ZX", i,i-1, "*w", Z,i-1, sep=""))
          cat("\n", paste("  w", Y, i, " ~ XY", i,i-1, "*w", X,i-1, " + YY", i,i-1, "*w", Y,i-1, " + ZY", i,i-1, "*w", Z,i-1, sep=""))
          cat("\n", paste("  w", Z, i, " ~ XZ", i,i-1, "*w", X,i-1, " + YZ", i,i-1, "*w", Y,i-1, " + ZZ", i,i-1, "*w", Z,i-1, sep=""))
        } # end (if Z/W)
      } # end (for i)
 
      if (lag == 2) {
        cat(rep("\n",2), "  # -- Estimate lagged effects between within-person centered variables (Lag = 2 waves) -- #")
        for (i in 3:no.waves) {
          if (Z == "NULL") {
            cat("\n", paste("  w", X,i, " ~ XX", i,i-2, "*w", X,i-2, " + YX", i,i-2, "*w", Y,i-2, sep=""))
            cat("\n", paste("  w", Y,i, " ~ XY", i,i-2, "*w", X,i-2, " + YY", i,i-2, "*w", Y,i-2, sep=""))
          } else if (W != "NULL") {
            cat("\n", paste("  w", X,i, " ~ XX", i,i-2, "*w", X,i-2, " + YX", i,i-2, "*w", Y,i-2, " + ZX", i,i-2, "*w", Z,i-2, " + WX", i,i-2, "*w", W,i-2, sep=""))
            cat("\n", paste("  w", Y,i, " ~ XY", i,i-2, "*w", X,i-2, " + YY", i,i-2, "*w", Y,i-2, " + ZY", i,i-2, "*w", Z,i-2, " + WY", i,i-2, "*w", W,i-2, sep=""))
            cat("\n", paste("  w", Z,i, " ~ XZ", i,i-2, "*w", X,i-2, " + YZ", i,i-2, "*w", Y,i-2, " + ZZ", i,i-2, "*w", Z,i-2, " + WZ", i,i-2, "*w", W,i-2, sep=""))
            cat("\n", paste("  w", W,i, " ~ XW", i,i-2, "*w", X,i-2, " + YW", i,i-2, "*w", Y,i-2, " + ZW", i,i-2, "*w", Z,i-2, " + WW", i,i-2, "*w", W,i-2, sep=""))
          } else if (Z != "NULL") {
            cat("\n", paste("  w", X,i, " ~ XX", i,i-2, "*w", X,i-2, " + YX", i,i-2, "*w", Y,i-2, " + ZX", i,i-2, "*w", Z,i-2, sep=""))
            cat("\n", paste("  w", Y,i, " ~ XY", i,i-2, "*w", X,i-2, " + YY", i,i-2, "*w", Y,i-2, " + ZY", i,i-2, "*w", Z,i-2, sep=""))
            cat("\n", paste("  w", Z,i, " ~ XZ", i,i-2, "*w", X,i-2, " + YZ", i,i-2, "*w", Y,i-2, " + ZZ", i,i-2, "*w", Z,i-2, sep=""))
          } # end (if Z/W)
        } # end (for i)
      } # end (if lag == 2)
 
      if (lag == 3) {
        cat(rep("\n",2), "  # -- Estimate lagged effects between within-person centered variables (Lag = 3 waves) -- #")
        for (i in 4:no.waves) {
          if (Z == "NULL") {
            cat("\n", paste("  w", X,i, " ~ XX", i,i-3, "*w", X,i-3, " + YX", i,i-3, "*w", Y,i-3, sep=""))
            cat("\n", paste("  w", Y,i, " ~ XY", i,i-3, "*w", X,i-3, " + YY", i,i-3, "*w", Y,i-3, sep=""))
          } else if (W != "NULL") {
            cat("\n", paste("  w", X,i, " ~ XX", i,i-3, "*w", X,i-3, " + YX", i,i-3, "*w", Y,i-3, " + ZX", i,i-3, "*w", Z,i-3, " + WX", i,i-3, "*w", W,i-3, sep=""))
            cat("\n", paste("  w", Y,i, " ~ XY", i,i-3, "*w", X,i-3, " + YY", i,i-3, "*w", Y,i-3, " + ZY", i,i-3, "*w", Z,i-3, " + WY", i,i-3, "*w", W,i-3, sep=""))
            cat("\n", paste("  w", Z,i, " ~ XZ", i,i-3, "*w", X,i-3, " + YZ", i,i-3, "*w", Y,i-3, " + ZZ", i,i-3, "*w", Z,i-3, " + WZ", i,i-3, "*w", W,i-3, sep=""))
            cat("\n", paste("  w", W,i, " ~ XW", i,i-3, "*w", X,i-3, " + YW", i,i-3, "*w", Y,i-3, " + ZW", i,i-3, "*w", Z,i-3, " + WW", i,i-3, "*w", W,i-3, sep=""))
          } else if (Z != "NULL") {
            cat("\n", paste("  w", X,i, " ~ XX", i,i-3, "*w", X,i-3, " + YX", i,i-3, "*w", Y,i-3, " + ZX", i,i-3, "*w", Z,i-3, sep=""))
            cat("\n", paste("  w", Y,i, " ~ XY", i,i-3, "*w", X,i-3, " + YY", i,i-3, "*w", Y,i-3, " + ZY", i,i-3, "*w", Z,i-3, sep=""))
            cat("\n", paste("  w", Z,i, " ~ XZ", i,i-3, "*w", X,i-3, " + YZ", i,i-3, "*w", Y,i-3, " + ZZ", i,i-3, "*w", Z,i-3, sep=""))
          } # end (if Z/W)
        } # end (for i)
      } # end (lag == 3)
 
      if (lag == 4) {
        cat(rep("\n",2), "  # -- Estimate lagged effects between within-person centered variables (Lag = 4 waves) -- #")
        for (i in 5:no.waves) {
          if (Z == "NULL") {
            cat("\n", paste("  w", X,i, " ~ XX", i,i-4, "*w", X,i-4, " + YX", i,i-4, "*w", Y,i-4, sep=""))
            cat("\n", paste("  w", Y,i, " ~ XY", i,i-4, "*w", X,i-4, " + YY", i,i-4, "*w", Y,i-4, sep=""))
          } else if (W != "NULL") {
            cat("\n", paste("  w", X,i, " ~ XX", i,i-4, "*w", X,i-4, " + YX", i,i-4, "*w", Y,i-4, " + ZX", i,i-4, "*w", Z,i-4, " + WX", i,i-4, "*w", W,i-4, sep=""))
            cat("\n", paste("  w", Y,i, " ~ XY", i,i-4, "*w", X,i-4, " + YY", i,i-4, "*w", Y,i-4, " + ZY", i,i-4, "*w", Z,i-4, " + WY", i,i-4, "*w", W,i-4, sep=""))
            cat("\n", paste("  w", Z,i, " ~ XZ", i,i-4, "*w", X,i-4, " + YZ", i,i-4, "*w", Y,i-4, " + ZZ", i,i-4, "*w", Z,i-4, " + WZ", i,i-4, "*w", W,i-4, sep=""))
            cat("\n", paste("  w", W,i, " ~ XW", i,i-4, "*w", X,i-4, " + YW", i,i-4, "*w", Y,i-4, " + ZW", i,i-4, "*w", Z,i-4, " + WW", i,i-4, "*w", W,i-4, sep=""))
          } else if (Z != "NULL") {
            cat("\n", paste("  w", X,i, " ~ XX", i,i-4, "*w", X,i-4, " + YX", i,i-4, "*w", Y,i-4, " + ZX", i,i-4, "*w", Z,i-4, sep=""))
            cat("\n", paste("  w", Y,i, " ~ XY", i,i-4, "*w", X,i-4, " + YY", i,i-4, "*w", Y,i-4, " + ZY", i,i-4, "*w", Z,i-4, sep=""))
            cat("\n", paste("  w", Z,i, " ~ XZ", i,i-4, "*w", X,i-4, " + YZ", i,i-4, "*w", Y,i-4, " + ZZ", i,i-4, "*w", Z,i-4, sep=""))
          } # end (if Z/W)
        } # end (for i)
      } # end (if lag == 4)
 
      cat(rep("\n",2), "  '")
 
      # -- Run Model ALTMLR -- #
      cat(rep("\n",2), "# -- Run Model ALTMLR -- #")
      cat(rep("\n",2), "  ALTMLR.fit <- suppressWarnings(lavaan::sem(ALT,")
      cat("\n", "   ", data.source, ",")
      cat("\n", "   missing = 'fiml',")
      cat("\n", "   meanstructure = TRUE,")
      cat("\n", "   information = 'observed',")
      cat("\n", "   estimator = 'MLR'")
      cat("\n", "))")
 
    sink()  # Stop writing to file ALT.txt
    ## ------------------------------- ##
 
 
    ## -- Execute ALT.txt and request summary outputs-- ##
    source('ALT.txt')
    print(lavaan::summary(ALTMLR.fit, fit.measure = TRUE, standardized = TRUE, rsq = TRUE))
    if (lavaan::lavInspect(ALTMLR.fit, what ="post.check") == FALSE) stop("The lavaan solution is non-admissible.")
    ## ------------------------------- ##
 
 
 
    ## -- Monte Carlo Simulation -- ##
 
    parEst <- lavaan::parameterEstimates(ALTMLR.fit, remove.nonfree = TRUE)
    pest2 <- parEst[,5]  # Estimated Parameters
    pest3 <- lavaan::lavTech(ALTMLR.fit, what = "vcov", add.labels = TRUE)  # Estimated Variance-Covariance of Estimated Parameters
 
    varI.eq = TRUE
    Invariance(parEst, pest2, pest3, no.path, MIset, no.compare, no.waves, lag, varI.eq, p, X, Y, Z, W)

    cat(rep("\n", 2))


    ## ----- Start writing script to ALT.txt ----- ##
    sink('ALT.txt')
      cat("\n", "# Specify the model (ALT)", "\n")
      cat("\n", "ALT <- '")
 
      cat(rep("\n",2), "  # -- Create between components (random intercepts) -- #")
      BX <- paste("  RI", X, " =~ 0*w", X, "1", sep="")
      BY <- paste("  RI", Y, " =~ 0*w", Y, "1", sep="")
      if (Z != "NULL") {
        BY <- paste("  RI", Z, " =~ 0*w", Z, "1", sep="")
      } # end (if Z != "NULL")
      if (W != "NULL") {
        BY <- paste("  RI", W, " =~ 0*w", W, "1", sep="")
      } # end (if W != "NULL")
 
      for (i in 2:no.waves) {
        BX <- paste(BX, " +1*w", X, i, sep="")
        BY <- paste(BY, " +1*w", Y, i, sep="")
        if (Z != "NULL") {
          BZ <- paste(BZ, " +1*w", Z, i, sep="")
        } # end (if Z != "NULL")
        if (W != "NULL") {
          BW <- paste(BW, " +1*w", W, i, sep="")
        } # end (if W != "NULL")
      } # end (for i)
 
      cat("\n", BX)
      cat("\n", BY)
      if (Z != "NULL") {
        cat("\n", BZ)
      } # end (if Z != "NULL")
      if (W != "NULL") {
        cat("\n", BW)
      } # end (if W != "NULL")
 
 
      cat(rep("\n",2), "  # -- Create between components (random slopes) -- #")
      BX <- paste("  RS", X, " =~ 0*w", X, "1", sep="")
      BY <- paste("  RS", Y, " =~ 0*w", Y, "1", sep="")
      if (Z != "NULL") {
        BY <- paste("  RS", Z, " =~ 0*w", Z, "1", sep="")
      } # end (if Z != "NULL")
      if (W != "NULL") {
        BY <- paste("  RS", W, " =~ 0*w", W, "1", sep="")
      } # end (if W != "NULL")
 
      for (i in 2:no.waves) {
        BX <- paste(BX, " +", (i-1), "*w", X, i, sep="")
        BY <- paste(BY, " +", (i-1), "*w", Y, i, sep="")
        if (Z != "NULL") {
          BZ <- paste(BZ, " +", (i-1), "*w", Z, i, sep="")
        } # end (if Z != "NULL")
        if (W != "NULL") {
          BW <- paste(BW, " +", (i-1), "*w", W, i, sep="")
        } # end (if W != "NULL")
      } # end (for i)
 
      cat("\n", BX)
      cat("\n", BY)
      if (Z != "NULL") {
        cat("\n", BZ)
      } # end (if Z != "NULL")
      if (W != "NULL") {
        cat("\n", BW)
      } # end (if W != "NULL")
 
      # -- Constrain Residual Variance of Observed Variables to Zero -- #
      cat(rep("\n",2), "  # -- Constrain residual variance of observed variables to zero -- #")
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
 
      # -- Create Latent Variables from Observed Variables -- #
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
 
      cat(rep("\n",2), "  # -- Estimate lagged effects between latent variables (Lag = 1 wave) -- #")
      cat("\n", "  #############################################")
      cat("\n", "  # Remove the subscripts for invariant paths #")
      cat("\n", "  #############################################")
      for (i in 2:no.waves) {
        if (Z == "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-1, "*w", X,i-1, " + YX", i,i-1, "*w", Y,i-1, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-1, "*w", X,i-1, " + YY", i,i-1, "*w", Y,i-1, sep=""))
        } else if (W != "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-1, "*w", X,i-1, " + YX", i,i-1, "*w", Y,i-1, " + ZX", i,i-1, "*w", Z,i-1, " + WX", i,i-1, "*w", W,i-1, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-1, "*w", X,i-1, " + YY", i,i-1, "*w", Y,i-1, " + ZY", i,i-1, "*w", Z,i-1, " + WY", i,i-1, "*w", W,i-1, sep=""))
          cat("\n", paste("  w", Z,i, " ~ XZ", i,i-1, "*w", X,i-1, " + YZ", i,i-1, "*w", Y,i-1, " + ZZ", i,i-1, "*w", Z,i-1, " + WZ", i,i-1, "*w", W,i-1, sep=""))
          cat("\n", paste("  w", W,i, " ~ XW", i,i-1, "*w", X,i-1, " + YW", i,i-1, "*w", Y,i-1, " + ZW", i,i-1, "*w", Z,i-1, " + WW", i,i-1, "*w", W,i-1, sep=""))
        } else if (Z != "NULL") {
          cat("\n", paste("  w", X,i, " ~ XX", i,i-1, "*w", X,i-1, " + YX", i,i-1, "*w", Y,i-1, " + ZX", i,i-1, "*w", Z,i-1, sep=""))
          cat("\n", paste("  w", Y,i, " ~ XY", i,i-1, "*w", X,i-1, " + YY", i,i-1, "*w", Y,i-1, " + ZY", i,i-1, "*w", Z,i-1, sep=""))
          cat("\n", paste("  w", Z,i, " ~ XZ", i,i-1, "*w", X,i-1, " + YZ", i,i-1, "*w", Y,i-1, " + ZZ", i,i-1, "*w", Z,i-1, sep=""))
        } # end (if Z)
      } # end (for i)
 
      if (lag == 2) {
        cat(rep("\n",2), "  # -- Estimate lagged effects between latent variables (Lag = 2 waves) -- #")
        cat("\n", "  #############################################")
        cat("\n", "  # Remove the subscripts for invariant paths #")
        cat("\n", "  #############################################")
        for (i in 3:no.waves) {
          if (Z == "NULL") {
            cat("\n", paste("  w", X,i, " ~ XX", i,i-2, "*w", X,i-2, " + YX", i,i-2, "*w", Y,i-2, sep=""))
            cat("\n", paste("  w", Y,i, " ~ XY", i,i-2, "*w", X,i-2, " + YY", i,i-2, "*w", Y,i-2, sep=""))
          } else if (W != "NULL") {
            cat("\n", paste("  w", X,i, " ~ XX", i,i-2, "*w", X,i-2, " + YX", i,i-2, "*w", Y,i-2, " + ZX", i,i-2, "*w", Z,i-2, " + WX", i,i-2, "*w", W,i-2, sep=""))
            cat("\n", paste("  w", Y,i, " ~ XY", i,i-2, "*w", X,i-2, " + YY", i,i-2, "*w", Y,i-2, " + ZY", i,i-2, "*w", Z,i-2, " + WY", i,i-2, "*w", W,i-2, sep=""))
            cat("\n", paste("  w", Z,i, " ~ XZ", i,i-2, "*w", X,i-2, " + YZ", i,i-2, "*w", Y,i-2, " + ZZ", i,i-2, "*w", Z,i-2, " + WZ", i,i-2, "*w", W,i-2, sep=""))
            cat("\n", paste("  w", W,i, " ~ XW", i,i-2, "*w", X,i-2, " + YW", i,i-2, "*w", Y,i-2, " + ZW", i,i-2, "*w", Z,i-2, " + WW", i,i-2, "*w", W,i-2, sep=""))
          } else if (Z != "NULL") {
            cat("\n", paste("  w", X,i, " ~ XX", i,i-2, "*w", X,i-2, " + YX", i,i-2, "*w", Y,i-2, " + ZX", i,i-2, "*w", Z,i-2, sep=""))
            cat("\n", paste("  w", Y,i, " ~ XY", i,i-2, "*w", X,i-2, " + YY", i,i-2, "*w", Y,i-2, " + ZY", i,i-2, "*w", Z,i-2, sep=""))
            cat("\n", paste("  w", Z,i, " ~ XZ", i,i-2, "*w", X,i-2, " + YZ", i,i-2, "*w", Y,i-2, " + ZZ", i,i-2, "*w", Z,i-2, sep=""))
          } # end (if Z/W)
        } # end (for i)
      } # end (if lag == 2)
 
      if (lag == 3) {
        cat(rep("\n",2), "  # -- Estimate lagged effects between latent variables (Lag = 3 waves) -- #")
        cat("\n", "  #############################################")
        cat("\n", "  # Remove the subscripts for invariant paths #")
        cat("\n", "  #############################################")
        for (i in 3:no.waves) {
          if (Z == "NULL") {
            cat("\n", paste("  w", X,i, " ~ XX", i,i-3, "*w", X,i-3, " + YX", i,i-3, "*w", Y,i-3, sep=""))
            cat("\n", paste("  w", Y,i, " ~ XY", i,i-3, "*w", X,i-3, " + YY", i,i-3, "*w", Y,i-3, sep=""))
          } else if (W != "NULL") {
            cat("\n", paste("  w", X,i, " ~ XX", i,i-3, "*w", X,i-3, " + YX", i,i-3, "*w", Y,i-3, " + ZX", i,i-3, "*w", Z,i-3, " + WX", i,i-3, "*w", W,i-3, sep=""))
            cat("\n", paste("  w", Y,i, " ~ XY", i,i-3, "*w", X,i-3, " + YY", i,i-3, "*w", Y,i-3, " + ZY", i,i-3, "*w", Z,i-3, " + WY", i,i-3, "*w", W,i-3, sep=""))
            cat("\n", paste("  w", Z,i, " ~ XZ", i,i-3, "*w", X,i-3, " + YZ", i,i-3, "*w", Y,i-3, " + ZZ", i,i-3, "*w", Z,i-3, " + WZ", i,i-3, "*w", W,i-3, sep=""))
            cat("\n", paste("  w", W,i, " ~ XW", i,i-3, "*w", X,i-3, " + YW", i,i-3, "*w", Y,i-3, " + ZW", i,i-3, "*w", Z,i-3, " + WW", i,i-3, "*w", W,i-3, sep=""))
          } else if (Z != "NULL") {
            cat("\n", paste("  w", X,i, " ~ XX", i,i-3, "*w", X,i-3, " + YX", i,i-3, "*w", Y,i-3, " + ZX", i,i-3, "*w", Z,i-3, sep=""))
            cat("\n", paste("  w", Y,i, " ~ XY", i,i-3, "*w", X,i-3, " + YY", i,i-3, "*w", Y,i-3, " + ZY", i,i-3, "*w", Z,i-3, sep=""))
            cat("\n", paste("  w", Z,i, " ~ XZ", i,i-3, "*w", X,i-3, " + YZ", i,i-3, "*w", Y,i-3, " + ZZ", i,i-3, "*w", Z,i-3, sep=""))
          } # end (if Z/W)
        } # end (for i)
      } # end (lag == 3)
 
     if (lag == 4) {
        cat(rep("\n",2), "  # -- Estimate lagged effects between latent variables (Lag = 4 waves) -- #")
        cat("\n", "  #############################################")
        cat("\n", "  # Remove the subscripts for invariant paths #")
        cat("\n", "  #############################################")
        for (i in 3:no.waves) {
          if (Z == "NULL") {
            cat("\n", paste("  w", X,i, " ~ XX", i,i-4, "*w", X,i-4, " + YX", i,i-4, "*w", Y,i-4, sep=""))
            cat("\n", paste("  w", Y,i, " ~ XY", i,i-4, "*w", X,i-4, " + YY", i,i-4, "*w", Y,i-4, sep=""))
          } else if (W != "NULL") {
            cat("\n", paste("  w", X,i, " ~ XX", i,i-4, "*w", X,i-4, " + YX", i,i-4, "*w", Y,i-4, " + ZX", i,i-4, "*w", Z,i-4, " + WX", i,i-4, "*w", W,i-4, sep=""))
            cat("\n", paste("  w", Y,i, " ~ XY", i,i-4, "*w", X,i-4, " + YY", i,i-4, "*w", Y,i-4, " + ZY", i,i-4, "*w", Z,i-4, " + WY", i,i-4, "*w", W,i-4, sep=""))
            cat("\n", paste("  w", Z,i, " ~ XZ", i,i-4, "*w", X,i-4, " + YZ", i,i-4, "*w", Y,i-4, " + ZZ", i,i-4, "*w", Z,i-4, " + WZ", i,i-4, "*w", W,i-4, sep=""))
            cat("\n", paste("  w", W,i, " ~ XW", i,i-4, "*w", X,i-4, " + YW", i,i-4, "*w", Y,i-4, " + ZW", i,i-4, "*w", Z,i-4, " + WW", i,i-4, "*w", W,i-4, sep=""))
          } else if (Z != "NULL") {
            cat("\n", paste("  w", X,i, " ~ XX", i,i-4, "*w", X,i-4, " + YX", i,i-4, "*w", Y,i-4, " + ZX", i,i-4, "*w", Z,i-4, sep=""))
            cat("\n", paste("  w", Y,i, " ~ XY", i,i-4, "*w", X,i-4, " + YY", i,i-4, "*w", Y,i-4, " + ZY", i,i-4, "*w", Z,i-4, sep=""))
            cat("\n", paste("  w", Z,i, " ~ XZ", i,i-4, "*w", X,i-4, " + YZ", i,i-4, "*w", Y,i-4, " + ZZ", i,i-4, "*w", Z,i-4, sep=""))
          } # end (if Z/W)
        } # end (for i)
      } # end (lag == 4)

 
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
 
 
      cat(rep("\n",2), "  # -- Estimate (residual) variance of latent variables -- #")
      for (i in 1:no.waves) {
        cat("\n", paste("  w", X, i, " ~~ ", "w", X, i, sep=""))
        cat("\n", paste("  w", Y, i, " ~~ ", "w", Y, i, sep=""))
        if (Z != "NULL") {
          cat("\n", paste("  w", Z, i, " ~~ ", "w", Z, i, sep=""))
        } # end (if Z)
        if (W != "NULL") {
          cat("\n", paste("  w", W, i, " ~~ ", "w", W, i, sep=""))
        } # end (if W)
      } # end (for i)

 
      cat(rep("\n",2), "  # -- Estimate grand means (intercepts) of observed variables -- #")
      for (i in 1:no.waves) {
        cat("\n", paste("  ", X, i, " ~ M", X, i, "*1", sep=""))
        cat("\n", paste("  ", Y, i, " ~ M", Y, i, "*1", sep=""))
        if (Z != "NULL") {
          cat("\n", paste("  ", Z, i, " ~ M", Z, i, "*1", sep=""))
        } # end (if Z)
        if (W != "NULL") {
          cat("\n", paste("  ", W, i, " ~ M", W, i, "*1", sep=""))
        } # end (if W)
      } # end (for i)
 
 
      # -- Estimate covariance among RIX RIY, RIZ, RIW and wx1, wy1, wz1, ww1 -- #
      cat(rep("\n",2), "  # -- Estimate covariance among RIx, RIy, RIz, RIw and wx1, wy1, wz1, ww1 -- #")
      cat("\n", "   RI", X, " ~~ w", X, "1", sep="")
      cat("\n", "   RI", X, " ~~ w", Y, "1", sep="")
      if (Z != "NULL") {
        cat("\n", "   RI", X, " ~~ w", Z, "1", sep="")
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", "   RI", X, " ~~ w", W, "1", sep="")
      } # end (if W)
      cat("\n", "   RI", Y, " ~~ w", X, "1", sep="")
      cat("\n", "   RI", Y, " ~~ w", Y, "1", sep="")
      if (Z != "NULL") {
        cat("\n", "   RI", Y, " ~~ w", Z, "1", sep="")
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", "   RI", Y, " ~~ w", W, "1", sep="")
      } # end (if W)
      if (Z != "NULL") {
        cat("\n", "   RI", Z, " ~~ w", X, "1", sep="")
        cat("\n", "   RI", Z, " ~~ w", Y, "1", sep="")
        cat("\n", "   RI", Z, " ~~ w", Z, "1", sep="")
      } # end (if Z)
      if (W != "NULL") {
        cat("\n", "   RI", Z, " ~~ w", W, "1", sep="")
        cat("\n", "   RI", W, " ~~ w", X, "1", sep="")
        cat("\n", "   RI", W, " ~~ w", Y, "1", sep="")
        cat("\n", "   RI", W, " ~~ w", Z, "1", sep="")
        cat("\n", "   RI", W, " ~~ w", W, "1", sep="")
      } # end (if W)
 
 
      # -- Estimate covariance among RSX RSY, RSZ, RSW and wx1, wy1, wz1, ww1 -- #
      cat(rep("\n",2), "  # -- Estimate covariance among RSx, RSy, RSz, RSw and wx1, wy1, wz1, ww1 -- #")
      cat("\n", "   RS", X, " ~~ w", X, "1", sep="")
      cat("\n", "   RS", X, " ~~ w", Y, "1", sep="")
      if (Z != "NULL") {
        cat("\n", "   RS", X, " ~~ w", Z, "1", sep="")
      } # end (if Z != "NULL")
      if (W != "NULL") {
        cat("\n", "   RS", X, " ~~ w", W, "1", sep="")
      } # end (if W != "NULL")
 
      cat("\n", "   RS", Y, " ~~ w", X, "1", sep="")
      cat("\n", "   RS", Y, " ~~ w", Y, "1", sep="")
      if (Z != "NULL") {
        cat("\n", "   RS", Y, " ~~ w", Z, "1", sep="")
      } # end (if Z != "NULL")
      if (W != "NULL") {
        cat("\n", "   RS", Y, " ~~ w", W, "1", sep="")
      } # end (if W != "NULL")
      if (Z != "NULL") {
        cat("\n", "   RS", Z, " ~~ w", X, "1", sep="")
        cat("\n", "   RS", Z, " ~~ w", Y, "1", sep="")
        cat("\n", "   RS", Z, " ~~ w", Z, "1", sep="")
      } # end (if Z != "NULL")
      if (W != "NULL") {
        cat("\n", "   RS", Z, " ~~ w", W, "1", sep="")
        cat("\n", "   RS", W, " ~~ w", X, "1", sep="")
        cat("\n", "   RS", W, " ~~ w", Y, "1", sep="")
        cat("\n", "   RS", W, " ~~ w", Z, "1", sep="")
        cat("\n", "   RS", W, " ~~ w", W, "1", sep="")
      } # end (if W != "NULL")
 
 
      cat(rep("\n",2), "  ##########################################")
      cat("\n", "  # Regression of observed variables on C1 #")
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
      cat(rep("\n",2), "  print(lavaan::summary(ALTMLR.fit, fit.measure = TRUE, standardized = TRUE, rsq = TRUE))", "\n")
 
      cat(rep("\n",2))
 
    sink() # Stop writing to file
 
  }  # end (Function ALT)

## ========================================================================================== ##

