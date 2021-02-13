#################################################################################
# Project: G-Null Paradox
# Date: February 7, 2021
#################################################################################

#### Load packages
library(data.table)
library(gfoRmula)
library(doParallel)
library(foreach)
setDTthreads(threads = 1)

#### Functions

## Generates longitudinal data sets
datagen <- function(i, K, alpha, beta, sigma1, theta, sigma2, binary_A_ind){
  A_lb <- 0
  A_ub <- 200
  Y_lb <- 0
  Y_ub <- 1000
  n_lags <- 10
  
  # Preallocate space in vectors
  id_             <- rep(as.numeric(i), K+n_lags-1)
  t0_             <- -(n_lags-1):(K-1)
  U_              <- rep(NA, K+n_lags-1)
  L_              <- rep(NA, K+n_lags-1)
  A_              <- rep(NA, K+n_lags-1)
  Y_              <- rep(NA, K+n_lags-1)
  
  # Generate U
  U <- runif(1, 0, 1)
  U_[1:(K+n_lags-1)] <- U

  # Generate data at time -n_lags + 1, ..., -1
  L_[1:(n_lags-1)] <- rbinom(n_lags-1, 1, plogis(alpha[1] + alpha[3] * U))
  A_[1:(n_lags-1)] <- rep(0, n_lags-1)
  
  for (j in n_lags:(K+n_lags-1)){
    L <- rbinom(1, 1, plogis(alpha[1] + alpha[2] * A_[j-1] + alpha[3] * U + alpha[4] * A_[j-1] * U))
    if (binary_A_ind){
      A <- rbinom(1, 1, plogis(beta[1] + beta[2] * A_[j-1] + beta[3] * L + beta[4] * A_[j-1] * L))
    } else {
      A <- rnorm(1, beta[1] + beta[2] * A_[j-1] + beta[3] * L + beta[4] * A_[j-1] * L, sigma1)
      while (A < A_lb | A > A_ub){
        A <- rnorm(1, beta[1] + beta[2] * A_[j-1] + beta[3] * L + beta[4] * A_[j-1] * L, sigma1)
      }
    }
    L_[j] <- L
    A_[j] <- A
  }
  Y <- rnorm(1, theta[1]+theta[2]*U, sigma2)
  while (Y < Y_lb | Y > Y_ub){
    Y <- rnorm(1, theta[1]+theta[2]*U, sigma2)
  }
  Y_[K+n_lags-1] <- Y
  
  # Consolidate data in a single data frame
  temp_data <- data.table(id = id_[1:(K+n_lags-1)],
                          t0 = t0_[1:(K+n_lags-1)],
                          U = U_[1:(K+n_lags-1)],
                          L = L_[1:(K+n_lags-1)],
                          A = A_[1:(K+n_lags-1)],
                          Y = Y_[1:(K+n_lags-1)])
  return(temp_data)
}

## Helper function for applying the parametric g-formula to simulated data
apply_gformula <- function(df, K_ind, binary_A_ind, repnum, flexibility){
  covnames <- c('L', 'A')
  intvars <- list(c('A'), c('A'))
  basecovs <- NA
  
  if (binary_A_ind){
    covtypes <- c('binary', 'binary')
    interventions <- list(list(c(static,rep(0, K[K_ind]))),
                          list(c(static, rep(1, K[K_ind]))))
    covlink <- c('logit', 'logit')
  } else {
    covtypes <- c('binary', 'normal')
    interventions <- list(list(c(static,rep(50, K[K_ind]))),
                          list(c(static, rep(150, K[K_ind]))))
    covlink <- c('logit', 'identity')
  }
  
  if (flexibility == 1){
    covmodels = c(L ~ lag1_A + lag_cumavg1_L,
                  A ~ 1)
    if (K_ind == 1){
      ymodel <- Y ~ cumavg_L + A + lag1_A
      histories <- c(lagged, lagavg, cumavg)
      histvars <- list(c('A'), c('L'), c('L'))
    } else {
      ymodel <- Y ~ cumavg_L + A + lag1_A + lag_cumavg2_A
      histories <- c(lagged, lagavg, cumavg)
      histvars <- list(c('A'), c('A', 'L'), c('L'))
    }
  } else if (flexibility == 2){
    covmodels = c(L ~ lag1_A + lag1_L + lag2_L + lag_cumavg3_L,
                  A ~ 1)
    if (K_ind == 1){
      ymodel <- Y ~ L + lag1_L + lag2_L + lag_cumavg3_L + A + lag1_A
      histories <- c(lagged, lagavg)
      histvars <- list(c('A', 'L'), c('L'))
    } else {
      ymodel <- Y ~ L + lag1_L + lag2_L + lag_cumavg3_L + A + lag1_A + lag_cumavg2_A
      histories <- c(lagged, lagavg)
      histvars <- list(c('A', 'L'), c('A', 'L'))
    }
  } else if (flexibility == 3){
    covmodels = c(L ~ lag1_A + lag1_L + lag2_L +
                    lag3_L + lag4_L + lag5_L + lag6_L
                  + lag7_L + lag8_L + lag9_L + lag10_L,
                  A ~ 1)
    if (K_ind == 1){
      ymodel <- Y ~ L + lag1_L + lag2_L + lag3_L + lag4_L + lag5_L + lag6_L + 
        lag7_L + lag8_L + lag9_L + lag10_L + A + lag1_A
      histories <- c(lagged)
      histvars <- list(c('A', 'L'))
    } else {
      ymodel <- Y ~ L + lag1_L + lag2_L + lag3_L + lag4_L + lag5_L + lag6_L + 
        lag7_L + lag8_L + lag9_L + lag10_L + A + lag1_A + lag_cumavg2_A
      histories <- c(lagged, lagavg)
      histvars <- list(c('A', 'L'), c('A'))
    }
  } else if (flexibility == 10){
    basecovs <- c('U')
    covmodels <- c(L ~ lag1_A + U + U * lag1_A,
                   A ~ 1)
    if (K_ind == 1){
      ymodel <- Y ~ U + A + lag1_A 
      histories <- c(lagged)
      histvars <- list(c('A'))
    } else {
      ymodel <- Y ~ U + A + lag1_A + lag_cumavg2_A
      histories <- c(lagged, lagavg)
      histvars <- list(c('A'), c('A'))
    }
  }
  covparams <- list(covlink = covlink, covmodels = covmodels)
  res <- gformula_continuous_eof(obs_data = df, id = 'id', time_name = 't0',
                                 outcome_name = 'Y', covnames = covnames,
                                 covtypes = covtypes, covparams = covparams,
                                 ymodel = ymodel, intvars = intvars,
                                 interventions = interventions,
                                 histories = histories, histvars = histvars,
                                 seed = 1234, nsamples = nsamples,
                                 ref_int = 2, show_progress = FALSE, 
                                 basecovs = basecovs)
  temp <- res$result
  output <- cbind(temp[Interv. == 1, 'g-form mean'],
                  temp[Interv. == 1, 'Mean SE'],
                  temp[Interv. == 1, 'Mean lower 95% CI'],
                  temp[Interv. == 1, 'Mean upper 95% CI'],
                  temp[Interv. == 2, 'g-form mean'],
                  temp[Interv. == 2, 'Mean SE'],
                  temp[Interv. == 2, 'Mean lower 95% CI'],
                  temp[Interv. == 2, 'Mean upper 95% CI'],
                  temp[Interv. == 1, 'Mean difference'],
                  temp[Interv. == 1, 'MD lower 95% CI'],
                  temp[Interv. == 1, 'MD upper 95% CI'],
                  repnum, K[K_ind])
  return(output)
}

## Driver function for the simulations
getsimres_singleiter <- function(repnum, K, alpha, beta, sigma1, theta, sigma2, 
                                 binary_A_ind, n, nsamples, parallel, seeds, flexibility){
  
  # K = 1
  set.seed(seeds)
  df <- lapply(as.list(1:n), FUN=function(ind){
    datagen(i = ind, K = K[1],
            alpha = alpha, beta = beta, sigma1 = sigma1, theta = theta, sigma2 = sigma2, 
            binary_A_ind = binary_A_ind)
    
  })
  df <- rbindlist(df)
  res1 <- apply_gformula(df = df, K_ind = 1, binary_A_ind = binary_A_ind, 
                         repnum = repnum, flexibility = flexibility)
  
  # K = 5
  set.seed(seeds+2)
  df <- lapply(as.list(1:n), FUN=function(ind){
    datagen(i = ind, K = K[2],
            alpha = alpha, beta = beta, sigma1 = sigma1, theta = theta, sigma2 = sigma2, 
            binary_A_ind = binary_A_ind)
  })
  df <- rbindlist(df)
  res2 <- apply_gformula(df = df, K_ind = 2, binary_A_ind = binary_A_ind, 
                         repnum = repnum, flexibility = flexibility)
  
  # K = 10
  set.seed(seeds+3)
  df <- lapply(as.list(1:n), FUN=function(ind){
    datagen(i = ind, K = K[3],
            alpha = alpha, beta = beta, sigma1 = sigma1, theta = theta, sigma2 = sigma2, 
            binary_A_ind = binary_A_ind)
  })
  df <- rbindlist(df)
  res3 <- apply_gformula(df = df, K_ind = 3, binary_A_ind = binary_A_ind, 
                         repnum = repnum, flexibility = flexibility)
  
  output <- rbind(res1, res2, res3)
  colnames(output) <- c("Int 1 Est", "Int 1 SE", "Int 1 lb", "Int 1 ub",
                        "Int 2 Est", "Int 2 SE", "Int 2 lb", "Int 2 ub",
                        "Mean difference", "MD lower 95% CI", "MD upper 95% CI",
                        "repetition", "K")
  
  temp <- ifelse(binary_A_ind, 'bin', 'cont')
  filename <- paste0('./Results/A=', temp, '_appl=', flexibility, '_iter=', 
                     repnum, '_date=', res_date, '.RData')
  save(output, file = filename)
  
  return(NULL)
}


## Function for getting all simulation results

getsimres_alliter <- function(flexibility, binary_A_ind, all_reps, ncores = 20){
  registerDoParallel(cores=ncores)
  if (binary_A_ind){
    temp <- foreach(repnum = all_reps) %dopar% {
      getsimres_singleiter(repnum = repnum, K = K, n = n,
                           alpha = alpha_bin, beta = beta_bin, sigma1 = NA,
                           theta = theta, sigma2 = sigma2,
                           binary_A_ind = binary_A_ind,
                           seeds = seeds2[repnum],
                           nsamples = nsamples,
                           flexibility = flexibility)
    }
    return(NULL)
  } else {
    temp <- foreach(repnum = all_reps) %dopar% {
      getsimres_singleiter(repnum = repnum, K = K, n = n,
                           alpha = alpha_cont, beta = beta_cont, sigma1 = sigma1,
                           theta = theta, sigma2 = sigma2,
                           binary_A_ind = binary_A_ind,
                           seeds = seeds1[repnum],
                           nsamples = nsamples,
                           flexibility = flexibility)
    }
    return(NULL)
  }
}


## Create results folder

if (!dir.exists('./Results')){
  dir.create('./Results')
}


## Set coefficients

K <- c(2, 6, 11)
alpha_cont <- c(1, -0.015, 1, 0.015)
beta_cont <- c(80, 0.1, 30, -0.05); sigma1 <- 25
alpha_bin <- c(0, -2.5, 1, 2.5)
beta_bin <- c(-1.25, 1, 1, 1)
theta <- c(350, 300); sigma2 <- 50
n <- 10000; 
nsamples <- 250 
reps <- 250 
seeds1 <- 2019 * c(1:reps)
seeds2 <- 2020 * c(1:reps)
res_date <- 'feb7_2021'