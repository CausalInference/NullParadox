#################################################################################
# Project: G-Null Paradox
# Date: February 7, 2021
#################################################################################

source('step0_helperfunctions.R')

# Binary: Least Flexible Application
start_time <- Sys.time()
temp <- getsimres_alliter(flexibility = 1, binary_A_ind = TRUE, all_reps = 1:reps)
end_time <- Sys.time()
time_bin_1 <- end_time-start_time
save(time_bin_1, file = paste0('./Results/time_bin_1_', res_date, '.RData'))


# Binary: Moderately Flexible Application
start_time <- Sys.time()
temp <- getsimres_alliter(flexibility = 2, binary_A_ind = TRUE, all_reps = 1:reps)
end_time <- Sys.time()
time_bin_2 <- end_time-start_time
save(time_bin_2, file = paste0('./Results/time_bin_2_', res_date, '.RData'))


# Binary: Benchmark
start_time <- Sys.time()
temp <- getsimres_alliter(flexibility = 10, binary_A_ind = TRUE, all_reps = 1:reps)
end_time <- Sys.time()
time_bin_check <- end_time-start_time
save(time_bin_check, file = paste0('./Results/time_bin_check_', res_date, '.RData'))
