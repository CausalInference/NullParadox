#################################################################################
# Project: G-Null Paradox
# Date: February 7, 2021
#################################################################################

source('step0_helperfunctions.R')

# Continuous: Least Flexible Application
start_time <- Sys.time()
temp <- getsimres_alliter(flexibility = 1, binary_A_ind = FALSE, all_reps = 1:reps)
end_time <- Sys.time()
time_cont_1 <- end_time-start_time
save(time_cont_1, file = paste0('./Results/time_cont_1', res_date, '.RData'))


# Continuous: Moderately Flexible Application
start_time <- Sys.time()
temp <- getsimres_alliter(flexibility = 2, binary_A_ind = FALSE, all_reps = 1:reps)
end_time <- Sys.time()
time_cont_2 <- end_time-start_time
save(time_cont_2, file = paste0('./Results/time_cont_2_', res_date, '.RData'))


# Continuous: Benchmark
start_time <- Sys.time()
temp <- getsimres_alliter(flexibility = 10, binary_A_ind = FALSE, all_reps = 1:reps)
end_time <- Sys.time()
time_cont_check <- end_time-start_time
save(time_cont_check, file = paste0('./Results/time_cont_check_', res_date, '.RData'))
