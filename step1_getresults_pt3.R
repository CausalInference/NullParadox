#################################################################################
# Project: G-Null Paradox
# Date: February 7, 2021
#################################################################################

source('step0_helperfunctions.R')

# Continuous: Most Flexible Application
start_time <- Sys.time()
temp <- getsimres_alliter(flexibility = 3, binary_A_ind = FALSE, all_reps = 1:reps)
end_time <- Sys.time()
time_cont_3 <- end_time-start_time
save(time_cont_3, file = paste0('./Results/time_cont_3_', res_date, '.RData'))


# Binary: Most Flexible Application
start_time <- Sys.time()
temp <- getsimres_alliter(flexibility = 3, binary_A_ind = TRUE, all_reps = 1:reps)
end_time <- Sys.time()
time_bin_3 <- end_time-start_time
save(time_bin_3, file = paste0('./Results/time_bin_3_', res_date, '.RData'))


