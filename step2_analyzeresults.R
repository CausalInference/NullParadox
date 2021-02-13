#################################################################################
# Project: G-Null Paradox
# Date: February 7, 2020
#################################################################################

library('xtable')
library('RColorBrewer')
library('vioplot')
source('step0_helperfunctions.R')

## Loading simulation results
for (flexibility in c(1,2,3,10)){
  for (outcome_type in c('bin', 'cont')){
    dfname <- paste('res', outcome_type, flexibility, sep = '.')
    assign(dfname, data.table())
    for (repnum in 1:reps){
      load(paste0('./Results/A=', outcome_type, '_appl=', flexibility, '_iter=', 
                  repnum, '_date=', res_date, '.RData'))
      
      #Note: Need to multiply by -1 because set ref_int = 2
      output$`Mean difference` <- -1 * output$`Mean difference`
      temp <- -1 * output$`MD lower 95% CI`
      output$`MD lower 95% CI` <- -1 * output$`MD upper 95% CI`
      output$`MD upper 95% CI` <- temp
      assign(dfname, rbind(get(dfname), output))
    }
  }
}




## Get performance table
get_performance <- function(res, truth, K){
  res$cov.yn <- ifelse(res$`MD lower 95% CI` <= 0 & res$`MD upper 95% CI` >= 0, 1, 0)
  res$cov.yn.int1 <- ifelse(res$`Int 1 lb` <= truth & res$`Int 1 ub` >= truth, 1, 0)
  res$cov.yn.int2 <- ifelse(res$`Int 2 lb` <= truth & res$`Int 2 ub` >= truth, 1, 0)
  
  f.table <- int1.table <- int2.table <- matrix(nrow = length(K), ncol = 3)
  colnames(f.table) <- colnames(int1.table) <- colnames(int2.table) <- c('Bias', 'SE', 'Coverage')
  rownames(f.table) <- rownames(int1.table) <- rownames(int2.table) <- paste0('K + 1 = ', K)
  
  f.table[,'Bias'] <- res[, mean(`Mean difference`), by = K]$V1 
  f.table[,'SE'] <- res[, sd(`Mean difference`), by = K]$V1
  f.table[,'Coverage'] <- res[, mean(cov.yn), by = K]$V1
  
  int1.table[,'Bias'] <- res[, mean(`Int 1 Est`) - true.y, by = K]$V1
  int1.table[,'SE'] <- res[, sd(`Int 1 Est`), by = K]$V1
  int1.table[,'Coverage'] <- res[, mean(cov.yn.int1), by = K]$V1
  
  int2.table[,'Bias'] <- res[, mean(`Int 2 Est`) - true.y, by = K]$V1
  int2.table[,'SE'] <- res[, sd(`Int 2 Est`), by = K]$V1
  int2.table[,'Coverage'] <- res[, mean(cov.yn.int2), by = K]$V1
  
  return(list(difference = round(f.table, 4), int1 = round(int1.table, 4), int2 = round(int2.table, 4)))
}
true.y <- theta[1] + 0.5 * theta[2]
performance.cont.1 <- get_performance(res.cont.1, true.y, K)
performance.cont.2 <- get_performance(res.cont.2, true.y, K)
performance.cont.3 <- get_performance(res.cont.3, true.y, K)
performance.cont.10 <- get_performance(res.cont.10, true.y, K)
performance.bin.1 <- get_performance(res.bin.1, true.y, K)
performance.bin.2 <- get_performance(res.bin.2, true.y, K)
performance.bin.3 <- get_performance(res.bin.3, true.y, K)
performance.bin.10 <- get_performance(res.bin.10, true.y, K)

## Create tables
table1 <- rbind(performance.cont.1$difference, 
                performance.cont.2$difference, 
                performance.cont.3$difference, 
                performance.cont.10$difference)
table2 <- rbind(performance.bin.1$difference, 
                performance.bin.2$difference, 
                performance.bin.3$difference, 
                performance.bin.10$difference)
xtable(cbind(table1, table2))


table1a <- rbind(performance.bin.1$int1,  
                 performance.bin.2$int1, 
                 performance.bin.3$int1, 
                 performance.bin.10$int1, 
                 performance.bin.1$int2,
                 performance.bin.2$int2, 
                 performance.bin.3$int2, 
                 performance.bin.10$int2)
xtable(table1a, digits = c(0, 2, 2, 2))

table2a <- rbind(performance.cont.1$int1,  
                performance.cont.2$int1, 
                performance.cont.3$int1, 
                performance.cont.10$int1, 
                performance.cont.1$int2,
                performance.cont.2$int2, 
                performance.cont.3$int2, 
                performance.cont.10$int2)
xtable(table2a, digits = c(0, 2, 2, 2))






# Figures

temp.bin <- rbind(res.bin.1, res.bin.2, res.bin.3, res.bin.10)
temp.bin$Method <- rep(1:4, each = length(K) * reps)
temp.bin$K <- temp.bin$K - 1

temp.cont <- rbind(res.cont.1,res.cont.2, res.cont.3, res.cont.10)
temp.cont$Method <- rep(1:4, each = length(K) * reps)
temp.cont$K <- temp.cont$K - 1

n <- 12
VEC <- c(1:3, 5:7, 9:11, 13:15)
at2 <- c(2, 6, 10, 14)
cols <- rep(brewer.pal(3, name = 'Blues'), 4)
labels <- c('Least\nFlexible', 'Moderately\nFlexible', 'Most\nFlexible', 'Benchmark\n')
cex.lab <- 1.25
cex.mtext <- 1.5

# ATE Res
pdf('SimRes.pdf', width = 12, height = 6)
par(mfrow = (c(1,2)))

par(mar=c(6.1, 6.1, 3.1, 2.1))
vioplot(`Mean difference` ~ K + Method,
        data = temp.bin, ylim = c(-35,35),
        main = '', ylab = '', xlab = '', col = cols, at = VEC, 
        drawRect = T, xaxt = 'n')
axis(side = 1, at = at2, labels = rep('', 4))
axis(side = 1, at = at2, labels = labels, line = 0.85, lty = 0)
mtext('A)', side = 3, line = 1.2, adj = -0.28, cex = cex.mtext)
title(xlab = 'G-Formula Application', cex.lab = cex.lab, line = 4)
title(ylab = 'Estimated ACE', cex.lab = cex.lab)
abline(h=0, lty = 3)

par(mar=c(6.1, 6.1, 3.1, 2.1))
vioplot(`Mean difference` ~ K + Method,
        data = temp.cont, ylim = c(-85,85), 
        main = '', ylab = '', xlab = '', col = cols, at = VEC, 
        drawRect = T, xaxt = 'n', yaxt = 'n')
axis(side = 1, at = at2, labels = rep('', 4))
axis(side = 1, at = at2, labels = labels, line = 0.85, lty = 0)
axis(side = 2, at = seq(from = -75, to = 75, by = 25))
mtext('B)', side = 3, line = 1.2, adj = -0.28, cex = cex.mtext)
title(xlab = 'G-Formula Application', cex.lab = cex.lab, line = 4)
title(ylab = 'Estimated ACE', cex.lab = cex.lab)
abline(h=0, lty = 3)
dev.off()


# Summary
pdf('SimRes_Appendix_Bin.pdf', width = 12, height = 6)
par(mfrow = (c(1,2)))

par(mar=c(6.1, 6.1, 3.1, 2.1))
vioplot(`Int 1 Est` ~ K + Method,
        data = temp.bin, ylim = c(480, 520), 
        main = '', ylab = '', xlab = '', col = cols, at = VEC, 
        drawRect = T, xaxt = 'n')
axis(side = 1, at = at2, labels = rep('', 4))
axis(side = 1, at = at2, labels = labels, line = 0.85, lty = 0)
mtext('A)', side = 3, line = 1.2, adj = -0.28, cex = cex.mtext)
title(xlab = 'G-Formula Application', cex.lab = cex.lab, line = 4)
title(ylab = expression('Estimated Mean Under ' ~ bar(a)[K] ~ '=' ~ bar(0)), cex.lab = cex.lab)
abline(h=500, lty = 3)

par(mar=c(6.1, 6.1, 3.1, 2.1))
vioplot(`Int 2 Est` ~ K + Method,
        data = temp.bin, ylim = c(480, 520), 
        main = '', ylab = '', xlab = '', col = cols, at = VEC, 
        drawRect = T, xaxt = 'n')
axis(side = 1, at = at2, labels = rep('', 4))
axis(side = 1, at = at2, labels = labels, line = 0.85, lty = 0)
mtext('B)', side = 3, line = 1.2, adj = -0.28, cex = cex.mtext)
title(xlab = 'G-Formula Application', cex.lab = cex.lab, line = 4)
title(ylab = expression('Estimated Mean Under ' ~ bar(a)[K] ~ '=' ~ bar(1)), cex.lab = cex.lab)
abline(h=500, lty = 3)

dev.off()


pdf('SimRes_Appendix_Cont.pdf', width = 12, height = 6)
par(mfrow = (c(1,2)))

par(mar=c(6.1, 6.1, 3.1, 2.1))
vioplot(`Int 1 Est` ~ K + Method,
        data = temp.cont, ylim = c(450, 550), 
        main = '', ylab = '', xlab = '', col = cols, at = VEC, 
        drawRect = T, xaxt = 'n')
axis(side = 1, at = at2, labels = rep('', 4))
axis(side = 1, at = at2, labels = labels, line = 0.85, lty = 0)
mtext('A)', side = 3, line = 1.2, adj = -0.28, cex = cex.mtext)
title(xlab = 'G-Formula Application', cex.lab = cex.lab, line = 4)
title(ylab = expression('Estimated Mean Under ' ~ bar(a)[K] ~ '=' ~ bar(50)), cex.lab = cex.lab)
abline(h=500, lty = 3)

par(mar=c(6.1, 6.1, 3.1, 2.1))
vioplot(`Int 2 Est` ~ K + Method,
        data = temp.cont, ylim = c(450, 550), 
        main = '', ylab = '', xlab = '', col = cols, at = VEC, 
        drawRect = T, xaxt = 'n')
axis(side = 1, at = at2, labels = rep('', 4))
axis(side = 1, at = at2, labels = labels, line = 0.85, lty = 0)
mtext('B)', side = 3, line = 1.2, adj = -0.28, cex = cex.mtext)
title(xlab = 'G-Formula Application', cex.lab = cex.lab, line = 4)
title(ylab = expression('Estimated Mean Under ' ~ bar(a)[K] ~ '=' ~ bar(150)), cex.lab = cex.lab)
abline(h=500, lty = 3)

dev.off()

