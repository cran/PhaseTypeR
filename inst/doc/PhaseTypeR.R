## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align = "center"
)

## ---- setup, include = FALSE--------------------------------------------------
library(PhaseTypeR)

## ----PH()_example-------------------------------------------------------------
subintensity_matrix <- matrix(c(-1.5, 0, 0,
                              1.5, -1, 0,
                              0, 1, -0.5), ncol = 3)
initial_probabilities <- c(0.9, 0.1, 0)
ph <- PH(subintensity_matrix, initial_probabilities)

ph


## -----------------------------------------------------------------------------
summary(ph)

## -----------------------------------------------------------------------------
subintensity_matrix <- matrix(c(-1.5, 0, 0,
                              1.5, -1, 0,
                              0, 1, -0.5), ncol = 3)
PH(subintensity_matrix)

## ----cont_example-------------------------------------------------------------
subintensity_matrix <- matrix(c(-1.5, 0, 0,
                              1.5, -1, 0,
                              0, 1, -0.5), ncol = 3)
initial_probabilities <- c(0.9, 0.1, 0)

ph <- PH(subintensity_matrix, initial_probabilities)

ph


## -----------------------------------------------------------------------------

net <- phase_type_to_network(ph)
set.seed(8)
plot(net, edge.curved=.1, edge.label = E(net)$weight,
       edge.color = ifelse(as_data_frame(net, what="edges")$from == 'V0',
                         'purple', 'grey'))


## -----------------------------------------------------------------------------
oldpar <- par()[c('mfrow', 'mar')]

## ----fig.width=5, fig.height=3------------------------------------------------
par(mfrow=c(2,3), mar=c(2,2,2,2))
set.seed(6)
for (i in c(0, 1, 3, 5, 10, 20)) {
  net <- phase_type_to_network(ph, i)
  plot(net, edge.arrow.size=0.5, edge.curved=.1, edge.label = E(net)$weight, main=paste0('time = ', i),
       edge.color = ifelse(as_data_frame(net, what="edges")$from == 'V0',
                         'purple', 'grey'))
}

## -----------------------------------------------------------------------------
par(oldpar)

## ----mean_var_cont, collapse = TRUE-------------------------------------------
cat('\n', 'Mean: ',mean(ph),'\n')
cat(' Variance: ',var(ph),'\n \n')

## -----------------------------------------------------------------------------
summary(ph)

## ----plot_cont----------------------------------------------------------------

x <- seq(0,8,length.out = 100)
pdf <- dPH(x, ph)

{plot(x, pdf, xlab = "Time to absorption", ylab = "PDF", col = "blue", 
     type = 'l', lwd=2)
title('Probability density function (PDF)')}

## ----plot_cont_2--------------------------------------------------------------
x <- seq(0,8,length.out = 100)
cdf <- pPH(x, ph)

{plot(x, cdf, xlab = "Time to absorption", ylab = "CDF", col = "orange", 
     type = 'l', lwd=2)
title('Cumulative density function (CDF)')
}

## -----------------------------------------------------------------------------
x <- seq(0.05, 0.95, 0.01)
{plot(x, qPH(x, ph), col = "green", type="l", lwd=2, xlim=c(0,1),
     ylab = "Time to absorbtion", xlab = "Quantile")
title('Quantile function')}

## -----------------------------------------------------------------------------
cat('10 random samples: \n', rPH(5, ph), '\n', rPH(5, ph))

## -----------------------------------------------------------------------------
hist(rPH(500, ph), breaks = 10)

## -----------------------------------------------------------------------------

set.seed(6)
tab_FullPH <- rFullPH(ph)
x <- c(tab_FullPH$time, sum(tab_FullPH$time)/length(tab_FullPH$time))
x <- cumsum(x)
y <- c(tab_FullPH$state, 4)

plot(x,y,type="n", xlim=c(0,x[length(x)]))
segments(c(0,x[-length(x)]),y,x,y)
points(c(0,x[-length(x)]),y,pch=16)
points(x[-length(x)],y[-length(x)],pch=1)
points(x[length(x)],y[length(x)],pch='>')



## ----disc_example-------------------------------------------------------------
subintensity_matrix <- matrix(c(0, 0.2, 0.8,
                              0.5, 0.5, 0,
                              0, 0, 0.4), ncol = 3, byrow = T)
initial_probabilities <- c(0.7, 0.3, 0)
dph <- DPH(subintensity_matrix, initial_probabilities)

dph


## -----------------------------------------------------------------------------

net <- phase_type_to_network(dph)
set.seed(8)
plot(net, edge.curved=.1, edge.label = E(net)$weight)


## ----mean_var_disc------------------------------------------------------------
cat('\n', 'Mean: ',mean(dph),'\n')
cat(' Variance: ',var(dph),'\n \n')

## ----plot_disc----------------------------------------------------------------
x <- seq(0,10,by=1)
pdf <- dDPH(x, dph)
{plot(x, pdf, xlab = "Time to absorption", ylab = "PDF and CDF", col = "blue", 
     type = 'l', lwd=2)
title('Probability function (PDF)')}

## ----plot_disc_2--------------------------------------------------------------
x <- seq(0,10,by=1)
cdf <- pDPH(x, dph)
{plot(x, cdf, xlab = "Time to absorption", ylab = "CDF", col = "orange", 
     type = 'l', lwd=2)
title('Probability function (CDF)')}

## ----plot_disc_quant----------------------------------------------------------
x <- seq(0.05, 0.95, 0.01)
{plot(x, qDPH(x, dph), col = "green", 
      ylab = "Time to absorption", 
      xlab = "Quantile",xlim=c(0,1),ylim=c(0,10),type="l")
title('Quantile function')}

## -----------------------------------------------------------------------------
hist(rDPH(500, dph), breaks = 10)

## ----reward_disc--------------------------------------------------------------
rwd.ph <- reward_phase_type(ph, c(1, 2, 3))
print(rwd.ph)

## -----------------------------------------------------------------------------
x <- seq( 0 , 20 ,length.out = 100)
pdf <- dPH(x,ph)
rwd.pdf <- dPH(x,rwd.ph)

{plot(x, pdf, xlab = "Time to absorption", ylab = "pdf", col = "orange", type = 'l')
lines(x, rwd.pdf, col = "blue")
title('Distribution functions before and after reward transformation')
legend("topright", legend=c('Original Phase-type', 'Reward transformed Phase-type'),
       col=c("orange", "blue"), lty=1,bty="n")}

## -----------------------------------------------------------------------------
rwd.dph <- reward_phase_type(dph, c(1, 2, 3))
print(rwd.dph)

## -----------------------------------------------------------------------------
x <- seq(0,20,by=1)
pdf <- dDPH(x, dph)
rwd.pdf <- dDPH(x, rwd.dph)
{plot(x, pdf, xlab = "Time to absorption", ylab = "pdf", col = "orange", 
     type = 'l', lwd=2)
lines(x, rwd.pdf, col = "blue",lwd=2)
title('Distribution functions before and after reward transformation')
legend("topright", bty="n", legend=c('Original DPH', 'Reward transformed DPH'),
       col=c("orange", "blue"), lty=1, lwd=2)}

## -----------------------------------------------------------------------------
reward_phase_type(ph, c(0, 2, 3))

## -----------------------------------------------------------------------------
rwd.ph <- reward_phase_type(ph, c(2, 0, 0))
print(rwd.ph)

## -----------------------------------------------------------------------------
x <- seq( 0 , 20 ,length.out = 100)
cdf <- pPH(x,ph)
rwd.cdf <- pPH(x,rwd.ph)

{plot(x, cdf, xlab = "Time to absorption", ylab = "pdf", col = "orange", type = 'l',ylim=c(0,1))
lines(x, rwd.cdf, col = "blue")
title('Cumulative distributions before and after reward transformation')
legend("bottomright", legend=c('Original Phase-type', 'Reward transformed Phase-type'),
       col=c("orange", "blue"), lty=1,bty="n")}

## -----------------------------------------------------------------------------
zero_dph <- reward_phase_type(dph, c(0, 1, 1))
print(zero_dph)

## -----------------------------------------------------------------------------
reward_phase_type(zero_dph, c(2, 3))

## ---- echo = F----------------------------------------------------------------
reward_phase_type(dph, c(0, 2, 3))

## -----------------------------------------------------------------------------
Rmat <- matrix(c(0, 1, 1,  2,
                 2, 1, 5,  2,
                 0, 1, 10, 2), nrow = 3, ncol=4, byrow=TRUE)
mph <- MPH(ph$subint_mat, ph$init_probs, reward_mat = Rmat)
print(mph)

## -----------------------------------------------------------------------------
mean(mph)

## -----------------------------------------------------------------------------
var(mph)

## -----------------------------------------------------------------------------
Rmat <- matrix(c(0, 1, 1,  2,
                 2, 1, 5,  2,
                 0, 1, 10, 2), nrow = 3, ncol=4, byrow=TRUE)
mdph <- MDPH(dph$subint_mat, dph$init_probs, reward_mat = Rmat)
print(mdph)

## -----------------------------------------------------------------------------
mean(mdph)

## -----------------------------------------------------------------------------
var(mdph)

## -----------------------------------------------------------------------------
set.seed(42)
rFullPH(ph)
set.seed(42)
rFullDPH(dph)
set.seed(42)
rFullMPH(mph)
set.seed(42)
rFullMDPH(mdph)

