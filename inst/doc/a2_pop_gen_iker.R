## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = F,
  error = F
)

## ----setup--------------------------------------------------------------------
library(PhaseTypeR)

## -----------------------------------------------------------------------------
# for any integer n > 2
subint_mat_t_mrca <- function(n) {
  subint_mat <- diag(-choose(n:2, 2))
  subint_mat[-(n-1), -1] <- subint_mat[-(n-1), -1] - subint_mat[-(n-1), -(n-1)]
  return(subint_mat)
}

## -----------------------------------------------------------------------------

n <- 4
subint_mat <- subint_mat_t_mrca(n)
init_probs <- matrix(c(1, rep(0,n-2)), 1, n-1)


## -----------------------------------------------------------------------------

t_mrca_4 <- PH(subint_mat, init_probs)

t_mrca_4


## -----------------------------------------------------------------------------

mean(t_mrca_4)


## -----------------------------------------------------------------------------

var(t_mrca_4)


## -----------------------------------------------------------------------------

subint_mat_t_total <- function(n) {
  subint_mat = matrix(c(0), nrow = n-1, ncol = n-1)
  for (i in 1:n-1) {
    subint_mat[i, i] = - 0.5 * (n - i)
    if (i < n-1) {
      subint_mat[i, i+1] = -subint_mat[i, i]
    }
  }
  subint_mat
}


## -----------------------------------------------------------------------------

n <- 4
subint_mat <- subint_mat_t_total(n)
init_probs <- matrix(c(1, rep(0,n-2)), 1, n-1)


## -----------------------------------------------------------------------------

t_total_4 <- PH(subint_mat, init_probs)

t_total_4


## -----------------------------------------------------------------------------

mean(t_total_4)
var(t_total_4)


## -----------------------------------------------------------------------------
t_total_4_bis <- reward_phase_type(t_mrca_4, c(4, 3, 2))
t_total_4_bis

## -----------------------------------------------------------------------------
dPH(0.5, t_mrca_4)
dPH(1:3, t_total_4)

## ----fig.width=6, fig.height=4, fig.align='center'----------------------------

x <- seq(0, 10, 0.1)
y <- dPH(x, t_mrca_4)
plot(x, y, type = 'l', col = 'orange')
y2 <- dPH(x, t_total_4)
lines(x, y2, col = 'blue')
legend(6, 0.5, legend=c(expression('T'[MRCA]), expression('T'[total])),
       col=c("orange", "blue"), lty=1)
title('Density function (n=4)')


## -----------------------------------------------------------------------------
qPH(0.5, t_mrca_4)
qPH(c(0.25, 0.75), t_total_4)

## ----fig.width=6, fig.height=4, fig.align='center'----------------------------

x <- seq(0,0.99,0.01)
y <- qPH(x, t_total_4)
plot(x, y, type = 'l', col = 'blue')
y2 <- qPH(x, t_mrca_4)
lines(x, y2, col = 'orange')
title('Quantile function (n=4)')
legend(0.1, 10, legend=c(expression('T'[MRCA]), expression('T'[total])),
       col=c("orange", "blue"), lty=1)


## -----------------------------------------------------------------------------
pPH(0.5, t_mrca_4)
pPH(c(0.25, 0.75), t_total_4)

## ----fig.width=6, fig.height=4, fig.align='center'----------------------------

x <- seq(0, 10, 0.1)
y <- pPH(x, t_mrca_4)
plot(x, y, type = 'l', col = 'orange')
y <- pPH(x, t_total_4)
lines(x, y, col = 'blue')
title('Probability function (n=4)')
legend(6, 0.4, legend=c(expression('T'[MRCA]), expression('T'[total])),
       col=c("orange", "blue"), lty=1)


## -----------------------------------------------------------------------------
set.seed(0)
rPH(3, t_mrca_4)
rPH(10, t_total_4)

## ----fig.width=6, fig.height=4, fig.align='center'----------------------------
set.seed(0)
x <- rPH(10000, t_total_4)
hist(x, main = '10,000 random draws (n=4)', breaks = seq(0, 30, 0.5), ylim = c(0, 3000), xlim=c(0, 20),
     col=rgb(0,0,1,0.5))
x <- rPH(10000, t_mrca_4)
hist(x, breaks = seq(0, 30, 0.5), add = T, col=rgb(1,0.5,0,0.5))
legend(15, 2000, legend=c(expression('T'[MRCA]), expression('T'[total])),
       col=c("orange", "blue"), lty=1)
box()

## -----------------------------------------------------------------------------
theta <- 3

subint_mat <- t_total_4$subint_mat
init_probs <- t_total_4$init_probs

disc_subint_mat <- solve(diag(nrow(subint_mat)) - 2/theta * subint_mat)

## -----------------------------------------------------------------------------

s_tot_4_3 <- DPH(disc_subint_mat, init_probs)
s_tot_4_3


## -----------------------------------------------------------------------------
mean(s_tot_4_3)
var(s_tot_4_3)

## -----------------------------------------------------------------------------
subint_itons_4 <- matrix(c(-6,  6,  0,  0,
                            0, -3,  2,  1,
                            0,  0, -1,  0,
                            0,  0,  0, -1), nrow = 4, byrow = T)

## -----------------------------------------------------------------------------
kingman_4 <- PH(subint_itons_4)
kingman_4

## -----------------------------------------------------------------------------
reward <- c(4, 3, 2, 2)
seg_sites_cont_4 <- reward_phase_type(kingman_4, reward)
seg_sites_cont_4

## -----------------------------------------------------------------------------
theta <- 3

subint_mat <- seg_sites_cont_4$subint_mat
init_probs <- seg_sites_cont_4$init_probs

disc_subint_mat <- solve(diag(nrow(subint_mat)) - 2/theta * subint_mat)

s_tot_4_3_bis <- DPH(disc_subint_mat, init_probs)
s_tot_4_3_bis

## -----------------------------------------------------------------------------
mean(s_tot_4_3_bis)
var(s_tot_4_3_bis)

## -----------------------------------------------------------------------------
reward <- c(4, 2, 1, 0)
xi1_cont_4 <- reward_phase_type(kingman_4, reward)
xi1_cont_4

## -----------------------------------------------------------------------------
theta <- 3

subint_mat <- xi1_cont_4$subint_mat
init_probs <- xi1_cont_4$init_probs

disc_subint_mat <- solve(diag(nrow(subint_mat)) - 2/theta * subint_mat)

xi1_4_3 <- DPH(disc_subint_mat, init_probs)
xi1_4_3

## -----------------------------------------------------------------------------
reward <- c(0, 1, 0, 2)
xi2_cont_4 <- reward_phase_type(kingman_4, reward)

theta <- 3
subint_mat <- xi2_cont_4$subint_mat
init_probs <- xi2_cont_4$init_probs

disc_subint_mat <- solve(diag(nrow(subint_mat)) - 2/theta * subint_mat)
xi2_4_3 <- DPH(disc_subint_mat, init_probs)
xi2_4_3

## -----------------------------------------------------------------------------
reward <- c(0, 0, 1, 0)
xi3_cont_4 <- reward_phase_type(kingman_4, reward)

theta <- 3
subint_mat <- xi3_cont_4$subint_mat
init_probs <- xi3_cont_4$init_probs

disc_subint_mat <- solve(diag(nrow(subint_mat)) - 2/theta * subint_mat)
xi3_4_3 <- DPH(disc_subint_mat, init_probs)
xi3_4_3

## -----------------------------------------------------------------------------

mean(xi1_4_3)-1; mean(xi2_4_3)-1; mean(xi3_4_3)-1
var(xi1_4_3);  var(xi2_4_3);  var(xi3_4_3)


## -----------------------------------------------------------------------------
dDPH(4, s_tot_4_3)
dDPH(1:3, xi1_4_3)

## ----fig.width=6, fig.height=4, fig.align='center'----------------------------

x <- 1:20
y <- dDPH(x, s_tot_4_3)
plot(x, y, type = 'l', col = 'green', ylim = c(0, 0.2))
y2 <- dDPH(x, xi1_4_3)
lines(x, y2, col = 'red')
legend(15, 0.15, legend=c(expression('S'[total]), expression('X'[1])),
       col=c("green", "red"), lty=1)
title('Density function (n=4)')


## -----------------------------------------------------------------------------
qDPH(0.5, s_tot_4_3)
qDPH(c(0.25, 0.75), xi1_4_3)

## ----fig.width=6, fig.height=4, fig.align='center'----------------------------

x <- seq(0,0.99,0.01)
y <- qDPH(x, s_tot_4_3)
plot(x, y, type = 'l', col = 'green')
y2 <- qDPH(x, xi1_4_3)
lines(x, y2, col = 'red')
title('Quantile function (n=4)')
legend(0.1, 18, legend=c(expression('S'[total]), expression('X'[1])),
       col=c("green", "red"), lty=1)


## -----------------------------------------------------------------------------
pDPH(0.5, s_tot_4_3)
pDPH(c(0.25, 0.75), xi1_4_3)

## ----fig.width=6, fig.height=4, fig.align='center'----------------------------

x <- 0:15
y <- pDPH(x, s_tot_4_3)
plot(x, y, type = 'l', col = 'green', ylim = c(0, 1))
y <- pDPH(x, xi1_4_3)
lines(x, y, col = 'red')
title('Probability function (n=4)')
legend(10, 0.4, legend=c(expression('S'[total]), expression('X'[1])),
       col=c("green", "red"), lty=1)


## -----------------------------------------------------------------------------
set.seed(0)
rDPH(3, s_tot_4_3)
rDPH(10, xi1_4_3)

## ----fig.width=6, fig.height=4, fig.align='center'----------------------------
set.seed(0)
x <- rDPH(10000, s_tot_4_3)-1
hist(x, main = '10,000 random draws (n=4)',
     breaks = seq(0, 40, 1), col=rgb(0,1,0,0.5),
     ylim = c(0, 3500), xlim = c(0, 25))
x <- rDPH(10000, xi1_4_3)-1
hist(x, breaks = seq(0, 40, 1), add = T, col=rgb(1,0,0,0.5))
legend(15, 2500, legend=c(expression('S'[total]), expression('X'[1])),
       col=c("green", "red"), lty=1)
box()

## -----------------------------------------------------------------------------

subint_itons_4 <- matrix(c(-6,  6,  0,  0,
                            0, -3,  2,  1,
                            0,  0, -1,  0,
                            0,  0,  0, -1), nrow = 4, byrow = T)


## -----------------------------------------------------------------------------

reward_mat_4 <- matrix(c(4, 0, 0,
                         2, 1, 0,
                         1, 0, 1,
                         0, 2, 0), nrow = 4, byrow = T)

## -----------------------------------------------------------------------------

Y_i_phase_type <- MPH(subint_mat = subint_itons_4, 
                      reward_mat = reward_mat_4)

Y_i_phase_type


## -----------------------------------------------------------------------------

mean(Y_i_phase_type)


## -----------------------------------------------------------------------------

mean(Y_i_phase_type, 1)
mean(Y_i_phase_type, 2)
mean(Y_i_phase_type, 3)



## -----------------------------------------------------------------------------

var(Y_i_phase_type)


