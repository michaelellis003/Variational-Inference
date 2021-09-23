source("univariate-normal.R")
library(ggplot2)
library(patchwork)

## simulate data
N <- 20
mu <- 100
sigmasq <- 225
X <- rnorm(N, mean = mu, sd = sqrt(sigmasq))

## fit using variational inference
output <- vi(X, max_iter = 10)

## data frame for plotting data
df <- dplyr::tibble(
    mu = seq(80, 120, length=1000),
    sigmasq = seq(90, 800, length=1000),
    h_mu = dnorm(mu, 
                  mean = output$mu_q_mu, sd = sqrt(output$sigmasq_q_mu)),
    h_sigmasq = dgamma(1/sigmasq, 
                        output$A_q_sigmasq, output$B_q_sigmasq)
                  
)

df_elbos <- dplyr::tibble(
    iteration = 1:10,
    elbos = output$elbos
)

p1 <- ggplot(df, aes(x=mu, y=h_mu)) +
    geom_line()

p2 <- ggplot(df, aes(x=sigmasq, y=h_sigmasq)) +
    geom_line()

p3 <- ggplot(df_elbos, aes(x = iteration, y = elbos)) + 
    geom_point() + 
    geom_line()

p3 + p1 + p2

mean(X)
