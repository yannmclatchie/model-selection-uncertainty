library(dplyr)
library(ggplot2)
library(ggExtra)
library(loo)
library(brms)

# sample the data
x1 = rnorm(0, 1)
x2 = rnorm(0, 1)
x3 = rnorm(0, 1)
x4 = rnorm(0, 1)
x5 = rnorm(0, 1)
x6 = rnorm(0, 1)
x7 = rnorm(0, 1)
x8 = rnorm(0, 1)
x9 = rnorm(0, 1)
x10 = rnorm(0, 1)
x11 = rnorm(0, 1)
y = rnorm(x11, 1)
data <- data.frame(
  cbind(y, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11)
)

# fit the models
fit1 <- brm(y ~ x1, family = gaussian(), data = data)
fit2 <- update(fit1, formula. = y ~ x2)
fit3 <- update(fit1, formula. = y ~ x3)
fit4 <- update(fit1, formula. = y ~ x4)
fit5 <- update(fit1, formula. = y ~ x5)
fit6 <- update(fit1, formula. = y ~ x6)
fit7 <- update(fit1, formula. = y ~ x7)
fit8 <- update(fit1, formula. = y ~ x8)
fit9 <- update(fit1, formula. = y ~ x9)
fit10 <- update(fit1, formula. = y ~ x10)
fit11 <- update(fit1, formula. = y ~ x11)

# compute model elpds
