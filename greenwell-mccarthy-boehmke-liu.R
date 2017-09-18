################################################################################
# Setup
################################################################################

# Load required packages
library(ggplot2)
library(MASS)
library(ordinal)
library(PResiduals)
library(rms)
library(VGAM)
library(sure)


################################################################################
# Detecting a misspecified mean structure
################################################################################

# Load the simulated quadratic data
data(df1)

# Fit a (correct) probit model
fit.polr <- polr(y ~ x + I(x ^ 2), data = df1, method = "probit")

# Probability-scale residuals
pres <- presid(fit.polr)

# Residual plots using the SBS residuals
p1 <- ggplot(data.frame(x = df1$x, y = pres), aes(x, y)) +
  geom_point(color = "#444444", shape = 19, size = 2, alpha = 0.5) +
  geom_smooth(se = FALSE, size = 1.2, color = "red") +
  ylab("Probability residual")
p2 <- ggplot(data.frame(y = pres), aes(sample = y)) +
  stat_qq(distribution = qunif, dparams = list(min = -1, max = 1), alpha = 0.5) +
  xlab("Sample quantile") +
  ylab("Theoretical quantile")

# Figure ?
pdf(file = "quadratic-correct-sbs.pdf", width = 8, height = 4)
grid.arrange(p1, p2, ncol = 2)
dev.off()

# Surrogate residuals
set.seed(101)  # for reproducibility
sres <- resids(fit.polr)

# Residual plots using the surrogate-based residuals
p1 <- autoplot(sres, what = "covariate", x = df1$x, xlab = "x")
p2 <- autoplot(sres, what = "qq", disttribution = pnorm)

# Figure ?
pdf(file = "quadratic-correct-surrogate.pdf", width = 8, height = 4)
grid.arrange(p1, p2, ncol = 2)
dev.off()


fit2.polr <- update(fit.polr, y ~ x)
p1 <- autoplot(fit2.polr, what = "covariate", x = df1$x, alpha = 0.5) +
  xlab("x") +
  ylab("Surrogate residual") +
  ggtitle("")
p2 <- ggplot(data.frame(x = df1$x, y = presid(fit2.polr)), aes(x, y)) +
  geom_point(color = "#444444", shape = 19, size = 2, alpha = 0.5) +
  geom_smooth(se = FALSE, size = 1.2, color = "red") +
  xlab("x") +
  ylab("Probability-scale residual")

pdf(file = "quadratic.pdf", width = 8, height = 4)
grid.arrange(p1, p2, ncol = 2)
dev.off()


################################################################################
# Detecting heteroscedasticty
################################################################################

# Fit a cumulative link model with probit link
fit.orm <- orm(y ~ x, data = df2, family = "probit", x = TRUE)

# Residual vs. covariate plots
set.seed(102)  # for reproducibility
p1 <- autoplot(fit.orm, what = "covariate", x = df2$x, xlab = "x")
p2 <- ggplot(data.frame(x = df2$x, y = presid(fit.orm)), aes(x, y)) +
  geom_point(color = "#444444", shape = 19, size = 2, alpha = 0.25) +
  geom_smooth(se = FALSE, size = 1.2, color = "red") +
  ylab("Probability scale residual")

# Figure ?
pdf(file = "heteroscedasticity.pdf", width = 8, height = 4)
grid.arrange(p1, p2, ncol = 2)
dev.off()

# Fit a VGAM (i.e., nonparametric model) with probit link
fit.vgam <- vgam(y ~ s(x), family = cumulative(link = probit, parallel = TRUE),
                 data = df2)

# Residual vs. covariate plots: jittering
set.seed(103)
p1 <- autoplot(fit.vgam, what = "covariate", x = df2$x, method = "jitter",
               xlab = "x")
p2 <- autoplot(fit.vgam, what = "covariate", x = df2$x, method = "jitter",
               jitter.scale = "response", xlab = "x")

# Figure?
pdf(file = "heteroscedasticity2.pdf", width = 8, height = 4)
grid.arrange(p1, p2, ncol = 2)
dev.off()


################################################################################
# Checking the proportionality assumption
################################################################################

ordinalize <- function(z, threshold) {
  sapply(z, FUN = function(x) {
    ordinal.value <- 1
    index <- 1
    while(index <= length(threshold) && x > threshold[index]) {
      ordinal.value <- ordinal.value + 1
      index <- index + 1
    }
    ordinal.value
  })
}

# Function to simulate the data from Example 5 in Dungang and Zhang (2017).
simProportionalityData <- function(n = 2000) {
  x <- runif(n, min = -3, max = 3)
  z1 <- 0 - 1 * x + rnorm(n)
  z2 <- 0 - 1.5 * x + rnorm(n)
  y1 <- ordinalize(z1, threshold = c(-1.5, 0))
  y2 <- ordinalize(z2, threshold = c(1, 3))
  data.frame("y" = as.ordered(c(y1, y2)), "x" = c(x, x))
}

# Simulate data
set.seed(977)
df4 <- simProportionalityData(n = 2000)
table(df4$y)

# Fit separate models to the df4 data set and genrate the difference in 
# surrogate values
fit1 <- vglm(y ~ x, data = df4[1:2000, ], 
             cumulative(link = probit, parallel = TRUE))
fit2 <- update(fit1, data = df4[2001:4000, ])
s1 <- surrogate(fit1)
s2 <- surrogate(fit2)
d <- data.frame(D = s1 - s2, x = df4[1:2000, ]$x)

# Scatterplot of D vs. x
p <- ggplot(d, aes(x = x, y = D)) +
  geom_point(color = "#444444", shape = 19, size = 2) +
  geom_smooth(se = FALSE, size = 1.2, color = "red")

# Figure ?
pdf(file = "proportionality.pdf", width = 7, height = 5)
print(p)
dev.off()


################################################################################
# Detecting a misspecified link function
################################################################################

# Fit models with various link functions to the simulated data
fit.probit <- polr(y ~ x + I(x ^ 2), data = df3, method = "probit")
fit.logistic <- polr(y ~ x + I(x ^ 2), data = df3, method = "logistic")
fit.loglog <- polr(y ~ x + I(x ^ 2), data = df3, method = "loglog")  # correct link
fit.cloglog <- polr(y ~ x + I(x ^ 2), data = df3, method = "cloglog")

# Construc Q-Q plots of the surrogate residuals for each model
p1 <- autoplot(fit.probit, nsim = 100, what = "qq")
p2 <- autoplot(fit.logistic, nsim = 100, what = "qq")
p3 <- autoplot(fit.loglog, nsim = 100, what = "qq")
p4 <-  autoplot(fit.cloglog, nsim = 100, what = "qq")

# Figure ?
pdf(file = "link.pdf", width = 7, height = 7)
grid.arrange(p1, p2, p3, p4, ncol = 2)  # bottom left plot is correct model
dev.off()

# Figure ?
pdf(file = "gof.pdf", width = 7, height = 7)
par(mfrow = c(2, 2), mar = c(4, 4, 2, 2) + 0.1) 
set.seed(8491)  # for reproducibility
plot(gof(fit.probit, nsim = 100, test = "ad"), main = "")
plot(gof(fit.logistic, nsim = 100, test = "ad"), main = "")
plot(gof(fit.loglog, nsim = 100, test = "ad"), main = "")
plot(gof(fit.cloglog, nsim = 100, test = "ad"), main = "")
dev.off()


################################################################################
# Interaction detection
################################################################################

# Function to simulate data from an ordered probit model with an interaction
# term
simInteractionData <- function(n = 2000) {
  threshold <- c(0, 20, 40)
  x1 <- runif(n, min = 1, max = 7)
  x2 <- gl(2, n / 2, labels = c("Control", "Treatment"))
  z <- 16 - 5 * x1 + 3 * (x2 == "Treatment") + 10 * x1 * (x2 == "Treatment") + 
    rnorm(n)
  y <- sapply(z, FUN = function(zz) {
    ordinal.value <- 1
    index <- 1
    while(index <= length(threshold) && zz > threshold[index]) {
      ordinal.value <- ordinal.value + 1
      index <- index + 1
    }
    ordinal.value
  })
  data.frame("y" = as.ordered(y), "x1" = x1, "x2" = x2)
}

# Simulate data
set.seed(977)
df4 <- simInteractionData(n = 2000)

library(ggplot2)
library(ordinal)
library(sure)

fit1 <- clm(y ~ x1, data = df4[df4$x2 == "Control", ], link = "probit")
fit2 <- clm(y ~ x1, data = df4[df4$x2 == "Treatment", ], link = "probit")
fit3 <- clm(y ~ x1*x2, data = df4, link = "probit")

set.seed(1105)
d1 <- cbind(df4[df4$x2 == "Control",], sur = surrogate(fit1, nsim = 25))
d2 <- cbind(df4[df4$x2 == "Treatment", ], sur = surrogate(fit2, nsim = 25))
p1 <- ggplot(d1, aes(x = x1, y = sur)) +
  geom_point(color = "#444444", shape = 19, size = 2, alpha = 0.5) +
  geom_smooth(se = FALSE, size = 1.2, color = "red") +
  ylab("Surrogate response") +
  xlab(expression(paste(x[1], " (control)")))
p2 <- ggplot(d2, aes(x = x1, y = sur)) +
  geom_point(color = "#444444", shape = 19, size = 2, alpha = 0.5) +
  geom_smooth(se = FALSE, size = 1.2, color = "red") +
  ylab("Surrogate response") + 
  xlab(expression(paste(x[1], " (treatment)")))
pdf(file = "interaction.pdf", width = 8, height = 4)
grid.arrange(p1, p2, ncol = 2)
dev.off()


################################################################################
# Quality of wine
################################################################################

data(wine, package = "ordinal")
wine.clm <- clm(rating ~ temp + contact, data = wine, link = "probit")

# Figure ?
pdf(file = "wine.pdf", width = 7, height = 5)
set.seed(1225)  # for reproducibility
grid.arrange(
  autoplot(wine.clm, nsim = 10, what = "qq"),
  autoplot(wine.clm, nsim = 10, what = "fitted", alpha = 0.5),
  autoplot(wine.clm, nsim = 10, what = "covariate", x = wine$temp, 
           xlab = "Temperature"),
  autoplot(wine.clm, nsim = 10, what = "covariate", x = wine$contact,
           xlab = "Contact"),
  ncol = 2
)
dev.off()
