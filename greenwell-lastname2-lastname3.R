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
  geom_point(alpha = 0.5) +
  geom_smooth(color = "red", se = FALSE) +
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
  geom_point(alpha = 0.5) +
  geom_smooth(color = "red", se = FALSE) +
  xlab("x") +
  ylab("Probability-scale residual")

pdf(file = "manuscript\\quadratic.pdf", width = 8, height = 4)
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
  geom_point(size = 2, alpha = 0.25) +
  geom_smooth(col = "red", se = FALSE) +
  ylab("Probability scale residual")

# Figure ?
pdf(file = "heteroscedasticity.pdf", width = 8, height = 4)
grid.arrange(p1, p2, ncol = 2)
dev.off()


################################################################################
# Checking the proportionality assumption
################################################################################

# Fit separate models to the df4 data set and genrate the difference in 
# surrogate values
library(VGAM)
fit1 <- vglm(y ~ x, data = df4[1:2000, ], 
             cumulative(link = probit, parallel = TRUE))
fit2 <- update(fit1, data = df4[2001:4000, ])
s1 <- surrogate(fit1)
s2 <- surrogate(fit2)
d <- data.frame(D = s1 - s2, X = df4[1:2000, ]$x)

# Scatterplot of D vs. X
p <- ggplot(d, aes(x = X, y = D)) +
  geom_point() +
  geom_smooth(col = "red", se = FALSE)

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


################################################################################
# Quality of wine
################################################################################

library(ordinal)
data(wine, package = "ordinal")
wine.clm <- clm(rating ~ temp * contact, data = wine)  # default logit link

# Figure ?
pdf(file = "manuscript\\wine.pdf", width = 8, height = 8)
set.seed(101)  # for reproducibility
grid.arrange(
  autoplot(wine.clm, nsim = 10, what = "qq"),
  autoplot(wine.clm, nsim = 10, what = "fitted"),
  autoplot(wine.clm, nsim = 10, what = "cov", x = wine$temp),
  autoplot(wine.clm, nsim = 10, what = "cov", x = wine$contact),
  ncol = 2
)
dev.off()
