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
