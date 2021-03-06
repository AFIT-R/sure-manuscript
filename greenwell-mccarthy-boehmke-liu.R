################################################################################
# Detecting a misspecified mean structure
################################################################################

# Fit a cumulative link model with probit link
library(sure)  # for residual function and sample data sets
library(MASS)  # for polr function
fit.polr <- polr(y ~ x + I(x ^ 2), data = df1, method = "probit")

# Obtain the SBS/probability-scale residuals
library(PResiduals)
pres <- presid(fit.polr)

# Residual vs. covariate plot and Q-Q plot
library(ggplot2)  # for plotting
p1 <- ggplot(data.frame(x = df1$x, y = pres), aes(x, y)) +
  geom_point(color = "#444444", shape = 19, size = 2, alpha = 0.5) +
  geom_smooth(color = "red", se = FALSE) +
  ylab("Probability-scale residual")
p2 <- ggplot(data.frame(y = pres), aes(sample = y)) +
  stat_qq(distribution = qunif, dparams = list(min = -1, max = 1), alpha = 0.5) +
  xlab("Sample quantile") +
  ylab("Theoretical quantile")
pdf(file = "quadratic-correct-sbs.pdf", width = 8, height = 4)
grid.arrange(p1, p2, ncol = 2)  # Figure 1
dev.off()

# Obtain surrogate residuals
set.seed(101)  # for reproducibility
sres <- resids(fit.polr)

# Residual vs. covariate plot and Q-Q plot
p1 <- autoplot(sres, what = "covariate", x = df1$x, xlab = "x")
p2 <- autoplot(sres, what = "qq", distribution = qnorm)
pdf(file = "quadratic-correct-surrogate.pdf", width = 8, height = 4)
grid.arrange(p1, p2, ncol = 2)  # Figure 2
dev.off()

set.seed(101)  # for reproducibility
autoplot(fit.polr, what = "qq")  # same as top right of Figure 1

# Residual vs. covariate plots when quadratic term is removed from the model
fit2.polr <- update(fit.polr, y ~ x)
set.seed(1055)
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
grid.arrange(p1, p2, ncol = 2)  # Figure 3
dev.off()


################################################################################
# Detecting heteroscedasticty
################################################################################

# Fit a cumulative link model with probit link
library(rms)  # for orm function
fit.orm <- orm(y ~ x, data = df2, family = "probit", x = TRUE)

pres <- presid(fit.orm)  # SBS residuals
set.seed(102)  # for reproducibility
sres <- resids(fit.orm)  # surrogate residuals

# Residual vs. covariate plots
p1 <- autoplot(sres, what = "covariate", x = df2$x, xlab = "x")
p2 <- ggplot(data.frame(x = df2$x, y = presid(fit.orm)), aes(x, y)) +
  geom_point(color = "#444444", shape = 19, size = 2, alpha = 0.25) +
  geom_smooth(col = "red", se = FALSE) +
  ylab("Probability scale residual")
pdf(file = "heteroscedasticity.pdf", width = 8, height = 4)
grid.arrange(p1, p2, ncol = 2)  # Figure 4
dev.off()


library(VGAM)  # for vgam and vglm functions
fit.vgam <- vgam(y ~ s(x), family = cumulative(link = probit, parallel = TRUE),
                 data = df2)

set.seed(103)  # for reproducibility
p1 <- autoplot(fit.vgam, what = "covariate", x = df2$x, method = "jitter",
               xlab = "x")
p2 <- autoplot(fit.vgam, what = "covariate", x = df2$x, method = "jitter",
               jitter.scale = "response", xlab = "x")
pdf(file = "heteroscedasticity2.pdf", width = 8, height = 4)
grid.arrange(p1, p2, ncol = 2)  # Figure 5
dev.off()


################################################################################
# Detecting a misspecified link function
################################################################################

# Fit models with various link functions to the simulated data
fit.probit <- polr(y ~ x + I(x ^ 2), data = df3, method = "probit")
fit.logistic <- polr(y ~ x + I(x ^ 2), data = df3, method = "logistic")
fit.loglog <- polr(y ~ x + I(x ^ 2), data = df3, method = "loglog")  # correct link
fit.cloglog <- polr(y ~ x + I(x ^ 2), data = df3, method = "cloglog")

# Construct Q-Q plots of the surrogate residuals for each model
set.seed(1056)  # for reproducibility
p1 <- autoplot(fit.probit, nsim = 100, what = "qq")
p2 <- autoplot(fit.logistic, nsim = 100, what = "qq")
p3 <- autoplot(fit.loglog, nsim = 100, what = "qq")
p4 <-  autoplot(fit.cloglog, nsim = 100, what = "qq")

# Figure 6
pdf(file = "link.pdf", width = 7, height = 7)
grid.arrange(p1, p2, p3, p4, ncol = 2)  # bottom left plot is correct model
dev.off()

# Figure 7
pdf(file = "gof.pdf", width = 7, height = 7)
par(mfrow = c(2, 2), mar = c(2, 4, 2, 2) + 0.1)
set.seed(8491)  # for reproducibility
plot(gof(fit.probit, nsim = 100, test = "ad"), main = "")
plot(gof(fit.logistic, nsim = 100, test = "ad"), main = "")
plot(gof(fit.loglog, nsim = 100, test = "ad"), main = "")
plot(gof(fit.cloglog, nsim = 100, test = "ad"), main = "")
dev.off()


################################################################################
# Checking the proportionality assumption
################################################################################

# Fit separate models (VGAM should already be loaded)
fit1 <- vglm(y ~ x, data = df4[1:2000, ],
             cumulative(link = probit, parallel = TRUE))
fit2 <- update(fit1, data = df4[2001:4000, ])

# Generate surrogate response values
set.seed(8671)  # for reproducibility
s1 <- surrogate(fit1)
s2 <- surrogate(fit2)

# Figure 8
pdf(file = "proportionality.pdf", width = 7, height = 5)
ggplot(data.frame(D = s1 - s2, x = df4[1:2000, ]$x) , aes(x = x, y = D)) +
  geom_point(color = "#444444", shape = 19, size = 2) +
  geom_smooth(se = FALSE, size = 1.2, color = "red")
dev.off()


################################################################################
# Interaction detection
################################################################################

library(ordinal)  # for clm function
fit1 <- clm(y ~ x1, data = df5[df5$x2 == "Control", ], link = "probit")
fit2 <- clm(y ~ x1, data = df5[df5$x2 == "Treatment", ], link = "probit")

set.seed(1105)  # for reproducibility
d1 <- cbind(df5[df5$x2 == "Control",], sur = surrogate(fit1, nsim = 25))
d2 <- cbind(df5[df5$x2 == "Treatment", ], sur = surrogate(fit2, nsim = 25))
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
grid.arrange(p1, p2, ncol = 2)  # Figure 9
dev.off()


################################################################################
# Quality of wine
################################################################################

data(wine, package = "ordinal")  # load wine data set
wine.clm <- clm(rating ~ temp + contact, data = wine, link = "probit")

pdf(file = "wine.pdf", width = 7, height = 5)
set.seed(1225)  # for reproducibility
grid.arrange(  # Figure 10
  autoplot(wine.clm, nsim = 10, what = "qq"),
  autoplot(wine.clm, nsim = 10, what = "fitted", alpha = 0.5),
  autoplot(wine.clm, nsim = 10, what = "covariate", x = wine$temp,
           xlab = "Temperature"),
  autoplot(wine.clm, nsim = 10, what = "covariate", x = wine$contact,
           xlab = "Contact"),
  ncol = 2
)
dev.off()
