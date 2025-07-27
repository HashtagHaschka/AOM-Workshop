# =============================================
# WORKSHOP: Copula-Based Endogeneity Correction
# =============================================
# https://github.com/HashtagHaschka/Copula-based-endogeneity-corrections
# =============================================


source("https://raw.githubusercontent.com/HashtagHaschka/AOM-Workshop/functions/CopReg2sCOPE-np_workshop.R")
source("https://raw.githubusercontent.com/HashtagHaschka/AOM-Workshop/functions/CopReg2sCOPE_workshop.R")
source("https://raw.githubusercontent.com/HashtagHaschka/AOM-Workshop/functions/CopRegBMW_workshop.R")
source("https://raw.githubusercontent.com/HashtagHaschka/AOM-Workshop/functions/CopRegIMA_workshop.R")
source("https://raw.githubusercontent.com/HashtagHaschka/AOM-Workshop/functions/CopRegJAMS_workshop.R")
source("https://raw.githubusercontent.com/HashtagHaschka/AOM-Workshop/functions/CopRegPG_workshop.R")
source("https://raw.githubusercontent.com/HashtagHaschka/AOM-Workshop/functions/ICAreg_workshop.R")


# --------------------------- #
#  ENDOGENEITY DEMONSTRATION  #
# --------------------------- #

set.seed(123)
n <- 1000
rho <- 0.8  # Endogeneity strength

# Generate endogenous data (DGP from Slide 4)
latent <- mvtnorm::rmvnorm(n, mean = c(0, 0), 
                           sigma = matrix(c(1, rho, rho, 1), ncol = 2))
P <- qlnorm(pnorm(latent[, 1]))  # Endogenous regressor
xi <- latent[, 2]                # Error term
X <- rnorm(n)                    # Exogenous control
Y <- 2 + 1.5*X + 3*P + xi        # True model
data1 <- as.data.frame(cbind(Y, P, X))

# Biased OLS estimation
ols_model <- lm(Y ~ X + P, data1)
summary(ols_model)  
# Coefficient of P is biased
# Intercept is not biased
# Coefficient of X is not biased



# ------------------------------ #
#  PARK-GUPTA (2012) CORRECTION  #
# ------------------------------ #

# Estimate with Park-Gupta
pg_model <- CopRegPG(formula = Y ~ P | X,
                     data = data1, cdf = "ecdf")
pg_model[[1]]
# not biased

# Results comparison
results_pg <- data.frame(
  Method = c("True", "Biased OLS", "Park-Gupta"),
  Intercept = c(2, coef(ols_model)[1], pg_model[[1]][1, 1]),
  Beta_X = c(1.5, coef(ols_model)[2], pg_model[[1]][3, 1]),
  Alpha_P = c(3, coef(ols_model)[3], pg_model[[1]][2, 1])
)
print(results_pg)

# Diagnostic plot: Original vs. Copula-transformed P
P_star <- qnorm(pobs(P))
par(mfrow = c(1,2))
hist(P, main = "Original P (Non-normal)", breaks = 30)
hist(P_star, main = "Copula-transformed P (Standard Normal)", breaks = 30)
par(mfrow = c(1,1))
plot(P_star ~ P, xlab = "Original P (Non-normal)", ylab = "Copula-transformed P (Standard Normal)")



# -------------------------------------------------- #
#  HASCHKA (2024) and YANG ET AL. (2025) CORRECTION  #
# -------------------------------------------------- #

# Simulate data with correlated regressors (Slide 17)
set.seed(456)
n <- 1000
rho <- 0.8  # Endogeneity strength
r <- 0.5  # P-X correlation
Sigma <- matrix(c(1, r, rho,
                  r, 1, 0,
                  rho, 0, 1), ncol=3)
latent <- mvtnorm::rmvnorm(n, sigma = Sigma)
P <- qexp(pnorm(latent[,1]), rate = 0.5)
X <- latent[, 2]
xi <- latent[, 3]
Y <- 2 + 1.5*X + 3*P + xi
data1 <- as.data.frame(cbind(Y, P, X))

# Biased OLS estimation
ols_model <- lm(Y ~ X + P, data1)
summary(ols_model)  # Coefficients of both X and P are biased

# Estimate with Haschka (2024) method
haschka_model <- CopRegIMA(formula = Y ~ P | X,
                          data = data1, cdf = "ecdf")
haschka_model[[1]]

# Estimate with 2sCOPE method
twoscope_model <- CopReg2sCOPE(formula = Y ~ P | X,
                              data = data1, cdf = "ecdf")
twoscope_model[[1]]

# Estimate with Park-Gupta
pg_model <- CopRegPG(formula = Y ~ P | X,
                     data = data1, cdf = "ecdf")
pg_model[[1]]


# Results comparison
results_two <- data.frame(
  Method = c("True", "Biased OLS", "Park-Gupta", "Haschka_IMA", "2sCOPE"),
  Intercept = c(2, coef(ols_model)[1], pg_model[[1]][1, 1], haschka_model[[1]][1, 1], twoscope_model[[1]][1, 1]), 
  Beta_X = c(1.5, coef(ols_model)[2], pg_model[[1]][3, 1], haschka_model[[1]][3, 1], twoscope_model[[1]][3, 1]), 
  Alpha_P = c(3, coef(ols_model)[3], pg_model[[1]][2, 1], haschka_model[[1]][2, 1], twoscope_model[[1]][2, 1])
)
print(results_two)



# ------------------------------------------------- #
#  LIENGAARD ET AL. (2025) - CATEGORICAL MODERATOR  #
# ------------------------------------------------- #

# Simulate data with endogeneity varying by group
set.seed(789)
n <- 1000
rho1 <- 0.4  # Endogeneity strength
rho2 <- 0.8  # Endogeneity strength
X <- sample(0:1, n, replace = TRUE)
P <- numeric(n)
xi <- numeric(n)

latent1 <- mvtnorm::rmvnorm(n/2, mean = c(0, 0), 
                           sigma = matrix(c(1, rho1, rho1, 1), ncol = 2))
latent2 <- mvtnorm::rmvnorm(n/2, mean = c(0, 0), 
                           sigma = matrix(c(1, rho2, rho2, 1), ncol = 2))

P <- c(qexp(pnorm(latent1[, 1]), rate = 0.5), 
       qunif(pnorm(latent2[, 1])))
X <- c(rep(0, n/2),
       rep(1, n/2))
xi <- c(latent1[, 2], 
        latent2[, 2])
Y <- 2 + 1.5*X + 3*P + xi
data1 <- as.data.frame(cbind(Y, P, X))

# Biased OLS estimation
ols_model <- lm(Y ~ X + P, data1)
summary(ols_model)  # Coefficients of both X and P are biased

# Estimate with Park-Gupta
pg_model <- CopRegPG(formula = Y ~ P | X,
                     data = data1, cdf = "ecdf")
pg_model[[1]]

# Estimate with Liengaard et al.
liengaard_model <- CopRegJAMS(formula = Y ~ P | as.factor(X),
                       data = data1, cdf = "ecdf")
liengaard_model[[1]]

# Estimate with Hu et al.
hu_model <- CopReg2sCOPEnp(formula = Y ~ P | as.factor(X),
                           data = data1)
hu_model[[1]]



# ------------------------- #
#  REAL DATA DEMONSTRATION  #
# ------------------------- #

library(bayesm)
library(dplyr)

data("orangeJuice")
dat1 <- orangeJuice[[1]]
help(orangeJuice)

dat1_FloridaNatural <- dat1 %>% filter(brand == 3)
dat1_FloridaNatural <- dat1_FloridaNatural %>%
  mutate(across(starts_with("price"), log))


# Fit least squares
mod1 <- lm(formula = logmove ~ price3 + feat + as.factor(deal),
           data = dat1_FloridaNatural)

# Fit Park & Gupta (2012)
mod_PG <- CopRegPG(...)
mod_PG[[1]]

# Fit Haschka (2024)
mod_IMA <- CopRegIMA(...)
mod_IMA[[1]]

# Fit 2sCOPE
mod_2sCOPE <- CopReg2sCOPE(...)
mod_2sCOPE[[1]]
  
# Fit Liengaard et al. (2025)
mod_JAMS <- CopRegJAMS(...)
mod_JAMS[[1]]

# Fit Hu et al. (2025)
mod_2sCOPEnp <- CopReg2sCOPEnp(...)
mod_2sCOPEnp[[1]]

# Fit ICA
mod1 <- ica_reg(formula = logmove ~ price3 | feat + as.factor(deal), 
                data = dat1_FloridaNatural)
mod1[[1]] # estimates



# ------------------- #
#  OMITTED VARIABLES  #
# ------------------- #

n <- 1000
W <- rnorm(n) # omitted variable
P <- W + (rgamma(n, shape = 1, rate = 1) - 1) # observed regressor
xi <- rnorm(n) # normal error

Y <- 3 + 1 * P + 1 * W + xi
data1 <- as.data.frame(cbind(Y, P, W))

mod1 <- lm(Y ~ P, data1) # variable W is omitted
summary(mod1) # bias because of omitted variable


# Estimate with Park-Gupta
pg_model <- CopRegPG(formula = Y ~ P,
                     data = data1, cdf = "ecdf")
pg_model[[1]] # still biased


mod1 <- ica_reg(formula = Y ~ P, data = data1)
mod1[[1]] # works


### AN EXAMPLE FOR 2sCOPE

n <- 1000
W <- rnorm(n) # omitted variable
X <- (rchisq(n, df = 1) - 1) # exogenous regressor
P <- X + W + (rgamma(n, shape = 1, rate = 1) - 1) # observed regressor
xi <- rnorm(n) # normal error

Y <- 3 + 1 * P + 1 * W + 1 * X + xi
data1 <- as.data.frame(cbind(Y, P, W, X))

mod1 <- lm(Y ~ P + X, data1) # variable W is omitted
summary(mod1) # bias because of omitted variable


# Estimate with 2sCOPE
scope_model <- CopReg2sCOPE(formula = Y ~ P | X,
                     data = data1, cdf = "ecdf")
scope_model[[1]] # still biased

mod1 <- ica_reg(formula = Y ~ P | X, data = data1)
mod1[[1]] # works

