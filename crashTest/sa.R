
#### Aim of prog: Test sensitivity analysis
## Comments:
# I create a particular case with a linear regression Y = b0 + b1 X1 + b2 X2
#	where the amplitude of Y is the same with respect to X1 and X2, but the
#	slopes and ranges differ

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(sensitivity)
library(data.table)
library(sensobol)

#### Create data
## X domain
x1_min = -1.02
x1_max = 12.4

x2_min = 5
x2_max = 5.25 # To force a value of beta2 = 4*alpha, much larger than beta1 in this case

## Regression coefficients
alpha = 4.56 # amplitude
beta0 = -1.25
(beta1 = alpha/(x1_max - x1_min))
(beta2 = alpha/(x2_max - x2_min))

## Y domain
y1_min = beta1*x1_min + beta0
y1_max = y1_min + alpha

y2_min = beta2*x2_min + beta0
y2_max = y2_min + alpha

all.equal(y1_max, beta1*x1_max + beta0)
all.equal(y2_max, beta2*x2_max + beta0)

#### Plots
x1 = seq(x1_min, x1_max, length = 30)
x2 = seq(x2_min, x2_max, length = 30)
y = outer(x1, x2, function(a, b) a*b^2)

# persp(x1, x2, y, theta = 30, phi = 30, expand = 0.5, col = "lightblue", ticktype = "detailed")

y1_max - y1_min
y2_max - y2_min

#? ---------------------------------------------------------------------------------------
#* ----------------------    PART I: Sensitivity analysis x1, x2    ----------------------
#? ---------------------------------------------------------------------------------------
#### Sensitivity analysis
## Common variables
explanatory_vars = c("x1", "x2")
matrices = c("A", "B", "AB")
order = "first"
N = 2^14
type = "QRN"
first = "saltelli"
total = "jansen"

## Sobol matrices
sobol_mat = sobol_matrices(matrices = matrices, N = N, order = order, type = type, params = explanatory_vars)

# Rescale the Sobol matrices
sobol_mat[, "x1"] = qunif(p = sobol_mat[, "x1"], min = x1_min, max = x1_max)
sobol_mat[, "x2"] = qunif(p = sobol_mat[, "x2"], min = x2_min, max = x2_max)

y = beta0 + beta1*sobol_mat[, "x1"] + beta2*sobol_mat[, "x2"]

ind_data = sobol_indices(matrices = matrices, Y = y, N = N, params = explanatory_vars, first = first, total = total, order = order)
ind_data$results[sensitivity == "Si"]
print(ind_data$results[sensitivity == "Si"], digits = 2)

#     original sensitivity parameters
# 1: 0.5049740          Si         x1
# 2: 0.5042093          Si         x2

#? ---------------------------------------------------------------------------------------------------
#* ----------------------    PART II: Sensitivity analysis x1, x2, reg coeff    ----------------------
#? ---------------------------------------------------------------------------------------------------
#### Sensitivity analysis, uncertainty parameters is a fraction of their absolute value
## Common variables
explanatory_vars = c("beta0", "beta1", "beta2", "x1", "x2")
matrices = c("A", "B", "AB")
order = "first"
N = 2^18
type = "QRN"
first = "saltelli"
total = "jansen"

## Sobol matrices
sobol_mat = sobol_matrices(matrices = matrices, N = N, order = order, type = type, params = explanatory_vars)

# Rescale the Sobol matrices
sobol_mat[, "beta0"] = qnorm(p = sobol_mat[, "beta0"], mean = beta0, sd = abs(beta0)/10)
sobol_mat[, "beta1"] = qnorm(p = sobol_mat[, "beta1"], mean = beta1, sd = abs(beta1)/10)
sobol_mat[, "beta2"] = qnorm(p = sobol_mat[, "beta2"], mean = beta2, sd = abs(beta2)/10)
sobol_mat[, "x1"] = qunif(p = sobol_mat[, "x1"], min = x1_min, max = x1_max)
sobol_mat[, "x2"] = qunif(p = sobol_mat[, "x2"], min = x2_min, max = x2_max)

y = sobol_mat[, "beta0"] + sobol_mat[, "beta1"]*sobol_mat[, "x1"] + sobol_mat[, "beta2"]*sobol_mat[, "x2"]

ind_all_frac = sobol_indices(matrices = matrices, Y = y, N = N, params = explanatory_vars, first = first, total = total, order = order)
ind_all_frac$results[sensitivity == "Si"]
print(ind_all_frac$results[sensitivity == "Si"], digits = 2)
ind_all_frac$results[sensitivity == "Si", sum(original)]

#        original sensitivity parameters
# 1: 0.0001735519          Si      beta0
# 2: 0.0004105937          Si      beta1
# 3: 0.9609832773          Si      beta2
# 4: 0.0190604162          Si         x1
# 5: 0.0190642330          Si         x2

#### Sensitivity analysis, uncertainty parameters is a set value
## Sobol matrices
sobol_mat = sobol_matrices(matrices = matrices, N = N, order = order, type = type, params = explanatory_vars)

# Rescale the Sobol matrices
sobol_mat[, "beta0"] = qnorm(p = sobol_mat[, "beta0"], mean = beta0, sd = 1)
sobol_mat[, "beta1"] = qnorm(p = sobol_mat[, "beta1"], mean = beta1, sd = 1)
sobol_mat[, "beta2"] = qnorm(p = sobol_mat[, "beta2"], mean = beta2, sd = 1)
sobol_mat[, "x1"] = qunif(p = sobol_mat[, "x1"], min = x1_min, max = x1_max)
sobol_mat[, "x2"] = qunif(p = sobol_mat[, "x2"], min = x2_min, max = x2_max)

y = sobol_mat[, "beta0"] + sobol_mat[, "beta1"]*sobol_mat[, "x1"] + sobol_mat[, "beta2"]*sobol_mat[, "x2"]

ind_all_set = sobol_indices(matrices = matrices, Y = y, N = N, params = explanatory_vars, first = first, total = total, order = order)
ind_all_frac$results[sensitivity == "Si"]
print(ind_all_set$results[sensitivity == "Si"], digits = 2)

#      original sensitivity parameters
# 1: 0.01282023          Si      beta0
# 2: 0.41445380          Si      beta1
# 3: 0.33623771          Si      beta2
# 4: 0.02219121          Si         x1
# 5: 0.02219246          Si         x2

#? --------------------------------------------------------------------------------------------------------------------------
#* ----------------------    PART III: Sensitivity analysis x1, x2, reg coeff, different amplitudes    ----------------------
#? --------------------------------------------------------------------------------------------------------------------------
#### Create data
## X domain
x1_min = -1.02
x1_max = 12.4

x2_min = 5
x2_max = 5.25 # Put 25.25 for instance, to get beta2 = 0.11, which is three times smaller than beta1 and see the effect on S.A.!

## Regression coefficients
alpha = 4.56 # amplitude
beta0 = -1.25
(beta1 = alpha/(x1_max - x1_min))
(beta2 = alpha/(2*(x2_max - x2_min)))

## Y domain
y1_min = beta1*x1_min + beta0
y1_max = y1_min + alpha

y2_min = beta2*x2_min + beta0
y2_max = y2_min + alpha/2

y1_max - y1_min
y2_max - y2_min

#### Sensitivity analysis
## Common variables
explanatory_vars = c("x1", "x2")
matrices = c("A", "B", "AB")
order = "first"
N = 2^14
type = "QRN"
first = "saltelli"
total = "jansen"

## Sobol matrices
sobol_mat = sobol_matrices(matrices = matrices, N = N, order = order, type = type, params = explanatory_vars)

# Rescale the Sobol matrices
sobol_mat[, "x1"] = qunif(p = sobol_mat[, "x1"], min = x1_min, max = x1_max)
sobol_mat[, "x2"] = qunif(p = sobol_mat[, "x2"], min = x2_min, max = x2_max)

y = beta0 + beta1*sobol_mat[, "x1"] + beta2*sobol_mat[, "x2"]

ind_data = sobol_indices(matrices = matrices, Y = y, N = N, params = explanatory_vars, first = first, total = total, order = order)
ind_data$results[sensitivity == "Si"]

#### Sensitivity analysis, uncertainty parameters is a fraction of their absolute value
## Common variables
explanatory_vars = c("beta0", "beta1", "beta2", "x1", "x2")
matrices = c("A", "B", "AB")
order = "first"
N = 2^18
type = "QRN"
first = "saltelli"
total = "jansen"

## Sobol matrices
sobol_mat = sobol_matrices(matrices = matrices, N = N, order = order, type = type, params = explanatory_vars)

# Rescale the Sobol matrices
sobol_mat[, "beta0"] = qnorm(p = sobol_mat[, "beta0"], mean = beta0, sd = abs(beta0)/10)
sobol_mat[, "beta1"] = qnorm(p = sobol_mat[, "beta1"], mean = beta1, sd = abs(beta1)/10)
sobol_mat[, "beta2"] = qnorm(p = sobol_mat[, "beta2"], mean = beta2, sd = abs(beta2)/10)
sobol_mat[, "x1"] = qunif(p = sobol_mat[, "x1"], min = x1_min, max = x1_max)
sobol_mat[, "x2"] = qunif(p = sobol_mat[, "x2"], min = x2_min, max = x2_max)

y = sobol_mat[, "beta0"] + sobol_mat[, "beta1"]*sobol_mat[, "x1"] + sobol_mat[, "beta2"]*sobol_mat[, "x2"]

ind_all_frac = sobol_indices(matrices = matrices, Y = y, N = N, params = explanatory_vars, first = first, total = total, order = order)
ind_all_frac$results[sensitivity == "Si"]
print(ind_all_frac$results[sensitivity == "Si"], digits = 2)
ind_all_frac$results[sensitivity == "Si", sum(original)]

#        original sensitivity parameters
# 1: 0.0006520513          Si      beta0
# 2: 0.0015511628          Si      beta1
# 3: 0.9070336207          Si      beta2
# 4: 0.0719540309          Si         x1
# 5: 0.0179940713          Si         x2

#### Sensitivity analysis, uncertainty parameters is a set value
## Sobol matrices
sobol_mat = sobol_matrices(matrices = matrices, N = N, order = order, type = type, params = explanatory_vars)

# Rescale the Sobol matrices
sobol_mat[, "beta0"] = qnorm(p = sobol_mat[, "beta0"], mean = beta0, sd = 1)
sobol_mat[, "beta1"] = qnorm(p = sobol_mat[, "beta1"], mean = beta1, sd = 1)
sobol_mat[, "beta2"] = qnorm(p = sobol_mat[, "beta2"], mean = beta2, sd = 1)
sobol_mat[, "x1"] = qunif(p = sobol_mat[, "x1"], min = x1_min, max = x1_max)
sobol_mat[, "x2"] = qunif(p = sobol_mat[, "x2"], min = x2_min, max = x2_max)

y = sobol_mat[, "beta0"] + sobol_mat[, "beta1"]*sobol_mat[, "x1"] + sobol_mat[, "beta2"]*sobol_mat[, "x2"]

ind_all_set = sobol_indices(matrices = matrices, Y = y, N = N, params = explanatory_vars, first = first, total = total, order = order)
ind_all_set$results[sensitivity == "Si"]

#       original sensitivity parameters
# 1: 0.013028621          Si      beta0
# 2: 0.421475138          Si      beta1
# 3: 0.341930697          Si      beta2
# 4: 0.022558875          Si         x1
# 5: 0.005641709          Si         x2

#? --------------------------------------------------------------------------------------------------------
#* ----------------------    PART VI: Comparing packages: sensobol vs sensitivity    ----------------------
#? ---------------------------------------------------------------------------------------------------------
#### Create data
## X domain
x1_min = -1.02
x1_max = 12.4

x2_min = 5
x2_max = 5.25 # Put 25.25 for instance, to get beta2 = 0.11, which is three times smaller than beta1 and see the effect on S.A.!

## Regression coefficients
alpha = 4.56 # amplitude
beta0 = -1.25
(beta1 = alpha/(x1_max - x1_min))
(beta2 = alpha/(2*(x2_max - x2_min)))

## Y domain
y1_min = beta1*x1_min + beta0
y1_max = y1_min + alpha

y2_min = beta2*x2_min + beta0
y2_max = y2_min + alpha/2

#### Sensitivity analysis with sensobol package
## Common variables
explanatory_vars = c("x1", "x2")
matrices = c("A", "B", "AB")
order = "first"
N = 2^14
type = "QRN"
first = "saltelli"
total = "jansen"

## Sobol matrices
sobol_mat = sobol_matrices(matrices = matrices, N = N, order = order, type = type, params = explanatory_vars)

# Rescale the Sobol matrices
sobol_mat[, "x1"] = qunif(p = sobol_mat[, "x1"], min = x1_min, max = x1_max)
sobol_mat[, "x2"] = qunif(p = sobol_mat[, "x2"], min = x2_min, max = x2_max)

y = beta0 + beta1*sobol_mat[, "x1"] + beta2*sobol_mat[, "x2"]

ind_data = sobol_indices(matrices = matrices, Y = y, N = N, params = explanatory_vars, first = first, total = total, order = order)
ind_data$results[sensitivity == "Si"]

#### Sensitivity analysis with sensitivity package
## Create the model
model_fct = function(X, beta0, beta1, beta2)
{
	varnames = c("x1", "x2")
	if (!all(varnames %in% colnames(X)))
		stop("The colnames of the matrix X mismatch the required names of the model")

	return (beta0 + beta1*X[, "x1"] + beta2*X[, "x2"])
}

## Generate the N x 2k matrix, with k the number of factors in the model (here, k = 2), and matrices X1 and X2
n = 2^14
k = 2
big_mat = matrix(data = runif(n = n*2*k), nrow = n, ncol = 2*k)

X1 = big_mat[, 1:2] # Correspond to matrix A in Saltelli 2008, p. 165
colnames(X1) = c("x1", "x2")

X2 = big_mat[, 3:4] # Correspond to matrix B in Saltelli 2008, p. 165
colnames(X2) = c("x1", "x2")

# Rescale the Sobol matrices
X1[, "x1"] = qunif(p = X1[, "x1"], min = x1_min, max = x1_max)
X1[, "x2"] = qunif(p = X1[, "x2"], min = x2_min, max = x2_max)

X2[, "x1"] = qunif(p = X2[, "x1"], min = x1_min, max = x1_max)
X2[, "x2"] = qunif(p = X2[, "x2"], min = x2_min, max = x2_max)

sobolmartinez(model = model_fct, X1 = X1, X2 = X2, beta0 = beta0, beta1 = beta1, beta2 = beta2)
