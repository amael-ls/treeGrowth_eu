
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(terra)

f1 = function(x)
	return (dnorm(x, mean = 0, sd = 1))

f2 = function(x)
	return (dnorm(x, mean = 1, sd = 2))

curve(f1, from = -8, to = 8, lwd = 2, col = "#CD212A", ylim = c(0, 0.4), xlim = c(-7, 9))
curve(f2, lwd = 2, col = "#A212DC", add = TRUE)

x = runif(n = 5e3, min = -30, max = 30)
x = sort(x)
x1 = f1(x)
x2 = f2(x)

tt = cbind(id = 1, part = 1, x, x1)
p1 = vect(tt, type = "polygons")

tt = cbind(id = 1, part = 1, x, x2)
p2 = vect(tt, type = "polygons")

terra::is.valid(p1)
terra::is.valid(p2)

big = union(p1, p2)
big = aggregate(big)
inter = intersect(p1, p2)

plot(p1, asp = NA)
plot(p2, asp = NA, add = TRUE)

plot(big, asp = NA)
plot(inter, asp = NA)

expanse(inter)
expanse(inter)/expanse(big)

#### Analytical solution: https://www.tandfonline.com/doi/abs/10.1080/03610928908830127
interGauss = function(mu1, mu2, sd1, sd2, is_sd)
{
	if (!is_sd)
		stop("Please use sd, not variance!")

	if (sd1 == sd2)
		stop("Not programmed when equal variances")

	v1 = sd1^2
	v2 = sd2^2

	inter_minus = (mu1*v2 - mu2*v1 - sd1*sd2*sqrt((mu1 - mu2)^2 + 2*(v2 - v1)*log(sd2/sd1)))/(v2 - v1)
	inter_plus = (mu1*v2 - mu2*v1 + sd1*sd2*sqrt((mu1 - mu2)^2 + 2*(v2 - v1)*log(sd2/sd1)))/(v2 - v1)

	inter1 = min(inter_minus, inter_plus)
	inter2 = max(inter_minus, inter_plus)

	int_f = function(x, mu1, mu2, sd1, sd2)
	{
		f1 = dnorm(x, mean = mu1, sd = sd1)
		f2 = dnorm(x, mean = mu2, sd = sd2)
		return (pmin(f1, f2))
	}
	
	area = integrate(int_f, -Inf, Inf, mu1, mu2, sd1, sd2)
	
	return (list(area = area, inter1 = inter1, inter2 = inter2))
}

aa = interGauss(mu1 = 0, mu2 = 1, sd1 = 1, sd2 = 2, is_sd = TRUE)

curve(f1, from = -8, to = 8, lwd = 2, col = "#CD212A", ylim = c(0, 0.4), xlim = c(-7, 9))
curve(f2, lwd = 2, col = "#A212DC", add = TRUE)
points(x = aa[["inter1"]], y = f2(aa[["inter1"]]), pch = 19)
points(x = aa[["inter2"]], y = f2(aa[["inter2"]]), pch = 19)
