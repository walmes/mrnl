library(lattice)
library(latticeExtra)
library(MASS)
library(rpanel)
library(nlstools)
source("../funções/viciobox.R")
da <- data.frame(x = 0:10, y = c(4.67, 6.2, 8.43, 8.95, 9.98, 9.82, 10, 9.07, 7.68, 
    5.92, 3.59))
p0 <- xyplot(y ~ x, da, col = 1)
p0
quado <- deriv3(~b0 + b1 * x + b2 * x^2, c("b0", "b1", "b2"), function(x, b0, b1, 
    b2) {
    NULL
})
quadc <- deriv3(~y0 + b2 * (x - x0)^2, c("y0", "x0", "b2"), function(x, y0, x0, b2) {
    NULL
})
no <- nls(y ~ quado(x, b0, b1, b2), data = da, start = list(b0 = 1, b1 = 1, b2 = 1))
summary(no)
nc <- nls(y ~ quadc(x, y0, x0, b2), data = da, start = list(y0 = 1, x0 = 1, b2 = 1))
summary(nc)
deviance(no)
deviance(nc)
rms.curv(no)
rms.curv(nc)
str(no)
attr(no$m$fitted(), "gradient")
attr(no$m$fitted(), "hessian")
attr(nc$m$fitted(), "gradient")
attr(nc$m$fitted(), "hessian")
p0 + layer(with(as.list(coef(nc)), panel.curve(y0 + b2 * (x - x0)^2, add = TRUE, 
    col = 2)))
da <- data.frame(y = c(51.03, 57.76, 26.6, 60.65, 87.07, 64.67, 91.28, 105.22, 72.74, 
    81.88, 97.62, 90.14, 89.88, 113.22, 90.91, 115.39, 112.63, 87.51, 104.69, 120.58, 
    114.32, 130.07, 117.65, 111.69, 128.54, 126.88, 127, 134.17, 149.66, 118.25, 
    132.67, 154.48, 129.11, 151.83, 147.66, 127.3), x = rep(c(15, 30, 45, 60, 75, 
    90, 120, 150, 180, 210, 240, 270), each = 3), r = rep(1:3, 12))
p0 <- xyplot(y ~ x, data = da, col = 1, xlim = c(0, NA), ylim = c(0, NA), xlab = "Dias após a inclubação", 
    ylab = "Conteúdo liberado acumulado de potássio")
p0
m1 <- deriv3(~A * (1 - exp(-C * x)), c("A", "C"), function(x, A, C) {
    NULL
})
m2 <- deriv3(~A * (1 - 2^(-x/V)), c("A", "V"), function(x, A, V) {
    NULL
})
m3 <- deriv3(~A * x/(V + x), c("A", "V"), function(x, A, V) {
    NULL
})
n1 <- nls(y ~ m1(x, A, C), data = da, start = list(A = 140, C = 0.1))
n2 <- nls(y ~ m2(x, A, V), data = da, start = list(A = 140, V = 40))
n3 <- nls(y ~ m3(x, A, V), data = da, start = list(A = 140, V = 40))
nn <- list(n1 = n1, n2 = n2, n3 = n3)
sapply(nn, deviance)
sapply(nn, coef)
lapply(nn, function(i) cov2cor(vcov(i)))
lapply(nn, rms.curv)
l1 <- function(A, C, y, x) {
    m <- A * (1 - exp(-C * x))
    s <- sqrt(sum((y - m)^2)/length(y))
    sum(dnorm(y, m, s, log = TRUE))
}
L1 <- Vectorize(l1, c("A", "C"))
l2 <- function(A, V, y, x) {
    m <- A * (1 - 2^(-x/V))
    s <- sqrt(sum((y - m)^2)/length(y))
    sum(dnorm(y, m, s, log = TRUE))
}
L2 <- Vectorize(l2, c("A", "V"))
l3 <- function(A, V, y, x) {
    m <- A * x/(V + x)
    s <- sqrt(sum((y - m)^2)/length(y))
    sum(dnorm(y, m, s, log = TRUE))
}
L3 <- Vectorize(l3, c("A", "V"))
lapply(nn, confint.default)
A.grid <- seq(120, 180, l = 100)
V.grid <- seq(20, 60, l = 100)
C.grid <- seq(0.013, 0.032, l = 100)
qchi <- qchisq(c(0.99, 0.95, seq(0.9, 0.2, -0.1)), df = 2)
surf1 <- outer(A.grid, C.grid, L1, y = da$y, x = da$x)
max(surf1)
logLik(n1)
surf1 <- -2 * (surf1 - logLik(n1))
surf2 <- outer(A.grid, V.grid, L2, y = da$y, x = da$x)
max(surf2)
logLik(n2)
surf2 <- -2 * (surf2 - logLik(n2))
surf3 <- outer(A.grid, V.grid, L3, y = da$y, x = da$x)
max(surf3)
logLik(n3)
surf3 <- -2 * (surf3 - logLik(n3))
par(mfcol = c(2, 2))
layout(matrix(c(1, 2, 0, 3), 2, 2))
contour(A.grid, C.grid, surf1, levels = qchi, xlab = "A", ylab = "C", col = "gray70", 
    main = "Superfície da função de log-verossimilhaça (m1)")
abline(v = coef(n1)["A"], h = coef(n1)["C"], lty = 2)
contour(A.grid, V.grid, surf2, levels = qchi, xlab = "A", ylab = "V", col = "gray70", 
    main = "Superfície da função de log-verossimilhaça (m2)")
abline(v = coef(n2)["A"], h = coef(n2)["V"], lty = 2)
contour(A.grid, V.grid, surf3, levels = qchi, xlab = "A", ylab = "V", col = "gray70", 
    main = "Superfície da função de log-verossimilhaça (m3)")
abline(v = coef(n3)["A"], h = coef(n3)["V"], lty = 2)
layout(1)
lapply(nn, function(i) cov2cor(vcov(i)))
lapply(nn, rms.curv)
lapply(nn, biasbox)
dados <- structure(list(dia = c(0L, 0L, 0L, 0L, 0L, 0L, 14L, 14L, 14L, 14L, 14L, 
    14L, 24L, 24L, 24L, 24L, 24L, 24L, 38L, 38L, 38L, 38L, 38L, 38L, 54L, 54L, 54L, 
    54L, 54L, 54L, 69L, 69L, 69L, 69L, 69L, 69L, 81L, 81L, 81L, 81L, 81L, 81L, 92L, 
    92L, 92L, 92L, 92L, 92L, 106L, 106L, 106L, 106L, 106L, 106L, 117L, 117L, 117L, 
    117L, 117L, 117L, 128L, 128L, 128L, 128L, 128L, 128L, 142L, 142L, 142L, 154L, 
    154L, 166L, 166L), inc2 = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.09, 
    0, 0.13, 0, 0, 0.08, 0.18, 0.05, 0.38, 0.25, 0.13, 0.15, 0.5, 0.14, 0.33, 0.4, 
    0.63, 0.16, 0.5, 0.62, 0.56, 0.73, 0.63, 0.3, 0.89, 0.57, 0.92, 0.62, 0.75, 0.52, 
    1, 0.67, 1, 1, 1, 0.53, 1, 0.85, 1, 1, 1, 0.93, 1, 1, 1, 1, 1, 0.96, 1, 1, 1, 
    1, 1, 1, 1, 1, 1)), .Names = c("dia", "inc2"), class = "data.frame")
p0 <- xyplot(inc2 ~ dia, data = dados, type = c("p", "a"), jitter.x = TRUE, col = 1, 
    xlab = "Dias após aparecimento da doença", ylab = "Incidência da doença")
p0
theta.a <- c(a1 = 73.1268069921413, a2 = 15.2365145641768)
theta.b <- c(b1 = 121.447513236392, b2 = -0.0656323168643211)
theta.c <- c(c1 = 4.79945886958956, c2 = -0.0656320017487536)
theta.d <- c(d1 = 0.00816670903740312, d2 = 0.0656324092024781)
f.a <- function(x, a1, a2) {
    1/(1 + exp((a1 - x)/a2))
}
f.b <- function(x, b1, b2) {
    1/(1 + b1 * exp(b2 * x))
}
f.c <- function(x, c1, c2) {
    1/(1 + exp(c1 + c2 * x))
}
f.d <- function(x, d1, d2) {
    1/(1 + (-1 + 1/d1) * exp(-d2 * x))
}
f.apanel <- function(p) {
    plot(inc2 ~ dia, dados)
    curve(f.a(x, p$a1, p$a2), add = TRUE)
    abline(v = p$a1 + c(-1, 0, 1) * p$a2 * 10/2, lty = 2)
    p
}
panel <- rp.control("modelo 1 - R")
rp.slider(panel, a1, 20, 160, initval = 80, showvalue = TRUE, action = f.apanel)
rp.slider(panel, a2, 0.1, 50, initval = 5, showvalue = TRUE, action = f.apanel)
f.bpanel <- function(p) {
    plot(inc2 ~ dia, dados)
    curve(f.b(x, p$b1, p$b2), add = TRUE)
    p
}
panel <- rp.control("modelo 2")
rp.slider(panel, b1, 0, 200, initval = 120, showvalue = TRUE, action = f.bpanel)
rp.slider(panel, b2, 0, -3, initval = -0.06, showvalue = TRUE, action = f.bpanel)
f.cpanel <- function(p) {
    plot(inc2 ~ dia, dados)
    curve(f.c(x, p$c1, p$c2), add = TRUE)
    p
}
panel <- rp.control("modelo 3")
rp.slider(panel, c1, 0, 20, initval = 5, showvalue = TRUE, action = f.cpanel)
rp.slider(panel, c2, 0, -3, initval = -0.06, showvalue = TRUE, action = f.cpanel)
f.dpanel <- function(p) {
    plot(inc2 ~ dia, dados)
    curve(f.d(x, p$d1, p$d2), add = TRUE)
    p
}
panel <- rp.control("modelo 4")
rp.slider(panel, d1, 0, 0.016, initval = 0.008, showvalue = TRUE, action = f.dpanel)
rp.slider(panel, d2, 0, 3, initval = 0.06, showvalue = TRUE, action = f.dpanel)
ll <- function(th, y, x, model, C = 1) {
    ex <- do.call(model, list(x = x, th = th))
    sd <- sqrt(crossprod(y - ex)/length(x))
    ll <- sum(dnorm(y, mean = ex, sd = sd, log = TRUE))
    ll * C
}
f.a <- function(x, th) {
    1/(1 + exp((th[1] - x)/th[2]))
}
f.b <- function(x, th) {
    1/(1 + th[1] * exp(th[2] * x))
}
f.c <- function(x, th) {
    1/(1 + exp(th[1] + th[2] * x))
}
f.d <- function(x, th) {
    1/(1 + (-1 + 1/th[1]) * exp(-th[2] * x))
}
y <- dados$inc
x <- dados$dia
init.list <- list(A = list(par = c(80, 13), model = f.a), B = list(par = c(120, -0.06), 
    model = f.b), C = list(par = c(5, -0.06), model = f.c), D = list(par = c(0.008, 
    0.065), model = f.d))
fixed.list <- list(fn = ll, x = x, y = y, method = "BFGS", control = list(fnscale = -1))
op.all <- lapply(init.list, function(i) {
    op <- do.call(optim, c(i, fixed.list))
    op
})
pars <- sapply(op.all, "[[", "par")
pars
ll0 <- sapply(op.all, "[[", "value")
ll0
update(p0, type = "p") + layer(panel.curve(f.a(x, th = pars[, 1]), 0, 180, add = TRUE, 
    lwd = 1.4))
leng <- 100
grid.a <- expand.grid(th1 = seq(67, 79, l = leng), th2 = seq(11, 20, l = leng))
grid.a$ll <- apply(grid.a, 1, ll, y = y, x = x, model = f.a)
grid.a$dev <- with(grid.a, -2 * (ll - ll0["A"]))
grid.b <- expand.grid(th1 = seq(30, 700, l = leng), th2 = seq(-0.09, -0.05, l = leng))
grid.b$ll <- apply(grid.b, 1, ll, y = y, x = x, model = f.b)
grid.b$dev <- with(grid.b, -2 * (ll - ll0["B"]))
grid.c <- expand.grid(th1 = seq(3.5, 6.7, l = leng), th2 = seq(-0.09, -0.05, l = leng))
grid.c$ll <- apply(grid.c, 1, ll, y = y, x = x, model = f.c)
grid.c$dev <- with(grid.c, -2 * (ll - ll0["C"]))
grid.d <- expand.grid(th1 = seq(5e-04, 0.027, l = leng), th2 = seq(0.049, 0.09, l = leng))
grid.d$ll <- apply(grid.d, 1, ll, y = y, x = x, model = f.d)
grid.d$dev <- with(grid.d, -2 * (ll - ll0["D"]))
grid.full <- rbind(grid.a, grid.b, grid.c, grid.d)
grid.full$para <- rep(c("A", "B", "C", "D"), each = leng^2)
expr <- c(expression(f[a](x) == 1/(1 + exp((a[1] - x)/a[2]))), expression(f[b](x) == 
    1/(1 + b[1] * exp(b[2] * x))), expression(f[c](x) == 1/(1 + exp(c[1] + c[2] * 
    x))), expression(f[d](x) == 1/(1 + (-1 + 1/d[1]) * exp(-d[2] * x))))
levels <- c(0.05, 0.5, 0.75, 0.9, 0.95, 0.99)
qchi <- qchisq(levels, df = 2)
contourplot(dev ~ th1 + th2 | para, data = grid.full, at = qchi, scales = "free", 
    labels = list(labels = as.character(levels), cex = 0.8), xlab = expression("(a,b,c,d)"[1]), 
    ylab = expression("(a,b,c,d)"[2]), strip = strip.custom(bg = "gray90", factor.levels = (expr)), 
    par.settings = list(layout.heights = list(strip = 2)), panel = function(...) {
        panel.contourplot(...)
        i <- which.packet()
        panel.abline(v = pars[1, i], h = pars[2, i], lty = 2)
    })
f.a <- deriv3(~1/(1 + exp((a1 - x)/a2)), c("a1", "a2"), function(x, a1, a2) NULL)
f.b <- deriv3(~1/(1 + b1 * exp(b2 * x)), c("b1", "b2"), function(x, b1, b2) NULL)
f.c <- deriv3(~1/(1 + exp(c1 + c2 * x)), c("c1", "c2"), function(x, c1, c2) NULL)
f.d <- deriv3(~1/(1 + (-1 + 1/d1) * exp(-d2 * x)), c("d1", "d2"), function(x, d1, 
    d2) NULL)
f.list <- c(inc2 ~ f.a(dia, a1, a2), inc2 ~ f.b(dia, b1, b2), inc2 ~ f.c(dia, c1, 
    c2), inc2 ~ f.d(dia, d1, d2))
i.list <- list(theta.a, theta.b, theta.c, theta.d)
n00 <- lapply(seq_along(f.list), function(i) {
    nls(f.list[[i]], start = i.list[[i]], data = dados)
})
lapply(n00, rms.curv)
lapply(n00, function(i) cov2cor(vcov(i)))
lapply(n00, biasbox)
