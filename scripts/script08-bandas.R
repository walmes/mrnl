library(lattice)
library(latticeExtra)
library(rootSolve)
library(nlme)
library(MASS)
source("http://fisher.osu.edu/~schroeder_9/AMIS900/blockdiag.R")
klib <- data.frame(k = c(51.03, 57.76, 26.6, 60.65, 87.07, 64.67, 91.28, 105.22, 
    72.74, 81.88, 97.62, 90.14, 89.88, 113.22, 90.91, 115.39, 112.63, 87.51, 104.69, 
    120.58, 114.32, 130.07, 117.65, 111.69, 128.54, 126.88, 127, 134.17, 149.66, 
    118.25, 132.67, 154.48, 129.11, 151.83, 147.66, 127.3), t = rep(c(15, 30, 45, 
    60, 75, 90, 120, 150, 180, 210, 240, 270), each = 3))
plot(k ~ t, data = klib)
expo.der <- deriv3(~A * (1 - exp(-log(2) * t/V)) + D * t, c("A", "V", "D"), function(t, 
    A, V, D) {
    NULL
})
str(expo.der)
p0 <- xyplot(k ~ t, data = klib, col = 1, xlab = "Período de incubacão (dias)", 
    ylab = "Potássio liberado acumulado (mg/kg de solo)")
p0
start <- list(A = 90, V = 20, D = 0.2)
p0 + layer(with(start, panel.curve(expo.der(x, A, V, D), add = TRUE, col = 2)))
n0 <- nls(k ~ expo.der(t, A, V, D), data = klib, start = start)
summary(n0)
confint(n0)
par(mfrow = c(2, 2))
plot(profile(n0))
layout(1)
rms.curv(n0)
str(n0)
str(n0$m$fitted())
c(n0$m$fitted())
attr(n0$m$fitted(), "gradient")
attr(n0$m$fitted(), "hessian")
f <- function(beta, t) {
    with(as.list(beta), A * (1 - exp(-log(2) * t/V)) + D * t)
}
gradient(f, x = coef(n0), t = klib$t)
pred <- data.frame(t = seq(0, 300, l = 100))
der <- do.call(expo.der, args = c(list(t = pred$t), as.list(coef(n0))))
str(der)
F <- attr(der, "gradient")
U <- chol(vcov(n0))
pred$se <- sqrt(apply(F %*% t(U), 1, function(x) sum(x^2)))
tval <- qt(p = c(lwr = 0.025, fit = 0.5, upr = 0.975), df = df.residual(n0))
me <- outer(pred$se, tval, "*")
pred <- cbind(pred, sweep(me, 1, der, "+"))
update(p0, xlim = c(0, NA), ylim = c(0, NA)) + as.layer(xyplot(fit + lwr + upr ~ 
    t, data = pred, type = "l", lty = c(1, 2, 2), col = 1))
umi <- read.table("http://www.leg.ufpr.br/~walmes/data/emr11.txt", header = TRUE, 
    sep = "\t")
str(umi)
p0 <- xyplot(umid ~ tempo | nome, data = umi, col = 1, ylab = "Umidade do solo (%)", 
    xlab = "Tempo da amostra no microondas (min)")
p0
p1 <- xyplot(umrel ~ tempo | nome, data = umi, col = 1, ylab = "Umidade relativa do solo", 
    xlab = "Tempo da amostra no microondas (min)")
p1
n0 <- nlsList(umrel ~ SSlogis(tempo, A, x0, S) | nome, data = umi)
summary(n0)
coef(n0)
pred <- expand.grid(tempo = 0:45, nome = levels(umi$nome))
pred$umrel <- predict(n0, newdata = pred)
p1 + as.layer(xyplot(umrel ~ tempo | nome, col = 2, data = pred, type = "l"))
n1 <- gnls(umrel ~ SSlogis(tempo, A, x0, S), params = list(A ~ -1 + nome, x0 ~ -1 + 
    nome, S ~ -1 + nome), data = umi, start = unlist(coef(n0)))
summary(n1)
t.grid <- seq(0, 50, by = 0.5)
pred <- expand.grid(tempo = t.grid, nome = levels(umi$nome))
pred$umrel <- predict(n1, newdata = pred)
model <- deriv3(~A/(1 + exp(-(x - x0)/S)), c("A", "x0", "S"), function(x, A, x0, 
    S) NULL)
coef(n0)
F <- lapply(split(coef(n0), f = 1:4), function(i) {
    attr(do.call(model, c(list(x = t.grid), c(i))), "gradient")
})
str(F)
bF <- blockdiag(F)
str(bF)
head(bF)
colnames(vcov(n1))
i <- ncol(bF)
j <- ncol(F[[1]])
bF <- bF[, rep(seq(1, i, j), t = j) + rep(0:(j - 1), e = i/j)]
head(bF)
U <- chol(vcov(n1))
pred$se <- sqrt(apply(bF %*% t(U), 1, function(x) sum(x^2)))
tval <- qt(p = c(lwr = 0.025, fit = 0.5, upr = 0.975), df = nrow(umi) - ncol(bF))
me <- outer(pred$se, tval, "*")
pred <- cbind(pred, sweep(me, 1, pred$umrel, "+"))
str(pred)
update(p1, xlim = c(0, 50), panel = function(...) {
    panel.abline(v = seq(0, 50, 5), h = seq(0, 1, 0.1), lty = 3, col = "gray70")
    panel.abline(h = 1, v = 40, lty = 3, lwd = 2)
    panel.xyplot(...)
}) + as.layer(xyplot(fit + lwr + upr ~ tempo | nome, data = pred, type = "l", col = 1, 
    lty = c(1, 2, 2), lwd = c(1.5, 1, 1)))
