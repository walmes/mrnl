library(lattice)
library(latticeExtra)
library(MASS)
library(bbmle)
library(nlstools)
library(nlme)
da <- read.table("http://www.leg.ufpr.br/~walmes/data/algodão.txt", header = TRUE, 
    sep = "\t", encoding = "latin1")
str(da)
da$desf <- da$desf/100
estagios <- c("Vegetativo", "Botão floral", "Florescimento", "Maçã", "Capulho")
xtabs(~estag + desf, da)
p0 <- xyplot(pcapu ~ desf | estag, data = da, layout = c(5, 1), strip = strip.custom(bg = "gray90", 
    factor.levels = estagios), col = 1, xlab = "Nível de desfolha artificial", ylab = "Peso de capulhos produzidos (g)")
p0
n0 <- lapply(split(da, da$estag), nls, formula = pcapu ~ f0 - f1 * (desf + 0.02)^exp(C), 
    start = list(f0 = 30, f1 = 15, C = 0))
lapply(n0, function(x) summary(x)$coeff)
n1 <- lapply(split(da, da$estag), nls, formula = pcapu ~ f0 - f1 * (desf + 0.02)^((log(5) - 
    log(f1))/log(xde)), start = list(f0 = 30, f1 = 10, xde = 0.7))
lapply(n1, function(x) summary(x)$coeff)
X0 <- do.call(rbind, lapply(n0, function(x) {
    cbind(summary(x)$coeff, confint.default(x))
}))
X1 <- do.call(rbind, lapply(n1, function(x) {
    cbind(summary(x)$coeff, confint.default(x))
}))
X0
X1
xs <- seq(0, 1, l = 30)
pred <- expand.grid(desf = xs, estag = levels(da$estag))
pred$y <- do.call(c, lapply(n0, predict, newdata = data.frame(desf = xs)))
p0 + as.layer(xyplot(y ~ desf | estag, data = pred, type = "l", col = 1))
myfun3 <- function(V) {
    x <- cov2cor(V)
    diag(x) <- diag(V)
    x[upper.tri(x)] <- NA
    x
}
myfun3(vcov(n1[[1]]))
myfun3(vcov(n1[[3]]))
vcovs <- lapply(lapply(n1, vcov), myfun3)
vcovs
prm0 <- deriv3(~f0 - f1 * x^exp(C), c("f0", "f1", "C"), function(x, f0, f1, C) {
    NULL
})
prm1 <- deriv3(~f0 - f1 * x^((log(5) - log(f1))/log(xde)), c("f0", "f1", "xde"), 
    function(x, f0, f1, xde) {
        NULL
    })
daS <- split(da, f = da$estag)
an0 <- lapply(daS, nls, formula = pcapu ~ prm0(desf + 0.02, f0, f1, C), start = list(f0 = 30, 
    f1 = 10, C = 1))
sapply(an0, coef)
an1 <- lapply(daS, nls, formula = pcapu ~ prm1(desf + 0.02, f0, f1, xde), start = list(f0 = 30, 
    f1 = 10, xde = 0.9))
sapply(an1, coef)
curv0 <- lapply(an0, MASS::rms.curv)
curv0 <- t(sapply(lapply(curv0, "[", 1:2), unlist))
curv1 <- lapply(an1, MASS::rms.curv)
curv1 <- t(sapply(lapply(curv1, "[", 1:2), unlist))
curvs <- rbind(curv0, curv1)
plot(curvs, pch = rep(c(1, 5), e = nrow(curv0)), xlab = "Curvatura devido efeito de parametrização", 
    ylab = "Curvatura intrínseca do modelo")
segments(curv0[, 1], curv0[, 2], curv1[, 1], curv1[, 2], col = ifelse(curv0[, 1] > 
    curv1[, 1], 1, 2))
legend("bottomright", legend = c("Prmtz original", "Prmtz xde"), pch = c(1, 5), bty = "n")
ll0 <- function(theta0, theta1, theta2, desf, pcapu) {
    x <- desf
    y <- pcapu
    ex <- theta0 - theta1 * (x + 0.02)^exp(theta2)
    sd <- sqrt(crossprod(y - ex)/(length(ex)))
    ll <- sum(dnorm(y, mean = ex, sd = sd, log = TRUE))
    -ll
}
ll1 <- function(theta0, theta1, varthetaq, desf, pcapu) {
    x <- desf
    y <- pcapu
    ex <- theta0 - theta1 * (x + 0.02)^((log(5) - log(theta1))/log(varthetaq))
    sd <- sqrt(crossprod(y - ex)/(length(ex)))
    ll <- sum(dnorm(y, mean = ex, sd = sd, log = TRUE))
    -ll
}
cl <- list(an0 = t(sapply(an0, coef)), an1 = t(sapply(an1, coef)))
colnames(cl[[1]]) <- c("theta0", "theta1", "theta2")
colnames(cl[[2]]) <- c("theta0", "theta1", "varthetaq")
al0 <- al1 <- vector(length = 5, mode = "list")
names(al0) <- names(al1) <- names(daS)
for (i in names(al0)) {
    al0[[i]] <- mle2(minuslogl = ll0, start = as.list(cl[[1]][i, ]), method = "BFGS", 
        data = c(as.list(subset(da, estag == i))))
    al1[[i]] <- mle2(minuslogl = ll1, start = as.list(cl[[2]][i, ]), method = "BFGS", 
        data = c(as.list(subset(da, estag == i))))
}
cbind(t(sapply(al0, coef)), t(sapply(al1, coef)))
cbind(sapply(al0, logLik), sapply(al1, logLik))
prof0 <- lapply(al0, function(i) plot(profile(i, maxsteps = 20), main = NULL))
paramsexpr1 <- c(expression(theta[0]), expression(theta[1]), expression(theta[2]))
paramsexpr2 <- c(expression(theta[0]), expression(theta[1]), expression(vartheta[q]))
layout(1)
plot(profile(al0[[1]], maxsteps = 20, which = 1, prof.lower = 30, prof.upper = 34))
limi <- list(c(25, 40), c(0, 25), c(-2, 4))
par(mfrow = c(5, 3), mar = c(4, 4, 2, 1))
for (estag in seq_along(al0)) {
    for (param in 1:3) {
        plot(profile(al0[[estag]], maxsteps = 20, which = param), xlab = "", ylab = "", 
            main = "", xlim = limi[[param]])
        mtext(side = 1, line = 2.5, text = paramsexpr1[param], cex = 1)
        mtext(side = 2, line = 2.5, text = "|z|", cex = 0.8)
        mtext(side = 3, line = 0.25, text = estagios[estag], cex = 1)
    }
}
layout(1)
limi <- list(c(25, 40), c(0, 25), c(0, 1))
par(mfrow = c(5, 3), mar = c(4, 4, 2, 1))
for (estag in seq_along(al0)) {
    for (param in 1:3) {
        plot(profile(al1[[estag]], maxsteps = 20, which = param), xlab = "", ylab = "", 
            main = "", xlim = limi[[param]])
        mtext(side = 1, line = 2.5, text = paramsexpr2[param], cex = 1)
        mtext(side = 2, line = 2.5, text = "|z|", cex = 0.8)
        mtext(side = 3, line = 0.25, text = estagios[estag], cex = 1)
    }
}
layout(1)
estg <- "4maça"
n0 <- nls(pcapu ~ prm0(desf + 0.01, f0, f1, C), data = subset(da, estag == estg), 
    start = list(f0 = 28, f1 = 21, C = 1.5))
n1 <- nls(pcapu ~ prm1(desf + 0.01, f0, f1, xde), data = subset(da, estag == estg), 
    start = list(f0 = 28, f1 = 21, xde = 0.63))
c0 <- nlsConfRegions(n0, exp = 3.5)
c1 <- nlsConfRegions(n1, exp = 2.7)
plot(c0)
x11()
plot(c1)
estg <- "1veg"
n0 <- nls(pcapu ~ prm0(desf + 0.01, f0, f1, C), data = subset(da, estag == estg), 
    start = list(f0 = 28, f1 = 21, C = 1.5))
n1 <- nls(pcapu ~ prm1(desf + 0.01, f0, f1, xde), data = subset(da, estag == estg), 
    start = list(f0 = 28, f1 = 21, xde = 0.63))
c0 <- nlsConfRegions(n0, exp = 3.5)
c1 <- nlsConfRegions(n1, exp = 2.7)
plot(c0)
x11()
plot(c1)
n0 <- lapply(split(da, da$estag), nls, formula = pcapu ~ f0 - f1 * (desf + 0.02)^exp(C), 
    start = list(f0 = 30, f1 = 15, C = 0))
n1 <- lapply(split(da, da$estag), nls, formula = pcapu ~ f0 - f1 * (desf + 0.02)^((log(5) - 
    log(f1))/log(xde)), start = list(f0 = 30, f1 = 10, xde = 0.7))
chute1 <- sapply(n1, coef)
chute2 <- cbind(chute1[, 1], sweep(chute1[, 2:5], 1, chute1[, 1], "-"))
c(t(chute2))
nn0 <- gnls(pcapu ~ f0 - f1 * desf^((log(5) - log(f1))/log(xde)), data = da, params = f0 + 
    f1 + xde ~ estag, start = c(t(chute2)))
summary(nn0)
anova(nn0, type = "marginal")
ctr <- cbind(1, contrasts(da$estag))
ctr <- kronecker(diag(3), ctr)
th <- ctr %*% coef(nn0)
th
vcov <- ctr %*% vcov(nn0) %*% t(ctr)
ep <- sqrt(diag(vcov))
ep
ztable <- data.frame(Estimate = th, Std.Err = ep, tvalue = th/ep, pvalue = 2 * pnorm(-abs(th/ep)))
ztable$star <- cut(ztable$pvalue, c(0, 0.01, 0.05, 0.1, 1), labels = c("**", "*", 
    ".", ""))
ztable$ics <- sweep(outer(ep, c(lwr = -1, est = 0, upr = 1) * qnorm(0.975), "*"), 
    1, th, "+")
ztable[, -(2:3)]
nn0 <- gnls(pcapu ~ f0 - f1 * desf^((log(5) - log(f1))/log(xde)), data = da, params = f0 + 
    f1 + xde ~ -1 + estag, start = c(t(chute1)))
coef(nn0)
desf.grid <- seq(0.01, 1, l = 30)
pred <- expand.grid(desf = desf.grid, estag = levels(da$estag))
pred$pcapu <- predict(nn0, newdata = pred)
ths <- split(th, f = rep(1:5, 3))
ths <- lapply(ths, function(x) {
    names(x) <- c("f0", "f1", "xde")
    x
})
model <- deriv3(~f0 - f1 * x^((log(5) - log(f1))/log(xde)), c("f0", "f1", "xde"), 
    function(x, f0, f1, xde) {
        NULL
    })
F <- lapply(ths, function(i) {
    attr(do.call(model, c(list(x = desf.grid), as.list(i))), "gradient")
})
str(F)
source("http://fisher.osu.edu/~schroeder_9/AMIS900/blockdiag.R")
bF <- blockdiag(F)
str(bF)
head(bF)
bF <- bF[, rep(seq(1, 15, 3), t = 3) + rep(0:2, e = 5)]
U <- chol(vcov(nn0))
pred$se <- sqrt(apply(bF %*% t(U), 1, function(x) sum(x^2)))
tval <- qt(p = c(lwr = 0.025, fit = 0.5, upr = 0.975), df = nrow(da) - ncol(bF))
me <- outer(pred$se, tval, "*")
pred <- cbind(pred, sweep(me, 1, pred$pcapu, "+"))
str(pred)
xyplot(pcapu ~ desf | estag, data = da, col = 1, layout = c(5, 1), xlim = c(-0.2, 
    1.2), xlab = "Desfolha artifical", ylab = "Peso de capulhos (g)", strip = strip.custom(bg = "gray90", 
    factor.levels = estagios), panel = function(...) {
    panel.xyplot(...)
    panel.abline(v = ztable$ics[10 + which.packet(), ], lty = c(2, 1, 2), col = 2)
}) + as.layer(xyplot(fit + lwr + upr ~ desf | estag, data = pred, type = "l", col = 1, 
    lty = c(1, 2, 2)))
