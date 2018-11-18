library(lattice)
library(latticeExtra)
library(nlme)
library(plyr)
source("../funções/bandas.R")
goi <- read.table("http://www.leg.ufpr.br/~walmes/data/goiaba.txt", header = TRUE, 
    sep = "\t")
str(goi)
p0 <- xyplot(peso ~ daa, data = goi, type = c("p", "a"), col = 1, xlab = "Dias após a antese", 
    ylab = "Massa fresca do fruto (g)")
p0
m <- with(goi, tapply(peso, daa, mean))
s <- with(goi, tapply(peso, daa, var))
scatter.smooth(s ~ m)
plot(s ~ m, log = "xy")
m0 <- nls(peso ~ A - (A - D) * exp(-exp(C * (daa - B))), data = goi, start = c(A = 200, 
    C = 0.09, B = 105, D = 7))
summary(m0)
m1 <- gnls(peso ~ A - (A - D) * exp(-exp(C * (daa - B))), data = goi, start = c(A = 225, 
    C = 0.05, B = 109, D = 7), weights = varPower(0.1, form = ~daa))
summary(m1)
m2 <- gnls(peso ~ A - (A - D) * exp(-exp(C * (daa - B))), data = goi, start = c(A = 225, 
    C = 0.05, B = 109, D = 7), weights = varPower(0.1, form = ~log(daa)))
summary(m2)
m3 <- gnls(peso ~ A - (A - D) * exp(-exp(C * (daa - B))), data = goi, start = c(A = 225, 
    C = 0.05, B = 109, D = 7), weights = varExp(0.03, form = ~daa))
summary(m3)
summary(m0)$coeff[, 1:2]
summary(m1)$tTable[, 1:2]
summary(m2)$tTable[, 1:2]
summary(m3)$tTable[, 1:2]
c(logLik(m0), logLik(m1), logLik(m2), logLik(m3))
anova(m1, m0)
anova(m2, m0)
anova(m3, m0)
residuos <- expand.grid(daa = goi$daa, modelo = c("Constante", "Potência", "Potência do logaritmo", 
    "Exponencial"))
residuos$respad = c(residuals(m0, type = "p"), residuals(m1, type = "p"), residuals(m2, 
    type = "p"), residuals(m3, type = "p"))
residuos$fitted = c(fitted(m0), fitted(m1), fitted(m2), fitted(m3))
qqmath(~respad | modelo, data = residuos)
xyplot(respad ~ daa | modelo, data = residuos)
qqmath(~respad | modelo, data = residuos, col = 1, xlab = "Quantis teóricos esperados da distribuição normal padrão", 
    ylab = "Resíduos padronizados", strip = strip.custom(bg = "gray90"), layout = c(2, 
        2), cex = 0.8)
xyplot(respad ~ daa | modelo, data = residuos, type = c("p", "smooth"), col = c(1), 
    xlab = "Dias após a antese", ylab = "Resíduos padronizados", strip = strip.custom(bg = "gray90"), 
    layout = c(2, 2), ylim = c(-4.5, 5.5))
model <- deriv3(~A - (A - D) * exp(-exp(C * (daa - B))), c("A", "C", "B", "D"), function(daa, 
    A, C, B, D) {
    NULL
})
l <- list(m0, m1, m2, m3)
names(l) <- levels(residuos$modelo)
L <- l
for (i in seq_along(L)) {
    m4 <- l[[i]]
    U <- chol(vcov(m4))
    pred <- data.frame(daa = seq(1, 140, 1))
    m <- model(daa = pred$daa, A = coef(m0)["A"], C = coef(m0)["C"], B = coef(m0)["B"], 
        D = coef(m0)["D"])
    F0 <- attr(m, "gradient")
    pred$y <- c(m)
    pred$se <- sqrt(apply(F0 %*% t(U), 1, function(x) sum(x^2)))
    z <- qnorm(0.975)
    pred <- transform(pred, lwr = y - z * se, upr = y + z * se)
    L[[i]] <- pred
}
str(L)
L <- ldply(L)
names(L)[1] <- "model"
L$model <- factor(L$model, levels = levels(residuos$model))
str(L)
icm <- ddply(goi, .(daa), summarise, m = mean(peso), lwr = mean(peso) - qt(0.975, 
    length(peso) - 1) * sd(peso) * sqrt(1/length(peso)), upr = mean(peso) + qt(0.975, 
    length(peso) - 1) * sd(peso) * sqrt(1/length(peso)))
ylim <- range(goi$peso)
xyplot(y ~ daa | model, data = L, type = "l", col = 1, ylab = "Massa fresca dos frutos (g)", 
    xlab = "Dias após a antese", strip = strip.custom(bg = "gray90"), prepanel = prepanel.cbH, 
    cty = "bands", ly = L$lwr, uy = L$upr, panel = panel.cbH) + layer(panel.arrows(icm$daa, 
    icm$lwr, icm$daa, icm$upr, code = 3, length = 0.05, angle = 90))
