#-----------------------------------------------------------------------
#                                            Prof. Dr. Walmes M. Zeviani
#                                leg.ufpr.br/~walmes · github.com/walmes
#                                        walmes@ufpr.br · @walmeszeviani
#                      Laboratory of Statistics and Geoinformation (LEG)
#                Department of Statistics · Federal University of Paraná
#                                       2018-nov-22 · Curitiba/PR/Brazil
#-----------------------------------------------------------------------

#=======================================================================
# Análise com estrutura experimental e efeito de blocos.

#-----------------------------------------------------------------------
# Pacotes.

library(lattice)
library(latticeExtra)
library(nlme)
library(rootSolve)
library(rpanel)
library(wzRfun)

#-----------------------------------------------------------------------

link <- "http://www.leg.ufpr.br/~walmes/data/soja.txt"
da <- read.table(link, header = TRUE, sep = "\t", dec = ",")
da <- subset(da, select = c(agua, potassio, rengrao, bloco))
da <- da[-74, ]
da$AG <- factor(da$agua)
names(da)[c(2:3)] <- c("x", "y")

db <- groupedData(y ~ x | bloco, data = da, order = FALSE)
str(db)

p0 <- xyplot(y ~ x | AG,
             data = db,
             layout = c(3, 1),
             col = 1,
             xlab = expression("Conteúdo de potássio no solo" ~ (mg ~ dm^{-3})),
             ylab = expression("Rendimento de grãos" ~ (g ~ vaso^{-1})))
print(p0)

#-----------------------------------------------------------------------
# Alguns modelos segmentados.

fx0 <- function(panel) {
    with(panel, {
        xb <- -b1/(2 * b2)
        curve(b0 + (b1 * x + b2 * x^2) * (x <= xb) + (b1 * xb + b2 * xb^2) * (x > xb),
              xlim[1], xlim[2], ylim = ylim)
        abline(a = b0, b = b1, lty = 2)
        abline(v = xb, lty = 2)
    })
    panel
}
panel <- rp.control(title = "Original", xlim = c(0, 10), ylim = c(0, 10))
rp.slider(panel, b0, 0, 10, initval = 0, showvalue = TRUE, action = fx0)
rp.slider(panel, b1, 0, 15, initval = 2, showvalue = TRUE, action = fx0)
rp.slider(panel, b2, -5, 0, initval = -0.2, showvalue = TRUE, action = fx0)

fx1 <- function(panel) {
    with(panel, {
        curve(yb + (x <= xb) * b2 * (x - xb)^2, xlim[1], xlim[2], ylim = ylim)
        abline(v = xb - 0:1, h = yb - c(0, -b2), lty = 2)
    })
    panel
}
panel <- rp.control(title = "Canônica", xlim = c(0, 10), ylim = c(0, 10))
rp.slider(panel, yb, 0, 10, initval = 8, showvalue = TRUE, action = fx1)
rp.slider(panel, xb, 0, 10, initval = 5, showvalue = TRUE, action = fx1)
rp.slider(panel, b2, -5, 0, initval = -0.2, showvalue = TRUE, action = fx1)

fx2 <- function(panel) {
    with(panel, {
        curve(b0 + (x < xb) * (b1 * x - 0.5 * b1 * x^2/xb) + (x >= xb) * (0.5 * b1 *
            xb), xlim[1], xlim[2], ylim = ylim)
        abline(v = xb - 0:1, lty = 2)
        abline(a = b0, b = b1, lty = 2)
    })
    panel
}
panel <- rp.control(title = "Mistura", xlim = c(0, 10), ylim = c(0, 10))
rp.slider(panel, b0, 0, 10, initval = 0, showvalue = TRUE, action = fx2)
rp.slider(panel, b1, 0, 10, initval = 3, showvalue = TRUE, action = fx2)
rp.slider(panel, xb, 0, 10, initval = 5, showvalue = TRUE, action = fx2)

#-----------------------------------------------------------------------

n1.lps <- nlme(y ~ th0 + th1 * x * (x <= thb) + th1 * thb * (x > thb),
               data = db,
               method = "ML",
               fixed = th0 + th1 + thb ~ AG,
               random = th0 ~ 1 | bloco,
               start = c(15, 0, 0, 0.22, 0, 0, 40, 20, 50))

# summary(n1.lps)
anova(n1.lps, type = "marginal")

n1.lpr <- update(n1.lps,
                 fixed = list(th0 + th1 ~ 1, thb ~ AG),
                 random = th0 ~ 1 | bloco,
                 start = c(15, 0.22, 40, 20, 50))
anova(n1.lps, n1.lpr)

pred <- expand.grid(x = 0:180, AG = levels(da$AG))
pred$y1s <- predict(n1.lps, newdata = pred, level = 0)
pred$y1r <- predict(n1.lpr, newdata = pred, level = 0)

print(p0) +
    as.layer(xyplot(y1s + y1r ~ x | AG, data = pred, type = "l"))

#-----------------------------------------------------------------------
# Para colocar as bandas de confiança.

# Tem parâmetros comuns aos níveis de água (AG).
numF1 <- function(theta, xi, AG) {
    theta[1] +
        (xi <= AG %*% theta[3:5]) * (theta[2] * xi) +
        (xi > AG %*% theta[3:5]) * (theta[2] * AG %*% theta[3:5])
}

pred1 <- expand.grid(x = seq(0, 180, 2), AG = levels(da$AG))
c1 <- fixef(n1.lpr)
F1 <- gradient(numF1,
               x = c1,
               xi = pred1$x,
               AG = model.matrix(~AG, pred1))

dim(F1)
colnames(F1) == names(c1)
head(F1)

U1 <- chol(vcov(n1.lpr))
pred1$se <- sqrt(apply(F1 %*% t(U1), 1, function(x) sum(x^2)))
zval <- qnorm(p = c(lwr = 0.025, fit = 0.5, upr = 0.975))
me <- outer(pred1$se, zval, "*")
fit <- predict(n1.lpr, newdata = pred1, level = 0)
pred1 <- cbind(pred1, sweep(me, 1, fit, "+"))
str(pred1)

print(p0) +
    as.layer(xyplot(fit ~ x | AG,
                    data = pred1,
                    type = "l",
                    col = 2,
                    prepanel = prepanel.cbH,
                    cty = "bands",
                    ly = pred1$lwr,
                    uy = pred1$upr,
                    panel = panel.cbH))

#-----------------------------------------------------------------------

n2.qos <- nlme(y ~ b0 + (b1 * x + b2 * x^2) * (x <= (-0.5 * b1/b2)) +
                   (b1 * (-0.5 * b1/b2) +
                    b2 * (-0.5 * b1/b2)^2) * (x > (-0.5 * b1/b2)),
               data = db,
               method = "ML",
               fixed = b0 + b1 + b2 ~ AG,
               random = b0 ~ 1 | bloco,
               start = c(15, 0, 0, 0.22, 0, 0, -0.001, 0, 0))

anova(n2.qos, type = "marginal")

n2.qoi <- update(n2.qos,
                 fixed = list(b0 ~ 1, b1 + b2 ~ AG),
                 random = b0 ~ 1 | bloco,
                 start = c(15, 0.22, 0, 0, -0.001, 0, 0))

anova(n2.qos, type = "marginal")
anova(n2.qos, n2.qoi)

anova(n2.qoi, type = "marginal")

n2.qor <- update(n2.qos,
                 fixed = list(b0 + b1 ~ 1, b2 ~ AG),
                 random = b0 ~ 1 | bloco,
                 start = c(15, 0.22, -0.001, 0, 0))

anova(n2.qoi, type = "marginal")
anova(n2.qos, n2.qoi, n2.qor)
anova(n2.qos, n2.qor)

pred$y2s <- predict(n2.qos, newdata = pred, level = 0)
pred$y2r <- predict(n2.qor, newdata = pred, level = 0)

print(p0) +
    as.layer(xyplot(y1s + y1r ~ x | AG,
                    data = pred,
                    type = "l")) +
    as.layer(xyplot(y2s + y2r ~ x | AG, data = pred, type = "l"))

n3.qms <- nlme(y ~ b0 +
                   (x < xb) * (b1 * x - 0.5 * b1 * x^2/xb) +
                   (x >= xb) * (0.5 * b1 * xb),
               data = db,
               method = "ML",
               fixed = b0 + b1 + xb ~ AG,
               random = b0 ~ 1 | bloco,
               start = c(15, 0, 0, 0.22, 0, 0, 88, 20, 40))

anova(n3.qms, type = "marginal")

n3.qmi <- update(n3.qms,
                 fixed = list(b0 ~ 1, b1 ~ xb ~ AG),
                 random = b0 ~ 1 | bloco,
                 start = c(15, 0.22, 0, 0, 88, 20, 40))

anova(n3.qms, type = "marginal")
anova(n3.qms, n3.qmi)

n3.qmr <- update(n3.qms,
                 fixed = list(b0 + b1 ~ 1, xb ~ AG),
                 random = b0 ~ 1 | bloco,
                 start = c(15, -0.001, 88, 20, 40))

anova(n3.qms, n3.qmr)

pred$y3s <- predict(n3.qms, newdata = pred, level = 0)
pred$y3r <- predict(n3.qmr, newdata = pred, level = 0)

print(p0) +
    as.layer(xyplot(y2s + y2r ~ x | AG, data = pred, type = "l")) +
    as.layer(xyplot(y3s + y3r ~ x | AG, data = pred, type = "l"))

logLik(n1.lpr)
logLik(n2.qor)
logLik(n3.qmr)

AIC(n1.lpr)
AIC(n2.qor)
AIC(n3.qmr)

n3.qmc <- update(n3.qms,
                 fixed = list(b0 + b1 ~ 1, xb ~ AG - 1),
                 random = b0 ~ 1 | bloco,
                 start = c(15, -0.001, 88, 90, 140))

summary(n3.qmc)$tTable
intervals(n3.qmc)

# Derivada analítica.
anaF <- function(x, b0, b1, xb) {
    cbind(
        1,
        (x <= xb) * (x - 0.5 * x^2/xb) + (x > xb) * (0.5 * xb),
        (x <= xb) * (0.5 * b1 * x^2 * xb^(-2)) + (x > xb) * (0.5 * b1))
}

# Para produzir derivada numérica.
numF <- function(theta, xi, AG) {
    theta[1] + (xi < AG %*% theta[3:5]) * (theta[2] * xi - 0.5 * theta[2] * xi^2/AG %*%
        theta[3:5]) + (xi >= AG %*% theta[3:5]) * (0.5 * theta[2] * AG %*% theta[3:5])
}

c3 <- fixef(n3.qmc)
s <- seq(0, 180, 30)

anaF(x = s, b0 = c3["b0"], b1 = c3["b1"], xb = c3[3])
gradient(numF, x = c3, xi = s, AG = cbind(rep(1, length(s)), 0, 0))[, 1:3]

pred3 <- expand.grid(x = seq(0, 180, 2), AG = levels(da$AG))
F3 <- gradient(numF,
               x = c3,
               xi = pred3$x,
               AG = model.matrix(~-1 + AG, pred3))

dim(F3)
colnames(F3) == names(c3)
head(F3)

U3 <- chol(vcov(n3.qmc))
pred3$se <- sqrt(apply(F3 %*% t(U3), 1, function(x) sum(x^2)))
zval <- qnorm(p = c(lwr = 0.025, fit = 0.5, upr = 0.975))
me <- outer(pred3$se, zval, "*")
fit <- predict(n3.qmc, newdata = pred3, level = 0)
pred3 <- cbind(pred3, sweep(me, 1, fit, "+"))
str(pred3)

print(p0) +
    as.layer(xyplot(fit ~ x | AG,
                    data = pred3,
                    type = "l",
                    col = 2,
                    prepanel = prepanel.cbH,
                    cty = "bands",
                    ly = pred3$lwr,
                    uy = pred3$upr,
                    panel = panel.cbH))

update(p0, strip = strip.custom(bg = "gray90")) +
    as.layer(xyplot(fit ~ x | AG, data = pred1, type = "l", col = 2,
                    prepanel = prepanel.cbH, cty = "bands",
                    ly = pred1$lwr, uy = pred1$upr, fill = "red",
                    panel = panel.cbH)) +
    as.layer(xyplot(fit ~ x | AG, data = pred3, type = "l", col = 4,
                    prepanel = prepanel.cbH, cty = "bands",
                    ly = pred3$lwr, uy = pred3$upr, fill = "blue",
                    panel = panel.cbH)) +
    layer(panel.abline(v = c(c1[3], c3[3]), col = c(2, 4)),
          packets = 1) +
    layer(panel.abline(v = c(sum(c1[3:4]), c3[4]), col = c(2, 4)),
          packets = 2) +
    layer(panel.abline(v = c(sum(c1[c(3, 5)]), c3[5]), col = c(2, 4)),
          packets = 3)

#-----------------------------------------------------------------------
