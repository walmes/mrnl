library(lattice)
library(latticeExtra)
library(rpanel)
library(nlme)
library(car)
library(plyr)
source("../funções/filleddensityplot.R")
cra <- read.table("http://www.leg.ufpr.br/~walmes/data/cra_manejo.txt", header = TRUE, 
    sep = "\t")
cras <- subset(cra, condi == "LVA3,5" & prof == 80)
cras$tens[cras$tens == 0] <- 0.1
cras <- with(cras, aggregate(cbind(umid), list(posi = posi, tens = tens), mean))
str(cras)
p0 <- xyplot(umid ~ log10(tens), groups = posi, data = cras, ylab = "Conteúdo de água do solo", 
    xlab = "Tensão matricial", auto.key = list(columns = 2))
p0
vg <- function(panel) {
    plot(umid ~ log10(tens), data = subset(cras, posi == panel$posi))
    with(panel, curve(tr + (ts - tr)/(1 + (a * 10^x)^n)^(1 - 1/n), add = TRUE, col = 2))
    panel
}
panel <- rp.control()
rp.slider(panel, ts, 0, 1, initval = 0.6, showvalue = TRUE, action = vg)
rp.slider(panel, tr, 0, 1, initval = 0.4, showvalue = TRUE, action = vg)
rp.slider(panel, a, 0.5, 5, initval = 1.3, showvalue = TRUE, action = vg)
rp.slider(panel, n, 1, 5, initval = 1.6, showvalue = TRUE, action = vg)
rp.radiogroup(panel, posi, vals = levels(cras$posi), action = vg)
n0 <- nls(umid ~ tr + (ts - tr)/(1 + (a * tens)^n)^(1 - 1/n), data = subset(cras, 
    posi == "EL"), start = list(tr = 0.37, ts = 0.59, a = 0.5, n = 1.47))
summary(n0)
n1 <- nls(umid ~ tr + (ts - tr)/(1 + (a * tens)^n)^(1 - 1/n), data = subset(cras, 
    posi == "L"), start = list(tr = 0.37, ts = 0.59, a = 0.5, n = 1.47))
summary(n1)
cols <- trellis.par.get()$superpose.line$col[1:2]
p0 + layer(with(as.list(coef(n0)), panel.curve(tr + (ts - tr)/(1 + (a * 10^x)^n)^(1 - 
    1/n), add = TRUE, col = cols[1]))) + layer(with(as.list(coef(n1)), panel.curve(tr + 
    (ts - tr)/(1 + (a * 10^x)^n)^(1 - 1/n), add = TRUE, col = cols[2])))
Scalc <- function(obj) {
    if (class(obj) == "nls") {
        with(as.list(coef(obj)), {
            -(ts - tr) * n * (1 + 1/(1 - 1/n))^(-(1 - 1/n) - 1)
        })
    }
    else {
        names(obj) <- c("tr", "ts", "a", "n")
        with(as.list(obj), {
            -(ts - tr) * n * (1 + 1/(1 - 1/n))^(-(1 - 1/n) - 1)
        })
    }
}
Scalc(n0)
Scalc(n1)
m <- coef(n0)
V <- vcov(n0)
round(V, 5)
boot0 <- mvrnorm(n = 1000, mu = m, Sigma = V)
str(boot0)
head(boot0)
S0 <- apply(boot0, 1, Scalc)
plot(density(S0, na.rm = TRUE))
abline(v = quantile(S0, prob = c(0.025, 0.975), na.rm = TRUE), col = 2)
m <- coef(n1)
V <- vcov(n1)
boot1 <- mvrnorm(n = 1000, mu = m, Sigma = V)
S1 <- apply(boot1, 1, Scalc)
plot(density(S1, na.rm = TRUE))
abline(v = quantile(S1, prob = c(0.025, 0.975), na.rm = TRUE), col = 4)
densityplot(~S0 + S1, n = 500)
ecdfplot(~S0 + S1, n = 500)
qqmath(~S0 + S1, n = 500)
Sdif <- S0 - S1
densityplot(~Sdif, n = 500)
quantile(Sdif, prob = c(0.025, 0.975), na.rm = TRUE)
B <- 500
res <- c(-residuals(n0), residuals(n0))
s <- summary(n0)$sigma
R <- as.data.frame(cbind(r1 = sample(res, B * length(res)/2, replace = TRUE), r2 = rnorm(B * 
    length(res)/2, 0, sd = s), tens = unique(cras$tens), f = fitted(n0)))
R$y1 <- R$f + R$r1
R$y2 <- R$f + R$r2
R$m <- rep(1:B, each = length(res)/2)
str(R)
xyplot(y1 ~ log10(tens), groups = m, data = R)
xyplot(y2 ~ log10(tens), groups = m, data = R)
xyplot(y2 + y1 ~ log10(tens), data = R, jitter.x = TRUE, key = list(columns = 1, 
    corner = c(0.05, 0.05), text = list(c("Paramétrico", "Não paramétrico")), 
    points = list(col = cols, pch = 1)))
n001 <- nlsList(y1 ~ tr + (ts - tr)/(1 + (a * tens)^n)^(1 - 1/n) | m, data = R, start = as.list(coef(n0)))
n002 <- nlsList(y2 ~ tr + (ts - tr)/(1 + (a * tens)^n)^(1 - 1/n) | m, data = R, start = as.list(coef(n0)))
boot01 <- coef(n001)
boot02 <- coef(n002)
str(boot01)
S01 <- apply(boot01, 1, Scalc)
S02 <- apply(boot02, 1, Scalc)
plot(density(S01, na.rm = TRUE), col = 2, main = NA)
abline(v = quantile(S01, prob = c(0.025, 0.975), na.rm = TRUE), col = 2, lty = 2)
lines(density(S02, na.rm = TRUE), col = 4)
abline(v = quantile(S02, prob = c(0.025, 0.975), na.rm = TRUE), col = 4, lty = 2)
lines(density(S0, na.rm = TRUE))
abline(v = quantile(S0, prob = c(0.025, 0.975), na.rm = TRUE), col = 1, lty = 1)
L <- do.call(rbind, list(metodo3 = data.frame(S = S0, M = "Dist. amos. theta"), metodo1 = data.frame(S = S02, 
    M = "Boot. paramétrico"), metodo2 = data.frame(S = S01, M = "Boot. não param.")))
str(L)
colseq <- brewer.pal(3, "Blues")[c(3:2, 1, 2:3)]
d0 <- densityplot(~S | M, data = L, layout = c(1, 3), ylab = "Densidade", xlab = "S", 
    probs = c(0.025, 0.05, 0.9, 0.975), cols = colseq, plot.points = "rug", panel = my.densityplot)
d0 + layer(panel.abline(v = Scalc(n0), col = 2))
g <- expression(-(ts - tr) * n * (1 + 1/(1 - 1/n))^(-(1 - 1/n) - 1))
pars <- c("tr", "ts", "a", "n")
lapply(pars, D, expr = g)
grad <- deriv3(~-(ts - tr) * n * (1 + 1/(1 - 1/n))^(-(1 - 1/n) - 1), c("tr", "ts", 
    "a", "n"), function(tr, ts, a, n) {
    NULL
})
gr <- do.call(grad, as.list(coef(n0)))
gr <- attr(gr, "gradient")
gr
Scalc(n0)
sqrt(gr %*% vcov(n0) %*% t(gr))
dm0 <- deltaMethod(n0, g = "-(ts-tr)*n*(1+1/(1-1/n))^(-(1-1/n)-1)")
dm0
dm1 <- deltaMethod(n1, g = "-(ts-tr)*n*(1+1/(1-1/n))^(-(1-1/n)-1)")
dm1
S0md <- with(dm0, c(Estimate = Estimate, SE = SE, LCL = Estimate - 1.96 * SE, UCL = Estimate + 
    1.96 * SE))
S0md
S1md <- with(dm1, c(Estimate = Estimate, SE = SE, LCL = Estimate - 1.96 * SE, UCL = Estimate + 
    1.96 * SE))
S1md
d0 + layer(panel.abline(v = S0md[c(1, 3, 4)], col = 2))
vgr <- function(panel) {
    plot(umid ~ log(tens), data = subset(cras, posi == panel$posi))
    with(panel, {
        curve(tr - Si * (1 + 1/(1 - 1/n))^((1 - 1/n) + 1)/(n * (1 + exp(n * (x - 
            psii))/(1 - 1/n))^(1 - 1/n)), add = TRUE, col = 2)
        abline(v = psii, col = 2, lty = 2)
    })
    panel
}
panel <- rp.control()
rp.slider(panel, Si, -0.02, -0.06, initval = -0.03, showvalue = TRUE, action = vgr)
rp.slider(panel, tr, 0, 1, initval = 0.4, showvalue = TRUE, action = vgr)
rp.slider(panel, psii, 0, 3, initval = 1, showvalue = TRUE, action = vgr)
rp.slider(panel, n, 1, 5, initval = 1.6, showvalue = TRUE, action = vgr)
rp.radiogroup(panel, posi, vals = levels(cra$posi), action = vgr)
n2 <- nls(umid ~ tr - Si * (1 + 1/(1 - 1/n))^((1 - 1/n) + 1)/(n * (1 + exp(n * (log(tens) - 
    psii))/(1 - 1/n))^(1 - 1/n)), data = subset(cras, posi == "EL"), start = list(Si = -0.03, 
    tr = 0.3, psii = 2.07, n = 1.2))
summary(n2)
n3 <- nls(umid ~ tr - Si * (1 + 1/(1 - 1/n))^((1 - 1/n) + 1)/(n * (1 + exp(n * (log(tens) - 
    psii))/(1 - 1/n))^(1 - 1/n)), data = subset(cras, posi == "L"), start = list(Si = -0.03, 
    tr = 0.3, psii = 2.07, n = 1.2))
summary(n3)
p0 + layer(with(as.list(coef(n2)), panel.curve(tr - Si * (1 + 1/(1 - 1/n))^((1 - 
    1/n) + 1)/(n * (1 + exp(n * (log(10^x) - psii))/(1 - 1/n))^(1 - 1/n)), add = TRUE, 
    col = cols[1]))) + layer(with(as.list(coef(n3)), panel.curve(tr - Si * (1 + 1/(1 - 
    1/n))^((1 - 1/n) + 1)/(n * (1 + exp(n * (log(10^x) - psii))/(1 - 1/n))^(1 - 1/n)), 
    add = TRUE, col = cols[2])))
summary(n2)
dm0
summary(n3)
dm1
confint(n2)
confint(n3)
par(mfrow = c(2, 2))
plot(profile(n1))
layout(1)
par(mfrow = c(2, 2))
plot(profile(n3))
layout(1)
start <- split(cbind(coef(n2), coef(n3)), f = 1:4)
names(start) <- names(coef(n2))
start
n00 <- nls(umid ~ tr[posi] - Si[posi] * (1 + 1/(1 - 1/n[posi]))^((1 - 1/n[posi]) + 
    1)/(n[posi] * (1 + exp(n[posi] * (log(tens) - psii[posi]))/(1 - 1/n[posi]))^(1 - 
    1/n[posi])), data = cras, start = start)
summary(n00)
confint.default(n00)
start[["Si"]] <- mean(start[["Si"]])
n01 <- nls(umid ~ tr[posi] - Si * (1 + 1/(1 - 1/n[posi]))^((1 - 1/n[posi]) + 1)/(n[posi] * 
    (1 + exp(n[posi] * (log(tens) - psii[posi]))/(1 - 1/n[posi]))^(1 - 1/n[posi])), 
    data = cras, start = start)
summary(n01)
anova(n00, n01)
pred <- expand.grid(posi = levels(cras$posi), tens = 10^seq(-1, 3.5, l = 30))
pred$y00 <- predict(n00, newdata = pred)
pred$y01 <- predict(n01, newdata = pred)
str(pred)
xyplot(umid ~ log10(tens), groups = posi, data = cras, auto.key = list(columns = 2, 
    type = "o", divide = 1, points = FALSE, lines = TRUE), xlab = "Tensão matricial do solo (log10 kPa)", 
    ylab = "Umidade do solo (g/g)") + as.layer(xyplot(y00 ~ log10(tens), groups = posi, 
    data = pred, type = "l")) + as.layer(xyplot(y01 ~ log10(tens), groups = posi, 
    data = pred, type = "l", lty = 2))
