library(lattice)
library(latticeExtra)
library(nlme)
klib <- data.frame(k = c(51.03, 57.76, 26.6, 60.65, 87.07, 64.67, 91.28, 105.22, 
    72.74, 81.88, 97.62, 90.14, 89.88, 113.22, 90.91, 115.39, 112.63, 87.51, 104.69, 
    120.58, 114.32, 130.07, 117.65, 111.69, 128.54, 126.88, 127, 134.17, 149.66, 
    118.25, 132.67, 154.48, 129.11, 151.83, 147.66, 127.3), t = rep(c(15, 30, 45, 
    60, 75, 90, 120, 150, 180, 210, 240, 270), each = 3))
klib$pote <- 1:3
p0 <- xyplot(k ~ t, data = klib, groups = pote, xlab = "Período de incubacão (dias)", 
    ylab = "Potássio liberado acumulado (mg/kg de solo)")
p0
xyplot(k ~ t, data = klib, groups = pote, type = "b", xlab = "Período de incubacão (dias)", 
    ylab = "Potássio liberado acumulado (mg/kg de solo)")
A <- 90
V <- 20
D <- 0.2
p0 + layer(panel.curve(A * (1 - exp(-log(2) * x/V)) + D * x, add = TRUE, col = 2))
start <- c(A = A, V = V, D = D)
klib <- groupedData(k ~ t | pote, data = klib)
str(klib)
n0 <- nlme(k ~ A * (1 - exp(-log(2) * t/V)) + D * t, data = klib, fixed = A + V + 
    D ~ 1, random = A ~ 1, start = start)
summary(n0)
ranef(n0)
n2 <- nlme(k ~ A * (1 - exp(-log(2) * t/V)) + D * t, data = klib, fixed = A + V + 
    D ~ 1, random = A + V ~ 1, start = start)
summary(n2)
anova(n0, n2)
n1 <- nls(k ~ A * (1 - exp(-log(2) * t/V)) + D * t, data = klib, start = start)
anova(n0, n1)
plot(augPred(n0, level = 0:1))
pred <- expand.grid(t = seq(0, 270, 3), pote = levels(klib$pote))
pred$y0 <- predict(n0, newdata = pred, level = 0)
pred$y1 <- predict(n0, newdata = pred, level = 1)
str(pred)
key <- list(text = list(c("fixo", "unidade experimental")), lines = list(col = trellis.par.get()$superpose.line$col[1:2]))
xyplot(k ~ t | pote, data = klib, key = key, xlab = "Dias após incubação no solo", 
    ylab = "Conteúdo acumulado libeardo de potássio", strip = strip.custom(bg = "gray90")) + 
    as.layer(xyplot(y0 + y1 ~ t | pote, data = pred, type = "l"))
da <- read.table("http://www.leg.ufpr.br/~walmes/data/soja.txt", header = TRUE, sep = "\t", 
    dec = ",")
str(da)
da <- da[-74, ]
da$k <- da$potassio
da$A <- factor(da$agua)
xyplot(rengrao ~ k, groups = A, data = da)
xyplot(rengrao ~ k | A, data = da, groups = bloco, type = "o")
da <- groupedData(rengrao ~ k | bloco, data = da)
nn0 <- nlme(rengrao ~ f0 + tx * k * (k < xde) + tx * xde * (k >= xde), data = da, 
    random = f0 ~ 1, fixed = f0 + tx + xde ~ A, start = c(10, 0, 0, 0.2, 0, 0, 50, 
        10, 40))
summary(nn0)
anova(nn0, type = "marginal")
nn1 <- nlme(rengrao ~ f0 + tx * k * (k < xde) + tx * xde * (k >= xde), data = da, 
    random = f0 ~ 1, fixed = list(f0 + tx ~ 1, xde ~ A), start = c(10, 0.2, 50, 10, 
        40))
summary(nn1)
nn3 <- nlme(rengrao ~ f0 + tx * k * (k < xde) + tx * xde * (k >= xde), data = da, 
    random = xde ~ 1, fixed = list(f0 + tx ~ 1, xde ~ A), start = c(10, 0.2, 50, 
        10, 40))
summary(nn3)
anova(nn1, nn3)
nn2 <- gnls(rengrao ~ f0 + tx * k * (k < xde) + tx * xde * (k >= xde), data = da, 
    params = list(f0 + tx ~ 1, xde ~ A), start = c(10, 0.2, 50, 10, 40))
anova(nn1, nn2)
pred <- expand.grid(k = seq(0, 180, 2), A = levels(da$A), bloco = levels(da$bloco))
pred$y <- predict(nn1, newdata = pred, level = 0)
pred$yi <- predict(nn1, newdata = pred, level = 1)
str(pred)
xyplot(rengrao ~ k | A, data = da, groups = bloco) + as.layer(xyplot(y ~ k | A, data = pred, 
    type = "a", lwd = 2, col = 1)) + as.layer(xyplot(yi ~ k | A, data = pred, groups = bloco, 
    type = "l", lty = 2))
