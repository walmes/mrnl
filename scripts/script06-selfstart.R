library(lattice)
library(latticeExtra)
library(car)
library(segmented)
library(nls2)
data(wtloss, package = "MASS")
p0 <- xyplot(Weight ~ Days, data = wtloss, col = 1, xlab = "Dias em dieta", ylab = "Peso (kg)")
print(p0)
p0 + layer(with(list(f0 = 81, f1 = 102, K = 141), panel.curve(f0 + f1 * 2^(-x/K), 
    add = TRUE, col = 2)))
pe <- with(wtloss, {
    f1 <- diff(range(Weight))
    f0 <- min(Weight)
    m <- f0 + 0.5 * f1
    mm <- which.min(abs(Weight - m))[1]
    K <- Days[mm]
    list(f0 = f0, f1 = f1, K = K)
})
p0 + layer(with(pe, panel.abline(v = K, h = c(f0 + c(0, 0.5, 1) * f1), lty = 2)))
model.init <- function(x, y) {
    f1 <- diff(range(y))
    f0 <- min(y)
    m <- f0 + 0.5 * f1
    mm <- which.min(abs(y - m))[1]
    K <- x[mm]
    list(f0 = f0, f1 = f1, K = K)
}
model.init(x = wtloss$Days, y = wtloss$Weight)
model.init <- function(mCall, data, LHS) {
    xy <- data.frame(sortedXyData(mCall[["input"]], LHS, data))
    xx <- xy[["x"]]
    yy <- xy[["y"]]
    f0 <- min(yy)
    f1 <- diff(range(yy))
    m <- f0 + 0.5 * f1
    mm <- which.min(abs(yy - m))[1]
    K <- xx[mm]
    value <- c(f0, f1, K)
    names(value) <- mCall[c("f0", "f1", "K")]
    value
}
SSmodel <- selfStart(model = ~f0 + f1 * 2^(-input/K), parameters = c("f0", "f1", 
    "K"), initial = model.init, template = function(input, f0, f1, K) {
    NULL
})
SSmodel
getInitial(Weight ~ SSmodel(Days, f0, f1, K), data = wtloss)
n0 <- nls(Weight ~ SSmodel(Days, f0, f1, K), data = wtloss, trace = TRUE)
summary(n0)
cm <- read.table("http://www.leg.ufpr.br/~walmes/data/cresmicelial.txt", header = TRUE, 
    sep = "\t", colClasses = c("factor", "numeric", "numeric"))
str(cm)
p1 <- xyplot(mmdia ~ temp | isolado, data = cm, type = c("p", "a"), col = 1, ylab = "Taxa de crescimento em diâmetro do disco micelial", 
    xlab = "Temperatuda da câmara de crescimento")
print(p1)
mqc.init <- function(mCall, data, LHS) {
    xy <- data.frame(sortedXyData(mCall[["input"]], LHS, data))
    m0 <- lm(y ~ poly(x, 2, raw = TRUE), data = xy)
    c0 <- coef(m0)
    cv <- c0[3]
    xc <- -0.5 * c0[2]/cv
    fc <- c0[1] + c0[2] * xc + cv * xc
    value <- c(fc, xc, cv)
    names(value) <- mCall[c("fc", "xc", "cv")]
    value
}
SSmqc <- selfStart(model = ~fc + cv * (input - xc)^2, parameters = c("fc", "xc", 
    "cv"), initial = mqc.init, template = function(input, fc, xc, cv) {
    NULL
})
getInitial(mmdia ~ SSmqc(temp, fc, xc, cv), data = subset(cm, isolado == "1"))
cms <- split(cm, f = cm$isolado)
sapply(cms, getInitial, object = mmdia ~ SSmqc(temp, fc, xc, cv))
n0 <- lapply(cms, nls, formula = mmdia ~ SSmqc(temp, fc, xc, cv))
c0 <- t(sapply(n0, coef))
c0
ic0 <- lapply(n0, confint.default)
ic0
res <- do.call(c, lapply(n0, residuals))
qqmath(~res | cm$isolado)
c0 <- lapply(n0, coef)
xyplot(mmdia ~ temp | isolado, data = cm, col = 1, ylab = "Taxa de crescimento em diâmetro do disco micelial", 
    xlab = "Temperatuda da câmara de crescimento", panel = function(...) {
        panel.xyplot(...)
        wp <- which.packet()
        with(as.list(c0[[wp]]), panel.curve(fc + cv * (x - xc)^2, col = 2, lwd = 2))
        panel.abline(v = c0[[wp]]["xc"], lty = 2, col = 2)
        panel.abline(v = ic0[[wp]]["xc", ], lty = 3, col = 2)
    })
ic0
ics <- t(sapply(ic0, function(i) i[2, ]))
ics <- cbind(i = rownames(ics), data.frame(ics))
names(ics) <- c("i", "lwr", "upr")
str(ics)
segplot(i ~ lwr + upr, data = ics, col = 1, xlab = "Temperatura ótima para o crescimento", 
    ylab = "Isolado", draw.bands = FALSE, segments.fun = panel.arrows, ends = "both", 
    angle = 90, length = 1, unit = "mm")
m0 <- lm(mmdia ~ temp + I(temp^2), cms[[1]])
n0 <- nls(mmdia ~ SSmqc(temp, fc, xc, cv), cms[[1]])
deltaMethod(m0, "-0.5*b1/b2", parameterNames = c("b0", "b1", "b2"))
summary(n0)$coeff["xc", 1:2]
apropos("^SS", where = TRUE)
soja <- read.table("http://www.leg.ufpr.br/~walmes/data/soja.txt", header = TRUE, 
    sep = "\t", dec = ",")
str(soja)
soja <- soja[-74, ]
xyplot(rengrao ~ potassio | factor(agua), data = soja, col = 1, xlab = "Conteúdo de potássio no solo", 
    ylab = "Produção de soja por vaso")
xyplot(rengrao ~ potassio, groups = agua, data = soja, type = c("p", "a"), xlab = "Conteúdo de potássio no solo", 
    ylab = "Produção de soja por vaso")
mrp.init <- function(mCall, data, LHS) {
    xy <- data.frame(sortedXyData(mCall[["input"]], LHS, data))
    m0 <- lm(y ~ x, data = xy)
    xi <- mean(xy[["x"]])
    s0 <- segmented(m0, seg.Z = ~x, psi = list(x = c(xi)), control = seg.control(display = FALSE))
    xde <- s0$psi[1, 2]
    f0 <- coef(s0)[1]
    tx <- coef(s0)[2]
    value <- c(f0, tx, xde)
    names(value) <- mCall[c("f0", "tx", "xde")]
    value
}
SSmrp <- selfStart(model = ~f0 + tx * input + xde, parameters = c("f0", "tx", "xde"), 
    initial = mrp.init, template = function(input, f0, tx, xde) {
        NULL
    })
body(SSmrp)
getInitial(rengrao ~ SSmrp(potassio, f0, tx, xde), data = subset(soja, agua == 50))
sojas <- split(soja, f = soja$agua)
lapply(sojas, getInitial, object = rengrao ~ SSmrp(potassio, f0, tx, xde))
n0 <- lapply(sojas, function(a) {
    start <- getInitial(rengrao ~ SSmrp(potassio, f0, tx, xde), data = a)
    n0 <- nls(rengrao ~ f0 + tx * potassio * (potassio < xde) + tx * xde * (potassio >= 
        xde), data = a, start = start)
    n0
})
sapply(n0, coef)
lapply(n0, summary)
res <- do.call(c, lapply(n0, residuals))
qqmath(~res | soja$agua)
c0 <- lapply(n0, coef)
ic0 <- lapply(n0, confint.default)
xyplot(rengrao ~ potassio | agua, data = soja, col = 1, panel = function(...) {
    panel.xyplot(...)
    wp <- which.packet()
    with(as.list(c0[[wp]]), panel.curve(f0 + tx * x * (x <= xde) + tx * xde * (x > 
        xde), col = 2))
    panel.abline(v = c0[[wp]]["xde"], lty = 2, col = 2)
    panel.abline(v = ic0[[wp]]["xde", ], lty = 3, col = 2)
})
milk <- read.table("http://www.leg.ufpr.br/~walmes/data/saxton_lactacao1.txt", header = TRUE, 
    sep = "\t", colClasses = c("factor", NA, NA))
str(milk)
xyplot(leite ~ dia, groups = vaca, data = milk, type = "o")
milks <- subset(milk, vaca == "3")
p4 <- xyplot(leite ~ dia, groups = vaca, data = milks, col = 1, type = "o")
p4
vi.grid <- expand.grid(A = seq(5, 30, l = 10), B = seq(0, 1, l = 10), C = seq(0, 
    0.1, l = 50))
nrow(vi.grid)
n0 <- nls2(leite ~ A * dia^B * exp(-C * dia), data = milks, start = vi.grid, algorithm = "brute-force")
coef(n0)
n1 <- nls(formula(n0), data = milks, start = coef(n0))
summary(n1)
p4 + layer(with(as.list(coef(n1)), panel.curve(A * x^B * exp(-C * x), add = TRUE, 
    col = 2)))
