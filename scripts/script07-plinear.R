library(lattice)
library(latticeExtra)
library(nls2)
data(wtloss, package = "MASS")
p0 <- xyplot(Weight ~ Days, data = wtloss, xlab = "Dias em dieta", ylab = "Peso (kg)")
print(p0)
p0 + layer(with(list(f0 = 81, f1 = 102, K = 141), panel.curve(f0 + f1 * 2^(-x/K), 
    add = TRUE, col = 2)))
n0 <- nls(Weight ~ f0 + f1 * 2^(-Days/K), data = wtloss, start = list(f0 = 81, f1 = 102, 
    K = 141), trace = TRUE)
summary(n0)
n1 <- nls(Weight ~ cbind(f0 = 1, f1 = 2^(-Days/K)), data = wtloss, start = list(K = 141), 
    trace = TRUE, algorithm = "plinear")
summary(n1)
confint.default(n1)
p0 + layer(with(as.list(coef(n0)), panel.curve(f0 + f1 * 2^(-x/K), add = TRUE, col = 2)))
data(turk0, package = "alr3")
str(turk0)
p1 <- xyplot(Gain ~ A, data = turk0, col = 1, xlab = "Nutriente", ylab = "Ganho de peso")
print(p1)
p1 + layer(with(list(f0 = 620, f1 = 160, k = 9), panel.curve(f0 + f1 * (1 - exp(-k * 
    x)), add = TRUE, col = 2)))
n0 <- nls(Gain ~ f0 + f1 * (1 - exp(-k * A)), data = turk0, start = list(f0 = 620, 
    f1 = 160, k = 9))
summary(n0)
D(expression(f0 + f1 * (1 - exp(-k * A))), "k")
n1 <- nls(Gain ~ cbind(1, 1 - exp(-k * A)), data = turk0, start = list(k = 9), trace = TRUE, 
    algorithm = "plinear")
summary(n1)
p1 + layer(with(as.list(coef(n0)), panel.curve(f0 + f1 * (1 - exp(-k * x)), add = TRUE, 
    col = 2)))
n0 <- nls(Gain ~ b0 + b1 * A * (A < xb) + b1 * xb * (A >= xb), data = turk0, start = list(b0 = 600, 
    b1 = 700, xb = 0.2), trace = TRUE)
n1 <- nls(Gain ~ cbind(b0 = 1, b1 = A * (A < xb) + xb * (A >= xb)), data = turk0, 
    start = list(xb = 0.2), trace = TRUE, algorithm = "plinear")
summary(n1)
p1 + layer(with(as.list(coef(n0)), panel.curve(b0 + b1 * x * (x < xb) + b1 * xb * 
    (x >= xb), add = TRUE, col = 2)))
cra <- read.table("http://www.leg.ufpr.br/~walmes/data/cra_manejo.txt", header = TRUE, 
    sep = "\t")
str(cra)
cra$tens[cra$tens == 0] <- 0.1
xyplot(umid ~ log10(tens) | condi * posi, groups = prof, data = cra, auto.key = TRUE)
cras <- droplevels(subset(cra, condi == "B"))
p2 <- xyplot(umid ~ tens, data = cras, col = 1, ylab = expression("Conteúdo de água do solo" ~ 
    (g ~ g^{
        -1
    })), xlab = expression("Tensão matricial (kPa)"), scales = list(x = list(log = 10)))
p2
start <- list(tr = 0.3, ts = 0.6, al = 0.5, n = 1.3)
p2 + layer(with(start, panel.curve(tr + (ts - tr)/(1 + (al * 10^x)^n)^(1 - 1/n), 
    add = TRUE, col = 2)))
D(expression(tr + (ts - tr)/(1 + (al * tens)^n)^(1 - 1/n)), "tr")
D(expression(tr + (ts - tr)/(1 + (al * tens)^n)^(1 - 1/n)), "ts")
str(cras)
n0 <- nls(umid ~ tr + (ts - tr)/(1 + (al * tens)^n)^(1 - 1/n), data = subset(cras), 
    start = start, trace = TRUE)
summary(n0)
n1 <- nls(umid ~ cbind(tr = 1 - 1/(1 + (al * tens)^n)^(1 - 1/n), ts = 1/(1 + (al * 
    tens)^n)^(1 - 1/n)), data = cras, start = list(al = 0.5, n = 1.3), trace = TRUE, 
    algorithm = "plinear")
summary(n1)
p2 + layer(with(as.list(coef(n0)), panel.curve(tr + (ts - tr)/(1 + (al * 10^x)^n)^(1 - 
    1/n), add = TRUE, col = 2)))
