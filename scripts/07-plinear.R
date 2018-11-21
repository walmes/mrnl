#-----------------------------------------------------------------------
#                                            Prof. Dr. Walmes M. Zeviani
#                                leg.ufpr.br/~walmes · github.com/walmes
#                                        walmes@ufpr.br · @walmeszeviani
#                      Laboratory of Statistics and Geoinformation (LEG)
#                Department of Statistics · Federal University of Paraná
#                                       2018-nov-21 · Curitiba/PR/Brazil
#-----------------------------------------------------------------------

#=======================================================================
# Modelo não lineares parcialmente lineares.

#-----------------------------------------------------------------------
# Pacotes.

library(lattice)
library(latticeExtra)
library(nls2)

data(wtloss, package = "MASS")
p0 <- xyplot(Weight ~ Days,
             data = wtloss,
             xlab = "Dias em dieta",
             ylab = "Peso (kg)")
print(p0)

start <- list(f0 = 81, f1 = 102, K = 141)

p0 + layer(panel.curve(f0 + f1 * 2^(-x/K), add = TRUE, col = 2),
           data = start)

n0 <- nls(Weight ~ f0 + f1 * 2^(-Days/K),
          data = wtloss,
          start = start,
          trace = TRUE)

summary(n0)

f <- expression(f0 + f1 * 2^(-Days/K))
D(f, "f0") # Linear em f0 (não envolve f0).
D(f, "f1") # Linear em f1 (não envolve f1).
D(f, "K")  # Não é linear.

n1 <- nls(Weight ~ cbind(f0 = 1, f1 = 2^(-Days/K)),
          data = wtloss,
          start = list(K = 141),
          trace = TRUE,
          algorithm = "plinear")

summary(n1)
confint.default(n1)

p0 + layer(panel.curve(f0 + f1 * 2^(-x/K), add = TRUE, col = 2),
           data = as.list(coef(n0)))

#-----------------------------------------------------------------------

data(turk0, package = "alr3")
str(turk0)

p1 <- xyplot(Gain ~ A,
             data = turk0,
             col = 1,
             xlab = "Nutriente",
             ylab = "Ganho de peso")
print(p1)

start <- list(f0 = 620, f1 = 160, k = 9)
p1 + layer(panel.curve(f0 + f1 * (1 - exp(-k * x)), add = TRUE, col = 2),
           data = start)

n0 <- nls(Gain ~ f0 + f1 * (1 - exp(-k * A)),
          data = turk0,
          start = start)

summary(n0)

f <- expression(f0 + f1 * (1 - exp(-k * A)))
D(f, "f0")
D(f, "f1")
D(f, "k")

n1 <- nls(Gain ~ cbind(1, 1 - exp(-k * A)),
          data = turk0,
          start = list(k = 9),
          trace = TRUE,
          algorithm = "plinear")
summary(n1)

p1 +
    layer(panel.curve(f0 + f1 * (1 - exp(-k * x)), add = TRUE, col = 2),
          data = as.list(coef(n0)))

n0 <- nls(Gain ~ b0 + b1 * A * (A < xb) + b1 * xb * (A >= xb),
          data = turk0,
          start = list(b0 = 600, b1 = 700, xb = 0.2),
          trace = TRUE)

f <- expression(b0 + b1 * A * (A < xb) + b1 * xb * (A >= xb))
D(f, "b0") # Derivadas analíticas não sabem tratar a função indicadora.

n1 <- nls(Gain ~ cbind(b0 = 1, b1 = A * (A < xb) + xb * (A >= xb)),
          data = turk0,
          start = list(xb = 0.2),
          trace = TRUE,
          algorithm = "plinear")
summary(n1)

# Não fica um bom ajuste.
p1 +
    layer(panel.curve(b0 + b1 * x * (x < xb) + b1 * xb * (x >= xb),
                      add = TRUE, col = 2),
          data = as.list(coef(n0)))

#-----------------------------------------------------------------------
# Modelo van Genuchten.

cra <- read.table("http://www.leg.ufpr.br/~walmes/data/cra_manejo.txt",
                  header = TRUE,
                  sep = "\t")
str(cra)
cra$tens[cra$tens == 0] <- 0.1

xyplot(umid ~ log10(tens) | condi * posi,
       groups = prof,
       data = cra,
       auto.key = TRUE)

cras <- droplevels(subset(cra, condi == "B"))

# Diagrama de dispersão.
p2 <- xyplot(umid ~ tens, data = cras, col = 1,
             ylab = expression("Conteúdo de água do solo" ~ (g ~ g^{-1})),
             xlab = expression("Tensão matricial (kPa)"),
             scales = list(x = list(log = 10)))
p2

# Valores iniciais.
start <- list(tr = 0.3, ts = 0.6, al = 0.5, n = 1.3)
p2 + layer(panel.curve(tr + (ts - tr)/(1 + (al * 10^x)^n)^(1 - 1/n),
                       add = TRUE, col = 2),
           data = start)

# Derivadas parciais.
f <- expression(tr + (ts - tr)/(1 + (al * tens)^n)^(1 - 1/n))
D(f, "tr") # Linear.
D(f, "ts") # Linear.
D(f, "al")
D(f, "n")

# Modelo que otimiza em 4 parâmetros.
n0 <- nls(umid ~ tr + (ts - tr)/(1 + (al * tens)^n)^(1 - 1/n),
          data = subset(cras),
          start = start,
          trace = TRUE)
summary(n0)

# Modelo que otimiza apenas nos 2 não lineares.
n1 <- nls(umid ~ cbind(tr = 1 - 1/(1 + (al * tens)^n)^(1 - 1/n),
                       ts = 1/(1 + (al * tens)^n)^(1 - 1/n)),
          data = cras,
          start = list(al = 0.5, n = 1.3),
          trace = TRUE,
          algorithm = "plinear")
summary(n1)

# Verifica o ajuste sobre os dados.
p2 + layer(panel.curve(tr + (ts - tr)/(1 + (al * 10^x)^n)^(1 - 1/n),
                       add = TRUE, col = 2),
           data = as.list(coef(n0)))

#-----------------------------------------------------------------------
