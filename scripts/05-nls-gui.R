#-----------------------------------------------------------------------
#                                            Prof. Dr. Walmes M. Zeviani
#                                leg.ufpr.br/~walmes · github.com/walmes
#                                        walmes@ufpr.br · @walmeszeviani
#                      Laboratory of Statistics and Geoinformation (LEG)
#                Department of Statistics · Federal University of Paraná
#                                       2018-nov-18 · Curitiba/PR/Brazil
#-----------------------------------------------------------------------

#=======================================================================
# Método gráfico interativo para ajustar modelos.

library(rpanel)
library(wzRfun)
library(latticeExtra)

# Se não possuir o pacote wzRfun, carregue a função da fonte.
source("https://raw.githubusercontent.com/walmes/wzRfun/master/R/rp.nls.R")

# Verifica a documentação.
help(np.nls, help_type = "html")

#-----------------------------------------------------------------------
# Ganho de peso dos perus.

data(turk0, package = "alr3")
str(turk0)

plot(Gain ~ A, data = turk0,
     xlab = "Metionine", ylab = "Weight gain")

turk.fit <- rp.nls(model = Gain ~ Int + (Asy - Int) * A/(Half + A),
                   data = turk0,
                   start = list(Int = c(600, 650),
                                Asy = c(750, 850),
                                Half = c(0, 0.2)),
                   extra = function(Int, Asy, Half) {
                       abline(h = c(Asy, Int), v = Half,
                              col = "green")
                   },
                   xlab = "Metionine", ylab = "Weight gain")

summary(turk.fit)
confint(turk.fit)

f <- formula(turk.fit)

plot(Gain ~ A, data = turk0,
     xlab = "Metionine", ylab = "Weight gain")
with(as.list(coef(turk.fit)), {
    eval(call("curve",
              f[[3]],
              add = TRUE,
              col = "purple",
              xname = intersect(all.vars(f[-2]), names(turk0))))
})

#-----------------------------------------------------------------------
# Ajuste de curva por grupos.

xyplot(rate ~ conc, groups = state, data = Puromycin,
       type = c("p", "smooth"), auto.key = TRUE)

Puro.fit <- rp.nls(model = rate ~ Int + (Top - Int) * conc/(Half + conc),
                   data = Puromycin,
                   start = list(Int = c(20, 70),
                                Top = c(100, 200),
                                Half = c(0, 0.6)),
                   subset = "state",
                   gpar = list(try = list(col = 2, lty = 2),
                               fit = list(col = "blue", lwd = 1.5)),
                   xlab = "Concentration", ylab = "Rate",
                   xlim = c(0, 1.2), ylim = c(40, 220))

length(Puro.fit)
sapply(Puro.fit, coef)
lapply(Puro.fit, confint.default)
sapply(Puro.fit, logLik)
sapply(Puro.fit, deviance)

#-----------------------------------------------------------------------
# Ajuste de curvas de retenção de água do solo.

cra <- read.table("http://www.leg.ufpr.br/~walmes/data/cra_manejo.txt",
                  header = TRUE, sep = "\t")

cra$tens[cra$tens == 0] <- 0.1
cras <- subset(cra, condi == "LVA3,5")
cras <- aggregate(umid ~ posi + prof + tens, data = cras, FUN = mean)
cras$caso <- with(cras, interaction(posi, prof))
cras$ltens <- log(cras$tens)

xyplot(umid ~ ltens | posi, groups = prof, data = cras,
       type = c("p", "a"))

# modelo : van Genuchten com retrição de Mualem.
# ltens  : representado por ltens (log da tensão matricial, psi).
# umid   : representado por umid, conteúdo de água do solo ( % ).
# Us     : assíntota inferior, mínimo da função, quando x -> +infinito.
# Ur     : assíntota superior, máximo da função, quando x -> -infinito.
# al     : locação, translada o ponto de inflexão.
# n      : forma, altera a taxa ao redor da inflexão.

model <- umid ~ Ur + (Us - Ur)/(1 + exp(n * (al + ltens)))^(1 - 1/n)

start <- list(Ur = c(init = 0.2, from = 0, to = 0.5),
              Us = c(init = 0.6, from = 0.4, to = 0.8),
              al = c(1, -2, 4),
              n = c(1.5, 1, 4))

cra.fit <- rp.nls(model = model,
                  data = cras,
                  start = start,
                  subset = "caso")

# Guarde os coeficientes para não precisar repetir o processo no
# futuro. Use `dput()`.
t(sapply(cra.fit, coef))

#-----------------------------------------------------------------------
# Curva de produção em função da desfolha do algodão.

cap <- read.table("http://www.leg.ufpr.br/~walmes/data/algodão.txt",
                  header = TRUE, sep = "\t", encoding = "latin1")

cap$desf <- cap$desf/100
cap <- subset(cap, select = c(estag, desf, pcapu))
cap$estag <- factor(cap$estag, labels = c("vegetativo", "botão floral",
                                          "florescimento", "maçã",
                                          "capulho"))
str(cap)

xyplot(pcapu ~ desf | estag, data = cap, layout = c(5, 1),
       xlab = "Nível de desfolha artificial", ylab = "Peso de capulhos")

# modelo: potência.
# desf  : representado por desf (nível de desfolha artifical).
# pcapu : representado por pcapu (peso de capulhos), produto do algodão.
# f0    : intercepto, valor da função quando x = 0 (situação ótima).
# delta : diferença no valor da função para x = 0 e x = 1 (situação
#         extrema).
# curv  : forma, indica como a função decresce, se th3 = 0 então função
#         linear.

model <- pcapu ~ f0 - delta * desf^exp(curv)
start <- list(f0 = c(30, 25, 35),
              delta = c(8, 0, 16),
              curv = c(0, -2, 4))

cap.fit <- rp.nls(model = model,
                  data = cap,
                  start = start,
                  subset = "estag")

length(cap.fit)
sapply(cap.fit, coef)
lapply(cap.fit, confint.default)

# Curvas sobre os valores observados.
xyplot(pcapu ~ desf | estag, data = cap, layout = c(5, 1),
       xlab = "Nível de desfolha artificial",
       ylab = "Peso de capulhos") +
    layer({
        n0 <- cap.fit[[which.packet()]]
        xx <- seq(0, 1, length = 21)
        yy <- predict(n0, newdata = list(desf = xx))
        panel.lines(xx, yy, col = 2)
    })

#-----------------------------------------------------------------------
# Curva de produção em função do nível de potássio no solo.

soja <- read.table("http://www.leg.ufpr.br/~walmes/data/soja.txt",
                   header = TRUE, sep = "\t", encoding = "latin1",
                   dec = ",")

soja$agua <- factor(soja$agua)
str(soja)

xyplot(rengrao ~ potassio | agua, data = soja)

# modelo : linear-plato.
# x      : representado por potássio, conteúdo de potássio do solo.
# y      : representado por rengrao, redimento de grãos por parcela.
# f0     : intercepto, valor da função quando x=0.
# tx     : taxa de incremento no rendimento por unidade de x.
# brk    : valor acima do qual a função é constante.

model <- rengrao ~ f0 + tx * potassio * (potassio < brk) +
    tx * brk * (potassio >= brk)
start <- list(f0 = c(15, 5, 25),
              tx = c(0.2, 0, 1),
              brk = c(50, 0, 180))

pot.fit <- rp.nls(model = model, data = soja, start = start,
                  subset = "agua")

sapply(pot.fit, coef)
lapply(pot.fit, confint.default)

# Curvas sobre os valores observados.
xyplot(rengrao ~ potassio | agua, data = soja) +
    layer({
        n0 <- pot.fit[[which.packet()]]
        xx <- seq(min(x), max(x), length = 51)
        yy <- predict(n0, newdata = list(potassio = xx))
        panel.lines(xx, yy, col = 2)
    })

#-----------------------------------------------------------------------
# Curva de lactação.

lac <- read.table("http://www.leg.ufpr.br/~walmes/data/saxton_lactacao1.txt",
                  header = TRUE, sep = "\t", encoding = "latin1")

lac$vaca <- factor(lac$vaca)
str(lac)

xyplot(leite ~ dia | vaca, data = lac)

# modelo : de Wood (nucleo da densidade gama).
# x      : representado por dia, dia após parto.
# y      : representado por leite, quantidade produzida.
# th1    : escala, desprovido de interpretação direta.
# th2    : forma, desprovido de interpretação direta.
# th3    : forma, desprovido de interpretação direta.

# É o núcleo da função de densidade da distribuição Gama.
model <- leite ~ th1 * dia^th2 * exp(-th3 * dia)
start <- list(th1 = c(15, 10, 20),
              th2 = c(0.2, 0.05, 0.5),
              th3 = c(0.0025, 0.001, 0.008))

lac.fit <- rp.nls(model = model,
                  data = lac,
                  start = start,
                  subset = "vaca",
                  xlim = c(0, 310))

sapply(lac.fit, coef)
lapply(lac.fit, confint.default)

xyplot(leite ~ dia | vaca, data = lac) +
    layer({
        n0 <- lac.fit[[which.packet()]]
        xx <- seq(min(x), max(x), length = 51)
        yy <- predict(n0, newdata = list(dia = xx))
        panel.lines(xx, yy, col = 2)
    })

#-----------------------------------------------------------------------
# Curvas de crescimento em placa de petri.

cre <- read.table("http://www.leg.ufpr.br/~walmes/data/cresmicelial.txt",
                  header = TRUE, sep = "\t", encoding = "latin1")

cre$isolado <- factor(cre$isolado)
cre$mmdia <- sqrt(cre$mmdia)
str(cre)

xyplot(mmdia ~ temp | isolado, data = cre)

# modelo : quadrático na forma canônica.
# x      : representado por temp, temperatura da câmara de crescimento.
# y      : representado por mmdia, taxa média de crescimento.
# thy    : valor da função no ponto de máximo.
# thc    : curvatura ou grau de especificidade à condição ótima.
# thx    : ponto de máximo, temperatura de crescimento mais rápido.

model <- mmdia ~ thy + thc * (temp - thx)^2
start <- list(thy = c(4, 0, 7),
              thc = c(-0.05, 0, -0.5),
              thx = c(23, 18, 30))

mic.fit <- rp.nls(model = model,
                  data = cre,
                  start = start,
                  subset = "isolado",
                  xlim = c(17, 31),
                  ylim = c(0, 6))

t(sapply(mic.fit, coef))

xyplot(mmdia ~ temp | isolado, data = cre) +
    layer({
        n0 <- mic.fit[[which.packet()]]
        xx <- seq(min(x), max(x), length = 51)
        yy <- predict(n0, newdata = list(temp = xx))
        panel.lines(xx, yy, col = 2)
    })

#-----------------------------------------------------------------------
# Curva de secagem do solo em microondas.

sec <- read.table("http://www.leg.ufpr.br/~walmes/data/emr11.txt",
                  header = TRUE, sep = "\t", encoding = "latin1")
str(sec)

xyplot(umrel ~ tempo | nome, data = sec)

# modelo : logístico.
# x      : representado por tempo, período da amostra dentro do
#          microondas.
# y      : representado por umrel, umidade relativa o conteúdo total
#          de água.
# th1    : assíntota superior.
# th2    : tempo para evaporar metade do conteúdo total de água.
# th3    : proporcional à taxa máxima do processo.

model <- umrel ~ th1/(1 + exp(-(tempo - th2)/th3))
start <- list(th1 = c(1, 0.8, 1.2),
              th2 = c(15, 0, 40),
              th3 = c(8, 2, 14))

sec.fit <- rp.nls(model = model,
                  data = sec,
                  start = start,
                  subset = "nome")

sapply(sec.fit, coef)
lapply(sec.fit, confint.default)

xyplot(umrel ~ tempo | nome, data = sec) +
    layer({
        n0 <- sec.fit[[which.packet()]]
        xx <- seq(min(x), max(x), length = 51)
        yy <- predict(n0, newdata = list(tempo = xx))
        panel.lines(xx, yy, col = 2)
    })

#-----------------------------------------------------------------------
