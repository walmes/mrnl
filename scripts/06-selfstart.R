#-----------------------------------------------------------------------
#                                            Prof. Dr. Walmes M. Zeviani
#                                leg.ufpr.br/~walmes · github.com/walmes
#                                        walmes@ufpr.br · @walmeszeviani
#                      Laboratory of Statistics and Geoinformation (LEG)
#                Department of Statistics · Federal University of Paraná
#                                       2018-nov-21 · Curitiba/PR/Brazil
#-----------------------------------------------------------------------

#=======================================================================
# Automação de valores iniciais.

#-----------------------------------------------------------------------
# Pacotes.

library(lattice)
library(latticeExtra)
library(car)
library(segmented)
library(nls2)

#-----------------------------------------------------------------------

data(wtloss, package = "MASS")

# Diagrama de dispersão dos dados.
p0 <- xyplot(Weight ~ Days, data = wtloss, col = 1,
             xlab = "Dias em dieta", ylab = "Peso (kg)")
print(p0)

# Adicionando nossa "adivinhação".
p0 +
    layer({
        panel.curve(f0 + f1 * 2^(-x/K),
                    add = TRUE, col = 2)
    },
    data = list(f0 = 81, f1 = 102, K = 141))

# Formalização matemática na nossa estimativa visual.
pe <- with(wtloss, {
    f1 <- diff(range(Weight))
    f0 <- min(Weight)
    m <- f0 + 0.5 * f1
    mm <- which.min(abs(Weight - m))[1]
    K <- Days[mm]
    list(f0 = f0, f1 = f1, K = K)
})

p0 +
    layer(panel.abline(v = K, h = c(f0 + c(0, 0.5, 1) * f1), lty = 2),
          data = pe)

# Cria uma função que determina os chutes a partir dos dados.
model.init <- function(x, y) {
    f1 <- diff(range(y))
    f0 <- min(y)
    m <- f0 + 0.5 * f1
    mm <- which.min(abs(y - m))[1]
    K <- x[mm]
    list(f0 = f0, f1 = f1, K = K)
}
model.init(x = wtloss$Days, y = wtloss$Weight)

# Cria a função de cabeçalho padrão para uso pela `nls()`.
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
    return(value)
}

# Cria uma função *self start* (SS).
SSmodel <- selfStart(model = ~f0 + f1 * 2^(-input/K),
                     parameters = c("f0", "f1", "K"),
                     initial = model.init,
                     template = function(input, f0, f1, K) {
                         NULL
                     })
SSmodel

# ATTENTION: o modelo tem o gradiente determinado. A função
# `model.init()` é armazenada como atributo.

# Usa a função para ter o valores inicias "estimados".
getInitial(Weight ~ SSmodel(Days, f0, f1, K),
           data = wtloss)

# O `trace = TRUE` mostra de onde partiu a otimização.
n0 <- nls(Weight ~ SSmodel(Days, f0, f1, K),
          data = wtloss,
          trace = TRUE)
summary(n0)

#-----------------------------------------------------------------------
# Um exemplo com modelo polinômio quadrático na forma canônica.

# Importa os dados.
cm <- read.table("http://www.leg.ufpr.br/~walmes/data/cresmicelial.txt",
                 header = TRUE,
                 sep = "\t",
                 colClasses = c("factor", "numeric", "numeric"))
str(cm)

# Diagrama de dispersão.
p1 <- xyplot(mmdia ~ temp | isolado, data = cm, type = c("p", "a"),
             col = 1,
             ylab = "Taxa de crescimento em diâmetro do disco micelial",
             xlab = "Temperatuda da câmara de crescimento")
print(p1)

# Função para automação de valores inicias.
mqc.init <- function(mCall, data, LHS) {
    xy <- data.frame(sortedXyData(mCall[["input"]], LHS, data))
    # Ajuste do modelo de segundo grau.
    m0 <- lm(y ~ poly(x, 2, raw = TRUE), data = xy)
    c0 <- coef(m0)
    cv <- c0[3]                        # Curvatura.
    xc <- -0.5 * c0[2]/cv              # Ponto estacionário.
    fc <- c0[1] + c0[2] * xc + cv * xc # Função no ponto estacionário.
    value <- c(fc, xc, cv)
    names(value) <- mCall[c("fc", "xc", "cv")]
    return(value)
}

# Modelo com *self start* como atributo.
SSmqc <- selfStart(model = ~fc + cv * (input - xc)^2,
                   parameters = c("fc", "xc", "cv"),
                   initial = mqc.init,
                   template = function(input, fc, xc, cv) {
                       NULL
                   })

# Verifica o funcionamento.
getInitial(mmdia ~ SSmqc(temp, fc, xc, cv),
           data = subset(cm, isolado == "1"))

# Partições por isolado para determinar valores iniciais.
cms <- split(cm, f = cm$isolado)
sapply(cms,
       FUN = getInitial,
       object = mmdia ~ SSmqc(temp, fc, xc, cv))

# Ajusta o modelo separado por isolado.
n0 <- lapply(cms,
             FUN = nls,
             formula = mmdia ~ SSmqc(temp, fc, xc, cv))

# Extração das estimativas.
c0 <- t(sapply(n0, coef))
c0

# Extração dos intervalos de confiança.
ic0 <- lapply(n0, confint.default)
ic0

# Extração dos resíduos.
res <- do.call(c, lapply(n0, residuals))
qqmath(~res | cm$isolado)

# Gráfico das curvas ajustadas.
c0 <- lapply(n0, coef)
xyplot(mmdia ~ temp | isolado, data = cm, col = 1,
       ylab = "Taxa de crescimento em diâmetro do disco micelial",
       xlab = "Temperatuda da câmara de crescimento") +
    layer({
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

segplot(i ~ lwr + upr,
        data = ics,
        col = 1,
        xlab = "Temperatura ótima para o crescimento",
        ylab = "Isolado",
        draw.bands = FALSE)

# Ajuste do modelo linear.
m0 <- lm(mmdia ~ temp + I(temp^2), cms[[1]])
n0 <- nls(mmdia ~ SSmqc(temp, fc, xc, cv), cms[[1]])

# Uso do método delta para estimar IC para funções dos parâmetros.
deltaMethod(m0,
            g = "-0.5*b1/b2",
            parameterNames = c("b0", "b1", "b2"))

summary(n0)$coeff["xc", 1:2]

# Funções *self start* disponíveis.
apropos("^SS", where = TRUE)

#-----------------------------------------------------------------------
# Modelo segmentado.

soja <- read.table("http://www.leg.ufpr.br/~walmes/data/soja.txt",
                   header = TRUE,
                   sep = "\t",
                   dec = ",")
soja <- soja[-74, ] # Remove outlier.
str(soja)

xyplot(rengrao ~ potassio | factor(agua), data = soja, col = 1,
       xlab = "Conteúdo de potássio no solo",
       ylab = "Produção de soja por vaso")

xyplot(rengrao ~ potassio, groups = agua, data = soja,
       type = c("p", "a"), xlab = "Conteúdo de potássio no solo",
       ylab = "Produção de soja por vaso")

mrp.init <- function(mCall, data, LHS) {
    xy <- data.frame(sortedXyData(mCall[["input"]], LHS, data))
    m0 <- lm(y ~ x, data = xy)
    xi <- mean(xy[["x"]])
    # Usa o pacote `segmented` para determinar o ponto de quebra.
    # ATTENTION: ele pode não convergir.
    s0 <- segmented(m0,
                    seg.Z = ~x,
                    psi = list(x = c(xi)),
                    control = seg.control(display = FALSE))
    xde <- s0$psi[1, 2]     # Ponto de quebra ou mudança de taxa.
    f0 <- coef(s0)[1]       # Intercepto: f(x = 0).
    tx <- coef(s0)[2]       # Taxa de incremento: f(x = u + 1) - f(x = u).
    value <- c(f0, tx, xde)
    names(value) <- mCall[c("f0", "tx", "xde")]
    return(value)
}

# Função *self start* para o modelo linear-platô.
SSmrp <- selfStart(model = ~f0 + tx * input + xde,
                   parameters = c("f0", "tx", "xde"),
                   initial = mrp.init,
                   template = function(input, f0, tx, xde) {
                       NULL
                   })
body(SSmrp)

# Verifica funcionamento.
getInitial(rengrao ~ SSmrp(potassio, f0, tx, xde),
           data = subset(soja, agua == 50))

# Ajusta o modelo em cada nível de água.
sojas <- split(soja, f = soja$agua)
lapply(sojas,
       FUN = getInitial,
       object = rengrao ~ SSmrp(potassio, f0, tx, xde))

n0 <- lapply(sojas,
             FUN = function(a) {
                 # Gera valores inicias.
                 start <- getInitial(rengrao ~ SSmrp(potassio, f0, tx, xde),
                                     data = a)
                 # Faz o ajuste.
                 n0 <- nls(rengrao ~ f0 + tx * potassio * (potassio < xde) +
                               tx * xde * (potassio >= xde),
                           data = a,
                           start = start)
                 n0
             })

# Extração de estimativas.
sapply(n0, coef)
lapply(n0, summary)

# Extração e inspeção dos resíduos.
res <- do.call(c, lapply(n0, residuals))
qqmath(~res | soja$agua)

# Intervalos de confiança.
c0 <- lapply(n0, coef)
ic0 <- lapply(n0, confint.default)

# Curvas ajustadas.
xyplot(rengrao ~ potassio | agua,
       data = soja,
       col = 1) +
    layer({
        wp <- which.packet()
        with(as.list(c0[[wp]]),
             panel.curve(f0 + tx * x * (x <= xde) +
                         tx * xde * (x > xde),
                         col = 2))
        panel.abline(v = c0[[wp]]["xde"], lty = 2, col = 2)
        panel.abline(v = ic0[[wp]]["xde", ], lty = 3, col = 2)
    })

#-----------------------------------------------------------------------
# Método força bruta.

# TIP: existem modelo não lineares que não possuem interpretação prática
# ou matemática/cartesiana tão simples. Nessa condição, inspeção visual
# para determinar valores iniciais é bem complicado. Então pode-se usar
# a força bruta, ou seja, passar um grid de valores inicias e usar o
# melhor ponto de suporte do grid como valores iniciais.

# Importação de dados de curva de lactação.
milk <- read.table("http://www.leg.ufpr.br/~walmes/data/saxton_lactacao1.txt",
                   header = TRUE, sep = "\t",
                   colClasses = c("factor", NA, NA))
str(milk)

# Diagrama de dispersão.
xyplot(leite ~ dia, groups = vaca, data = milk, type = "o")

milks <- subset(milk, vaca == "3")
p4 <- xyplot(leite ~ dia, groups = vaca, data = milks, col = 1, type = "o")
p4

# Grid de valores iniciais.
vi.grid <- expand.grid(A = seq(5, 30, l = 10),
                       B = seq(0, 1, l = 10),
                       C = seq(0, 0.1, l = 50))
nrow(vi.grid)

# Ajuste para uma sintonia grossa (aproximação de valores iniciais).
n0 <- nls2(leite ~ A * dia^B * exp(-C * dia),
           data = milks,
           start = vi.grid,
           algorithm = "brute-force")
coef(n0)

# ATTENTION: `n0` estimou os parâmetros dentro do grid. Não
# necessariamente é o ótimo verdadeiro mas espera-se que esteja próximo
# para ser um chute inicial bem sucedido.

# Ajuste para *sintonia fina* (fine tunning).
n1 <- nls(formula(n0),
          data = milks,
          start = coef(n0))

summary(n1)

p4 +
    layer({
        panel.curve(A * x^B * exp(-C * x), add = TRUE, col = 2)
    },
    data = as.list(coef(n1)))

#-----------------------------------------------------------------------
# Mais sobre valores inicias.
#
# REVIEW: revisar o nome dos pacotes.
#
# library(drc)
# library(HydroME)

#-----------------------------------------------------------------------
