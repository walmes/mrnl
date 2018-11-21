#-----------------------------------------------------------------------
#                                            Prof. Dr. Walmes M. Zeviani
#                                leg.ufpr.br/~walmes · github.com/walmes
#                                        walmes@ufpr.br · @walmeszeviani
#                      Laboratory of Statistics and Geoinformation (LEG)
#                Department of Statistics · Federal University of Paraná
#                                       2018-nov-21 · Curitiba/PR/Brazil
#-----------------------------------------------------------------------

#=======================================================================
# Bandas de confiança.

#-----------------------------------------------------------------------
# Pacotes.

library(lattice)
library(latticeExtra)
library(rootSolve)
library(MASS)
library(nlme)

# Não está mais disponível. Usar `Matrix::bdiag()`.
# source("http://fisher.osu.edu/~schroeder_9/AMIS900/blockdiag.R")

#-----------------------------------------------------------------------
# Exemplo de uma curva apenas.

# Dados de liberação de potássio: conteúdo acumulado de potássio no
# tempo.
klib <- data.frame(k = c(51.03, 57.76, 26.6, 60.65, 87.07, 64.67, 91.28,
                         105.22, 72.74, 81.88, 97.62, 90.14, 89.88,
                         113.22, 90.91, 115.39, 112.63, 87.51, 104.69,
                         120.58, 114.32, 130.07, 117.65, 111.69, 128.54,
                         126.88, 127, 134.17, 149.66, 118.25, 132.67,
                         154.48, 129.11, 151.83, 147.66, 127.3),
                   t = rep(c(15, 30, 45, 60, 75, 90, 120, 150, 180, 210,
                             240, 270), each = 3))

# Diagrama de dispersão.
plot(k ~ t, data = klib)

expo.der <- deriv3(expr = ~A * (1 - exp(-log(2) * t/V)) + D * t,
                   namevec = c("A", "V", "D"),
                   function.arg = function(t, A, V, D) {
                       NULL
                   })
str(expo.der)

# Diagramas de dispersão.
p0 <- xyplot(k ~ t, data = klib, col = 1,
             xlab = "Período de incubacão (dias)",
             ylab = "Potássio liberado acumulado (mg/kg de solo)")
p0

# Verifica se os chutes iniciais são bons.
start <- list(A = 90, V = 20, D = 0.2)
p0 + layer(panel.curve(expo.der(x, A, V, D), add = TRUE, col = 2),
           data = start)

# Ajusta o modelo aos dados.
n0 <- nls(k ~ expo.der(t, A, V, D), data = klib, start = start)
summary(n0)

# Intervalo de confiança (de perfil de verossimilhanças).
confint(n0)

# Gráfico dos perfils de verossimilhanças.
par(mfrow = c(2, 2))
plot(profile(n0))
layout(1)

# Medidas de curvatura de Bates & Watts.
rms.curv(n0)

# Estrutura do modelo.
str(n0$m$fitted())

# Extração dos componentes.
c(n0$m$fitted())                # f(x, \theta).
attr(n0$m$fitted(), "gradient") # f'(x, \theta).
attr(n0$m$fitted(), "hessian")  # f''(x, \theta).

# Função que tem o vetor de parâmetros como primeiro argumento.
f <- function(theta, t) {
    with(as.list(theta),
         A * (1 - exp(-log(2) * t/V)) + D * t)
}

# Avalia o gradiente numericamente.
gradient(f, x = coef(n0), t = klib$t)

# Grid para obter os valores preditos.
pred <- data.frame(t = seq(0, 300, l = 100))

# Para construir os preditos e o gradiente avaliação em x.
f012 <- do.call(expo.der,
                args = c(list(t = pred$t),
                         as.list(coef(n0))))
str(f012)

# Extrai o gradiente.
f1 <- attr(f012, "gradient")

# Cholesky da matriz de covariância dos parâmetros (performance).
U <- chol(vcov(n0))

# Quantil da t para fazer o intervalo de confiança para predição de x.
tval <- qt(p = c(lwr = 0.025, fit = 0.5, upr = 0.975),
           df = df.residual(n0))

# Erro padrão de predição e margem de erro em cada ponto.
pred$se <- sqrt(apply(f1 %*% t(U), MARGIN = 1,
                      FUN = function(x) sum(x^2)))
pred$me <- outer(X = pred$se, Y = tval, FUN = "*")

# Concatena.
pred <- cbind(pred, sweep(x = pred$me, MARGIN = 1, STATS = f012, FUN = "+"))

# Gráfico com as bandas de confiança para a média \mu_i = f(x_i, \theta).
update(p0, xlim = c(0, NA), ylim = c(0, NA)) +
    as.layer(xyplot(fit + lwr + upr ~ t,
                    data = pred,
                    type = "l",
                    lty = c(1, 2, 2),
                    col = 1))

#-----------------------------------------------------------------------
# O caso de várias curvas simultaneamente ou curvas dentro de um
# delineamento experimental.

# Amostras de solo independentes foram colocadas para "secar" (perder
# água) no microondas. Foram usados intervalos de tempo de 4
# minutos. Foram estudados 4 solos diferentes com duas repetições por
# solo em cada tempo de permanência no microondas. A resposta é a
# umidade extraída da amostra de solo que é expressa de forma relativa,
# ou seja, proporcional a umidade conhecida adicionada ao solo. A
# pergunta da pesquisa é: 40 minutos é tempo suficiente para extrair
# toda a umidade da amostra de solo?

# Trabalho apresentado na Escola de Modelos de Regressão em 2011.
umi <- read.table("http://www.leg.ufpr.br/~walmes/data/emr11.txt",
                  header = TRUE, sep = "\t")
str(umi)

p1 <- xyplot(umrel ~ tempo | nome, data = umi, col = 1,
             ylab = "Umidade relativa do solo",
             xlab = "Tempo da amostra no microondas (min)")
p1 + layer(panel.abline(h = 1, lty = 2))

# Ajuste usando a logística com *self start*.
n0 <- nlsList(umrel ~ SSlogis(tempo, A, x0, S) | nome,
              data = umi)

# `nlsList()`: faz o ajuste dos modelos separadado pela variável
# condicionante.

# Estimativas.
coef(n0)
summary(n0)

# Grid para predição e valores preditos.
pred <- expand.grid(tempo = 0:45, nome = levels(umi$nome))
pred$umrel <- predict(n0, newdata = pred)
str(pred)

# Curvas ajustadas sobre os dados.
p1 +
    as.layer(xyplot(umrel ~ tempo | nome,
                    col = 2,
                    data = pred,
                    type = "l"))

# É recomendado fazer primeiros ajustes separadas para ganhar
# aprendizado sobre o modelo. Estratégia "dividir para
# conquistar". Depois de fazer os ajustes separados, use as estimativas
# como valores iniciais para o modelo conjunto. Pelo fato da dimensão
# aumentar (maior número de parâmetros), o ajuste de modelo conjuntos é
# mais complicado.

# Valores que serão usados como valores iniciais.
unlist(coef(n0))

# ATTENTION: os valores iniciais correspondem a estimativas por nível do
# fator condicionante, ou seja, correspondem a estimativas para fórmulas
# do tipo . ~ 0 + nome e são diferentes das estimativas sob a fórmula
# . ~ nome.

# Ajuste de modelo conjunto.
n1 <- gnls(umrel ~ SSlogis(tempo, A, x0, S),
           params = list(A ~ -1 + nome ,
                         x0 ~ -1 + nome,
                         S ~ -1 + nome),
           data = umi,
           start = unlist(coef(n0)))
summary(n1)

# IMPORTANT: veja a ordem dos parâmetros. Para obter a matriz de
# derivadas parciais, o efeito da variável condicionante tem que ser
# acomodado. A solução é usar operações de soma direta e matrizes e
# produto de Kronecker.

# Grid de combinações dos fatores para obter a predição.
pred <- expand.grid(tempo = seq(0, 50, length.out = 31),
                    nome = levels(umi$nome),
                    KEEP.OUT.ATTRS = FALSE)
pred$umrel <- predict(n1, newdata = pred)
str(pred)

# Gera o modelo para obter as derivadas parciais do modelo.
model <- deriv3(expr = ~A/(1 + exp(-(x - x0)/S)),
                namevec = c("A", "x0", "S"),
                function.arg = function(x, A, x0, S) NULL)

coef(n0)

# Faz uma lista com o vetor de estimativas de cada solo.
coef_list <- split(coef(n0), f = rownames(coef(n0)))
coef_list

# Determina a matriz de derivadas parciais.
i <- 0
f1 <- lapply(coef_list,
             FUN = function(theta) {
                 i <<- i + 1
                 x <- pred$tempo[pred$nome == levels(pred$nome)[i]]
                 f012 <- do.call(what = model,
                                 args = c(list(x = x),
                                          c(theta)))
                 # Retorna as derivadas parciais para o solo i.
                 attr(f012, "gradient")
             })

# É uma lista de matrizes conforme esperado da `lapply()`.
str(f1)

# Soma direta de matrizes.
bf1 <- as.matrix(Matrix::bdiag(f1))
str(bf1)

# DANGER: Cuidado com a order dos parâmetros.
by(bf1, pred$nome, FUN = function(x) x[1:3, ])

levelplot(t(bf1)[ncol(bf1):1, ]) +
    layer({
        n <- nlevels(pred$nome)
        u <- nrow(pred)/n
        v <- ncol(bf1)/n
        panel.abline(v = seq(from = v, by = v, length.out = n - 1) + 0.5,
                     h = seq(from = u, by = u, length.out = n - 1) + 0.5,
                     col = 1, lty = 2)
    })

# Ordem dos parâmetros no modelo ajustado.
colnames(vcov(n1))

# Índices para reordenar as colunas da matriz para a ordem correta.
# TIP: use uma matriz!
ro <- matrix(1:ncol(bf1), nrow = length(f1), byrow = TRUE)
c(ro)

# Muda para a ordem correta dos parâmetros.
bf1 <- bf1[, c(ro)]
head(bf1)

# Quantil da t.
tval <- qt(p = c(lwr = 0.025, fit = 0.5, upr = 0.975),
           df = nrow(umi) - ncol(bf1))

# Cholesky da covariância das estimativas.
U <- chol(vcov(n1))

pred$se <- sqrt(apply(bf1 %*% t(U),
                      MARGIN = 1,
                      FUN = function(x) sum(x^2)))
pred$me <- outer(X = pred$se, Y = tval, FUN = "*")
pred <- cbind(pred, sweep(x = pred$me,
                          MARGIN = 1,
                          STATS = pred$umrel,
                          FUN = "+"))
pred$me <- NULL
str(pred)

update(p1, xlim = c(0, 50)) +
    as.layer(
        xyplot(fit + lwr + upr ~ tempo | nome,
               data = pred,
               type = "l",
               col = 1,
               lty = c(1, 2, 2),
               lwd = c(1.5, 1, 1))) +
    layer({
        panel.abline(v = seq(0, 50, 5),
                     h = seq(0, 1, 0.1),
                     lty = 3,
                     col = "gray90")
        panel.abline(h = 1, v = 40, lty = 3, lwd = 2, col = "orange")
    }, under = TRUE)

# Intervalos de confiança para a predição no tempo 40.
subset(pred, tempo == 40)

#-----------------------------------------------------------------------
