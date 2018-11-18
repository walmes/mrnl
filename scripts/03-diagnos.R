#-----------------------------------------------------------------------
#                                            Prof. Dr. Walmes M. Zeviani
#                                leg.ufpr.br/~walmes · github.com/walmes
#                                        walmes@ufpr.br · @walmeszeviani
#                      Laboratory of Statistics and Geoinformation (LEG)
#                Department of Statistics · Federal University of Paraná
#                                       2018-nov-17 · Curitiba/PR/Brazil
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
# Pacotes.

library(tidyverse)
library(rpanel)
# library(proto)
library(wzRfun) # R2() e as.lm().

# Funções úteis para o Curso.
source("_funs.R")

#-----------------------------------------------------------------------
# Análise exploratória.

data(turk0, package = "alr3")
str(turk0)

turk0 %>%
    group_by(A) %>%
    summarise(n = n(),
              Gain_m = mean(Gain))

# QUESTION: Em quais aspectos fazer a análise com as médias amostrais
# difere da análise com todas as observações?  Se o peso `n` for
# incorporado a análise, existe ainda diferença em algum aspecto?

p0 <- ggplot(data = turk0, mapping = aes(x = A, y = Gain)) +
    geom_point() +
    xlab("Metionina (% da dieta)") +
    ylab("Ganho de peso (g)")

p0 +
    stat_summary(geom = "line", fun.y = "mean") +
    geom_smooth()

#-----------------------------------------------------------------------
# Encontrar valores iniciais e compreender a função.

fit <- rp.nls(model = Gain ~ Int + (Asy - Int) * A/(Half + A),
              data = turk0,
              start = list(Int = c(600, 650),
                           Asy = c(750, 850),
                           Half = c(0, 0.2)),
              extra = function(Int, Asy, Half) {
                  abline(h = c(Asy, Int),
                         v = Half,
                         col = "green")
              },
              xlab = "Metionine",
              ylab = "Weight gain")

summary(fit)

#-----------------------------------------------------------------------
# Usando a nls() mas com outra parametrização.

n0 <- nls(Gain ~ f0 + f1 * (1 - exp(-k * A)),
          data = turk0,
          start = list(f0 = 600, f1 = 200, k = 10))
summary(n0)

# k é proporcional a taxa em f(0) pois f'(0) = k * f1.
D(expression(f0 + f1 * (1 - exp(-k * A))), "A")

plot(Gain ~ A,
     data = turk0,
     xlab = "Metionina (% da dieta)",
     ylab = "Gano de peso (g)")
with(as.list(coef(n0)),
     curve(f0 + f1 * (1 - exp(-k * x)),
           add = TRUE,
           col = 2,
           lwd = 3))

#-----------------------------------------------------------------------
# Existe falta de ajuste?

m0 <- lm(Gain ~ factor(A), data = turk0)

# Teste da soma de quadrados residual extra.
anova(n0, m0)

p0 +
    stat_function(fun = function(x) predict(n0, newdata = list(A = x)),
                  color = "cyan") +
    stat_summary(geom = "point", fun.y = "mean", color = "orange")

#-----------------------------------------------------------------------
# Inspeção nos resíduos.

plot(residuals(n0) ~ fitted(n0))
qqnorm(residuals(n0))

par(mfrow = c(2, 2))
plot(as.lm(n0))
layout(1)

# Coeficiente de determinação.
1 - deviance(n0)/deviance(lm(Gain ~ 1, turk0))
R2nls(n0)

# cor(fitted(n0), turk0$Gain)^2
cor(x = fitted(n0),
    y = get(as.character(n0$data))[, all.vars(formula(n0))[1]])^2

# A função `as.lm()` está simplesmente ajustando um modelo não linear.
# Isso é feito usando a matriz de derivadas parciais em relação aos
# parâmetros.

# Para obter a matriz de gradiente (derivadas parciais).
der <- deriv3(expr = formula(n0)[-2],
                 namevec = names(coef(n0)),
                 function.arg = function(A, f0, f1, k) { NULL })

c0 <- coef(n0)
m <- der(A = turk0$A,
         f0 = c0["f0"],
         f1 = c0["f1"],
         k = c0["k"])
X <- attr(m, "gradient")
tail(X)

# Ajuste de modelo linear usando a matriz X.
# m0 <- lm(Gain ~ 0 + X, data = turk0)
m0 <- lm(Gain ~ X[, -1], data = turk0)

# ATTENTION: o primeiro parâmetro e o intercepto, por isso a primeira
# coluna de X foi removida.

par(mfrow = c(2, 2))
plot(m0)
layout(1)

anova(m0)
summary(m0)

# ATTENTION: O último parâmetro não aparece linearmente no modelo.
cbind(coef(n0), coef(m0))

#-----------------------------------------------------------------------
# Mais um exemplo.

data(wtloss, package = "MASS")

ggplot(wtloss, aes(Days, Weight)) +
    geom_point()

fit <- rp.nls(model = Weight ~ f0 + f1 * 2^(-Days/K),
              data = wtloss,
              start = list(f0 = c(70, 140),
                           f1 = c(20, 200),
                           K = c(0, 250)))

summary(fit)

#-----------------------------------------------------------------------
# Fazendo ajuste com a `nls()`.

n0 <- nls(Weight ~ f0 + f1 * 2^(-Days/K),
          data = wtloss,
          start = list(f0 = 77, f1 = 107, K = 150))
summary(n0)

par(mfrow = c(2, 2))
plot(as.lm(n0))
layout(1)

plot(acf(residuals(n0)))
plot(pacf(residuals(n0)))
R2nls(n0)

#-----------------------------------------------------------------------
# Modelo hipsométrico.

# Importação dos dados.
dap <- read.table("http://www.leg.ufpr.br/~walmes/data/dap.txt",
                  header = TRUE)
names(dap) <- c("d", "h")

# Ordena e separa os casos completos.
dap <- dap[order(dap$d), ]
dapcc <- dap[complete.cases(dap), ]
rownames(dapcc) <- NULL

p0 <- ggplot(dapcc, aes(d, h)) +
    geom_point() +
    xlab("Diâmetro a altura do peito") +
    ylab("Altura da árvore")
p0

# Ajuste do modelo.
n0 <- nls(h ~ b0 * (1 - exp(-b1 * d))^b2,
          data = dapcc,
          start = list(b0 = 35, b1 = 0.1, b2 = 1.05))
summary(n0)

# Verificando o ajuste.
p0 +
    stat_function(fun = function(d) predict(n0, newdata = list(d = d)),
                  color = "purple")

# Para diagnósticos, transformar em modelo linear.
m0 <- as.lm(n0)
m0$call <- NULL # Elimina elemento para evitar prints verbosos.

par(mfrow = c(2, 2))
plot(m0)
layout(1)

im <- influence.measures(m0)
summary(im)

# str(im)
rem <- im$is.inf[, "dffit"]

# Exibe os pontos influentes.
p0 +
    geom_point(data = dapcc[rem, ],
               color = "red", pch = 1, cex = 3)

# Reajusta com remoção dos outliers.
n1 <- update(n0, data = dapcc[!rem, ])
summary(n1)

# Intervalos de confiança para os parâmetros.
confint.default(n1)
confint(n1)

# Ajuste de um modelo mais simples.
n2 <- nls(h ~ b0 * (1 - exp(-b1 * d)),
          data = dapcc[!rem, ],
          start = list(b0 = 35, b1 = 0.1))

summary(n2)
anova(n1, n2)

# Uma patologia bem comum.
cov2cor(vcov(n1))
cov2cor(vcov(n2))

#-----------------------------------------------------------------------
# Predições: cenas para os próximos capítulos.

grid <- data.frame(d = seq(min(dapcc$d),
                           max(dapcc$d),
                           length.out = 101))

a <- confbands(object = n1, newdata = grid[1])
colnames(a) <- paste(colnames(a), "n1", sep = "_")
grid <- cbind(grid, a)

a <- confbands(object = n2, newdata = grid[1])
colnames(a) <- paste(colnames(a), "n2", sep = "_")
grid <- cbind(grid, a)
str(grid)

ggplot(grid) +
    geom_point(data = dapcc[!rem, ], mapping = aes(d, h)) +
    xlab("Diâmetro a altura do peito") +
    ylab("Altura da árvore") +
    geom_line(aes(d, lwr_n1), col = "orange") +
    geom_line(aes(d, upr_n1), col = "orange") +
    geom_line(aes(d, lwr_n2), col = "purple") +
    geom_line(aes(d, upr_n2), col = "purple")

#-----------------------------------------------------------------------
