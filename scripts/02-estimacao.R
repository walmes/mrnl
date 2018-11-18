#-----------------------------------------------------------------------
#                                            Prof. Dr. Walmes M. Zeviani
#                                leg.ufpr.br/~walmes · github.com/walmes
#                                        walmes@ufpr.br · @walmeszeviani
#                      Laboratory of Statistics and Geoinformation (LEG)
#                Department of Statistics · Federal University of Paraná
#                                       2018-nov-18 · Curitiba/PR/Brazil
#-----------------------------------------------------------------------

#=======================================================================
# Estimação e inferência em modelos não lineares.

#-----------------------------------------------------------------------
# Pacotes.

library(bbmle)
library(rootSolve)

source("_funs.R")

#-----------------------------------------------------------------------
# Dados artificiais simulados do modelo Michaelis-Menten.

set.seed(111)
x <- 1:10
y <- 10 * x/(3 + x) + rnorm(x, 0, 0.3)
plot(y ~ x)

# Derivadas da função em relação à \theta.
f <- expression(A * x/(V + x))

# Derivadas parciais de primeira ordem
D(f, "A")
D(f, "V")

# Derivadas parciais de segunda ordem.
D(D(f, "A"), "A")
D(D(f, "A"), "V")
D(D(f, "V"), "A")
D(D(f, "V"), "V")

# `deriv3()` retorna \eta(.), \eta'(.) e \eta''(.).
d3 <- deriv3(expr = ~A * x/(V + x),
             namevec = c("A", "V"),
             function.arg = function(x, A, V) { NULL })

# Uso.
str(d3(x, A = 10, V = 3))

#-----------------------------------------------------------------------
# Estimação de \theta com método numérico.

# Valor inicial para A e V.
theta0 <- c(7, 5)

sqr0 <- crossprod(y)      # Soma de quadrados total (não corrigida).
dif <- sqr0               # Valor inicial para diferença.
i <- 0                    # Contador de iterações.
maxiter <- 50             # Número máximo de iterações.
tol <- 1e-8;              # DIferença mínima para convergência.
while (i <= 50 & dif > tol) {
    #----------------------------------------
    txt <- paste(formatC(theta0, digits = 6, format = "f"),
                         collapse = " ")
    txt <- sprintf("i: %d\t SQR: %0.8f\t theta: %s\n",
                   i, sqr0, txt)
    cat(txt)
    #----------------------------------------
    f012 <- d3(x, A = theta0[1], V = theta0[2])
    f0 <- c(f012)                # f(.)
    f1 <- attr(f012, "gradient") # f'(.)
    res <- (y - f0)              # Resíduos.
    sqr1 <- crossprod(res)       # SQR.
    dif <- abs(sqr0 - sqr1)      # Redução em SQR.
    theta1 <- theta0 + solve(crossprod(f1), crossprod(f1, res))
    theta0 <- c(theta1)
    sqr0 <- sqr1
    i <- i + 1
    #----------------------------------------
}

theta0

#-----------------------------------------------------------------------
# Conferindo com a nls().

n0 <- nls(y ~ A * x/(V + x),
          start = c(A = 7, V = 5),
          trace = TRUE)
summary(n0)

#-----------------------------------------------------------------------
# Visualiando o ajuste.

plot(y ~ x)
with(as.list(coef(n0)),
     curve(A * x/(V + x),
           col = 2,
           add = TRUE))

#-----------------------------------------------------------------------
# Quadro de partição de soma de quadrados.

sqtot <- sum(y^2)
sqtcm <- sum((y - mean(y))^2)
sqres <- sum((y - f0)^2)
sqreg <- sqtot - sqres
p <- length(theta0)
r <- length(y) - length(theta0)

data.frame(GL = c(reg = p, res = r),
           SQ = c(sqreg, sqres),
           QM = c(sqreg/p, sqres/r),
           F = c((sqreg/p)/(sqres/r), NA))

# Mais simles! A fórmula é y ~ 0 por que o modelo não tem intercepto.
anova(n0, lm(y ~ 0))

#-----------------------------------------------------------------------
# Obtendo estimativa da variância resídual.

s2 <- sum((y - f0)^2)/(length(x) - length(theta0))
s2

# Extração do objeto contendo o ajuste.
summary(n0)$sigma^2

#-----------------------------------------------------------------------
# Obtendo erro-padrão da estimativa.

# Mínimos quadrados ordinários (OLS).
v <- s2 * solve(t(f1) %*% f1)
v
sqrt(diag(v))

# Extrator da matriz de covariâncias das estimativas do modelo ajustado.
vcov(n0)

# Usando W = I * \sigma^2: Forma geral (OLS, WLS, GLS).
v <- solve(t(f1) %*% solve(diag(s2, length(y))) %*% f1)
sdt <- sqrt(diag(v))

#-----------------------------------------------------------------------
# Intervalo de confiança assintótico para os parâmetros.

qnt <- c(-1, 1) * qnorm(0.975)
theta0[1] + qnt * sdt[1]
theta0[2] + qnt * sdt[2]

# Método disponível para o objeto.
confint.default(n0)

# ATTENTION: está sendo usando o quantil da Normal e não o da t.

#-----------------------------------------------------------------------
# Quadro de somas de quadrados.

anova(n0, lm(y ~ 0))[2:1, ]

#-----------------------------------------------------------------------
# Mensagens de erro durante convergência com a `nls()`.

# 1) Matriz de gradiente singular (devido A = 0) nos valores iniciais.
n0 <- nls(y ~ A * x/(V + x),
          start = c(A = 0, V = 5),
          trace = TRUE)

f012 <- d3(x, A = 0, V = 5)
f1 <- attr(f012, "gradient")
f1

crossprod(f1)

# TIP: estude os espaço paramétrico dos parâmetros e as propriedades
# matemáticas do modelo não linear.

# 2) Indeterminações matemáticas como  NaN e Inf.
# (devido a V = -1 e ter x = 1, então V + x = 0 no denominador)
n0 <- nls(y ~ A * x/(V + x),
          start = c(A = 7, V = -1),
          trace = TRUE)

f012 <- d3(x, A = 7, V = -1)
f1 <- attr(f012, "gradient")
f1

# TIP: estude os espaço paramétrico dos parâmetros e as propriedades
# matemáticas do modelo não linear e, além disso, verifique se o
# conjunto de valores das preditoras está no domínio correto.

# 3) Convergência para mínimo local.
n0 <- nls(y ~ A * x/(V + x),
          start = c(A = 7, V = -1.1),
          trace = TRUE)

# TIP: conhecimento preliminar sobre o modelo e os dados, uma boa
# análise exploratória ou métodos gráficos para obter chutes iniciais
# podem ser úteis.

# 4) Gradiente singular por valores iniciais ruins.
n0 <- nls(y ~ A * x/(V + x),
          start = c(A = -17, V = 100),
          trace = TRUE)

f012 <- d3(x, A = -1044.548, V = 7358.790)
f1 <- attr(f012, "gradient")
f1

# TIP: o segredo é estar com o modelo correto e fornecer bons chutes
# iniciais.

# 5) Modelo completamente não identificável por parametrização.
n0 <- nls(y ~ A * B * x/(V + x),
          start = c(A = 5, B = 2, V = 1),
          trace = TRUE)

n0 <- nls(y ~ A^C * x/(V + x),
          start = c(A = 3, C = 0.1, V = 1),
          trace = TRUE)

# TIP: existem infinitas combinações de A * B ou A^C que poduzem o mesmo
# valor. O que é estimável é U = A * B = A^C, mas isoladamete não é
# estimável. O modelo deve ser simplificado na expressão.

# Dessa forma, fator de passos reduzido, pouca curvatura ao redor das
# estimativas, são fatores relacionados à baixa identificabilidade do
# modelo, correlações altas entre estimativas, etc.

# 6) Modelo incorreto para os dados.
d <- list(y = c(0.9, 2, 3, 4, 5, 6, 6.9, 8.1, 9, 10.1),
          x = c(2.6, 3.9, 4.9, 5, 6.2, 6.7, 6.6, 7, 7.2, 7.5))
plot(y ~ x, data = d)

n0 <- nls(y ~ A * x/(V + x),
          data = d,
          start = c(A = 10, V = 1),
          trace = TRUE)

# TIP: os dados exibem um sinal de concavidade para cima e o modelo só
# acomoda concavidade para baixo. Os parâmetros são mandados para o
# infinito na tentativa de f(.) ficar uma reta, mas o gradiente fica
# singular.  Estude a função para saber se ela acomoda o sinal dos
# dados.

#-----------------------------------------------------------------------------
# Valores iniciais e trajetória sobre a superfície da função objetivo.

# Função objetivo RSS - residual sum of squares.
rss <- function(A, V, x, y) {
    sum((y - (A * x)/(V + x))^2)
}
RSS <- Vectorize(FUN = rss, vectorize.args = c("A", "V"))

# Grid dos valores dos parâmetros para plotar a função.
A.grid <- seq(0, 40, length.out = 100)
V.grid <- seq(0, 20, length.out = 100)
rss.surf <- outer(A.grid, V.grid, RSS, x = x, y = y)

# Contornos.
contour(A.grid, V.grid, rss.surf,
        levels = (1:35)^2, xlab = "A", ylab = "V",
        col = "gray70", main = "Superfície da função RSS")

persp(A.grid, V.grid, log(rss.surf), col = "red")

#-----------------------------------------------------------------------
# Visualização interativa com RGL.

library(rgl)

# z <- log(rss.surf)
z <- rss.surf
ylim <- range(z)
ylen <- ylim[2] - ylim[1] + 1
colorlut <- heat.colors(ylen)
col <- colorlut[z - ylim[1] + 1]

persp3d(A.grid, V.grid, z, color = col, back = "lines")

#-----------------------------------------------------------------------
# Diferentes rotas até a convergência.

# Lista de chutes iniciais para mostrar diferentes trajetórias.
start.list <- list(s1 = c(A = 0.1, V = 0.1),
                   s2 = c(A = 40, V = 20),
                   s3 = c(A = 35, V = 2.5),
                   s4 = c(A = 18, V = 18))

# oldpar <- par()
# par(mfrow = c(2, 2), mar = c(4, 4, 1, 1))

for (lis in 1:4) {
    contour(A.grid, V.grid, rss.surf,
            levels = (seq(1, 35, 2))^2,
            xlab = "A", ylab = "V", col = "gray70")
    sink("trace.txt")
    n0 <- nls(y ~ A * x/(V + x),
              start = start.list[[lis]],
              trace = TRUE)
    sink()
    trace <- read.table("trace.txt")
    for (i in seq(nrow(trace) - 1)) {
        arrows(trace[i, "V3"],
               trace[i, "V4"],
               trace[i + 1, "V3"],
               trace[i + 1, "V4"],
               col = 2,
               length = 0.1)
        abline(v = trace[i + 1, "V3"],
               h = trace[i + 1, "V4"],
               col = "orange",
               lty = 3)
        Sys.sleep(1)
        print(c(i, trace[i + 1, "V3"], trace[i + 1, "V4"]))
    }
}

sink()

n0 <- nls(y ~ A * x/(V + x),
          start = c(A = 7, V = 5))
summary(n0)

#-----------------------------------------------------------------------
# Agora é conveniente interpretar a matriz de covariância das
# estimativas!

vcov(n0)
cov2cor(vcov(n0))

# S(\theta) = \sum_{i = 1}^{n} (y_i - f(x_i, \theta))^2.
#
# TODO: fazer S*(\theta) pela aproximação em série de Taylor de
# S(\theta) ao redor de \hat{\theta}.

# Derivadas parciais de primeira ordem de S(\theta) em \theta.
D(expression((y - (A * x/(V + x)))^2), "A")
D(expression((y - (A * x/(V + x)))^2), "V")

# Vetor gradiente.
G <- with(as.list(coef(n0)), {
    # ATTENTION: as expressões pode ser simplicada para redução de custo
    # computacional.
    s <- cbind(-(2 * (x/(V + x) * (y - (A * x/(V + x))))),
          2 * (A * x/(V + x)^2 * (y - (A * x/(V + x)))))
    matrix(colSums(s), ncol = ncol(s))
})
G

# IMPORTANT! G não é necessário porque G vale 0 em \hat{\theta}.

# Derivadas parciais de segunda ordem de S(\theta) em \theta.
D(D(expression((y - (A * x/(V + x)))^2), "A"), "A")
D(D(expression((y - (A * x/(V + x)))^2), "V"), "V")
D(D(expression((y - (A * x/(V + x)))^2), "A"), "V")

# Matriz hessiana.
H <- with(as.list(coef(n0)), {
    # ATTENTION: as expressões pode ser simplicada para redução de custo
    # computacional.
    mAA <- sum(2 * (x/(V + x) * (x/(V + x))))
    mVV <- sum(2 * (A * x/(V + x)^2 * (A * x/(V + x)^2) - A * x * (2 * (V + x))/((V + x)^2)^2 * (y - (A * x/(V + x)))))
    mAV <- sum(-(2 * (x/(V + x) * (A * x/(V + x)^2) - x/(V + x)^2 * (y - (A * x/(V + x))))))
    matrix(c(mAA, mAV, mAV, mVV), 2, 2)
})
H

# Função que retorna a aproximação da SQR.
rss_a <- function(A, V, deviance = deviance(n0), b = coef(n0)) {
    # u <- rbind(A, V) - rbind(coef(n0))
    u <- c(A, V) - b
    u <- matrix(u, nrow = 2)
    deviance + G %*% u + t(u) %*% H %*% u
}
RSS_approx <- Vectorize(FUN = rss_a, vectorize.args = c("A", "V"))

# Intervalo de confiança assintótico.
confint.default(n0)

# Região de confiança para os parâmetros.
np <- length(coef(n0))
nl <- length(y)
conf_reg <- ellipse::ellipse(n0, t = sqrt(np * qf(0.95, np, nl)))
rss_cont <- c(rss_a(A = conf_reg[1, 1], V = conf_reg[1, 2],
                    deviance = deviance(n0), b = coef(n0)))

# Determinação da superfície.
er <- extendrange(conf_reg[, 1], f = 0.1)
A.grid <- seq(er[1], er[2], length.out = 51)
er <- extendrange(conf_reg[, 2], f = 0.1)
V.grid <- seq(er[1], er[2], length.out = 51)
rss.surf <- outer(A.grid, V.grid, RSS, x = x, y = y)
rss_a.surf <- outer(A.grid,
                    V.grid,
                    FUN = RSS_approx,
                    deviance = deviance(n0),
                    b = coef(n0))
rss_a.surf[rss_a.surf > max(rss.surf)] <- NA

# Gráfico.
persp3d(A.grid, V.grid, rss.surf,
        col = "orange", back = "lines", front = "fill")
persp3d(A.grid, V.grid, rss_a.surf,
        col = "cyan", back = "lines", front = "fill", add = TRUE)
segments3d(rbind(c(coef(n0), deviance(n0)),
                 c(coef(n0), max(rss.surf))), col = "red", lwd = 3)
points3d(cbind(conf_reg, rss_cont), col = "green")
lines3d(cbind(conf_reg, rss_cont), col = "green")

# Ponto de corte na deviance para região de confiança.
np <- length(coef(n0))
nl <- length(y)
dev95 <- deviance(n0) *
    (1 + (np/(nl - np)) *
     qf(p = 0.95, df1 = np, df2 = (nl - np), lower.tail = TRUE))

# Adiciona NA para dar ênfase.
rss.surf2 <- rss.surf
rss.surf2[rss.surf2 > 2 * dev95] <- NA

persp3d(A.grid, V.grid, rss.surf2, col = "purple",
        back = "lines", front = "fill",
        zlim = extendrange(rss.surf2, f = 1))
persp3d(A.grid, V.grid, rep(dev95, prod(dim(rss.surf))),
        col = "yellow", add = TRUE)
lines3d(cbind(conf_reg, rss_cont), col = "green")
segments3d(rbind(c(coef(n0), deviance(n0)),
                 c(coef(n0), max(rss.surf2, na.rm = TRUE))),
           col = "red", lwd = 3)

#-----------------------------------------------------------------------------
# Valores iniciais inadequados podem levar para mínimos locais.  Clique
# sobre a superfície de contornos para começar a iteração a partir da
# escolha.

A.grid <- seq(-5, 20, l = 100)
V.grid <- seq(-3, 10, l = 100)
rss.surf <- outer(A.grid, V.grid, RSS, y = y, x = x)

contour(A.grid, V.grid, rss.surf,
        levels = (1:35)^2, xlab = "A",
        ylab = "V", col = "gray70",
        main = "Superfície da função RSS")
start.click <- locator(n = 1)
names(start.click) <- c("A", "V")
sink("trace.txt")
n0 <- nls(y ~ A * x/(V + x),
          start = start.click,
          trace = TRUE)
sink()
trace <- read.table("trace.txt")
for (i in seq(nrow(trace) - 1)) {
    arrows(trace[i, "V3"],
           trace[i, "V4"],
           trace[i + 1, "V3"],
           trace[i + 1, "V4"],
           col = 2,
           length = 0.1)
    Sys.sleep(0.3)
}

#-----------------------------------------------------------------------
# Aspectos da inferência baseada em verossimilhança.

set.seed(321)
x <- 1:10
# y <- 1 - exp(-log(2) * x/2) + rnorm(x, 0, 0.05)
y <- 1 - 2^(-x/2) + rnorm(x, 0, 0.05)
plot(y ~ x)

# Parametrizações do modelo monomolecular (assíntota fixada em 1):
#  * com parâmetro theta_c (taxa na origem): 1 - exp(-thc * x)
#  * com parâmetro theta_v (tempo de meia vida): 1 - 2^(-x/thv).

#-----------------------------------------------------------------------
# Funções de log-verossimilhança.

llc <- function(thc, x, y) {
    m <- 1 - exp(-thc * x)
    s <- sqrt(sum((y - m)^2)/length(y))
    sum(dnorm(y, m, s, log = TRUE))
}

llv <- function(thv, x, y) {
    m <- 1 - 2^(-x/thv)
    s <- sqrt(sum((y - m)^2)/length(y))
    sum(dnorm(y, m, s, log = TRUE))
}

#-----------------------------------------------------------------------
# Visualizando a função de ll.

len <- 100

thv.g <- seq(1.7, 2.3, l = 40)
llv.g <- sapply(thv.g, llv, x = x, y = y)

thc.g <- seq(log(2)/max(thv.g), log(2)/min(thv.g), l = 40)
llc.g <- sapply(thc.g, llc, x = x, y = y)

par(mfrow = c(1, 2))
plot(llc.g ~ thc.g, type = "l")
plot(llv.g ~ thv.g, type = "l")
layout(1)

#-----------------------------------------------------------------------------
# Estimando os parâmetros por máxima verossimilhança.

opc <- optim(c(log(2)/2), llc, x = x, y = y, method = "BFGS",
             control = list(fnscale = -1), hessian = TRUE)
opv <- optim(c(2), llv, x = x, y = y, method = "BFGS",
             control = list(fnscale = -1), hessian = TRUE)

(thc <- opc$par)
(thv <- opv$par)
log(2)/thv

(llthc <- opc$value)
(llthv <- opv$value)

(vthc <- -1/opc$hessian)
(vthv <- -1/opv$hessian)
qchi <- qchisq(0.95, df = 1)
qnor <- qnorm(0.975)

#-----------------------------------------------------------------------
# Visualizando a log-verossimilhança e a deviance.

par(mfrow = c(2, 2))
# 1 -----------------------------------------
plot(llc.g ~ thc.g,
     type = "l",
     ylab = "log-verossimilhança",
     xlab = expression(th[c]))
abline(v = thc, h = llthc - c(0, qchi), lty = 2)
# 2 -----------------------------------------
devc <- -2 * (llc.g - llthc)
plot(devc ~ thc.g,
     type = "l",
     ylab = "Deviance",
     xlab = expression(th[c]))
abline(v = thc, h = c(0, qchi), lty = 2)
curve((x - thc)^2/c(vthc), add = TRUE, col = 2)
# 3 -----------------------------------------
plot(llv.g ~ thv.g,
     type = "l",
     ylab = "log-verossimilhança",
     xlab = expression(th[v]))
abline(v = thv, h = llthv - c(0, qchi), lty = 2)
# 4 -----------------------------------------
devv <- -2 * (llv.g - llthv)
plot(devv ~ thv.g,
     type = "l",
     ylab = "Deviance",
     xlab = expression(th[v]))
abline(v = thv, h = c(0, qchi), lty = 2)
curve((x - thv)^2/c(vthv), add = TRUE, col = 2)
layout(1)

#-----------------------------------------------------------------------
# Funções de verossimilhança na escala da deviance para obter IC.

dev.c <- function(thc0, llthc, qchi, x, y) {
    dev <- sapply(thc0, function(thc00) {
        m0 <- 1 - exp(-thc00 * x)
        s0 <- sqrt(sum((y - m0)^2)/length(y))
        ll0 <- sum(dnorm(y, m0, s0, log = TRUE))
        dev <- -2 * (ll0 - llthc)
    })
    return(dev - qchi)
}

dev.v <- function(thv0, llthv, qchi, x, y) {
    dev <- sapply(thv0, function(thv00) {
        m0 <- 1 - exp(-log(2) * x/thv00)
        s0 <- sqrt(sum((y - m0)^2)/length(y))
        ll0 <- sum(dnorm(y, m0, s0, log = TRUE))
        dev <- -2 * (ll0 - llthv)
    })
    return(dev - qchi)
}

#-----------------------------------------------------------------------
# Obtendo os IC baseados na deviance e na aproximação quadrática.

thc.lim <- uniroot.all(dev.c, interval = c(0.3, 0.4), llthc = llthc,
                       qchi = qchi, x = x, y = y)
thc.lim

thv.lim <- uniroot.all(dev.v, interval = c(1, 3), llthv = llthv,
                       qchi = qchi, x = x, y = y)
thv.lim

thc.piv <- thc + c(-1, 1) * qnor * sqrt(c(vthc))
thv.piv <- thv + c(-1, 1) * qnor * sqrt(c(vthv))

#-----------------------------------------------------------------------
# Verificando as diferenças com relação aos extremos do IC.

par(mfrow = c(1, 2))
# 1
devv <- -2 * (llv.g - llthv)
plot(devv ~ thv.g, type = "l",
     ylab = "Deviance",
     xlab = expression(th[v]))
abline(v = thv, h = c(0, qchi), lty = 2)
curve(((x - thv)^2/c(vthv)), add = TRUE, col = 2)
abline(v = thv.lim, lty = 2)
abline(v = thv.piv, lty = 2, col = 2)
# 2
devc <- -2 * (llc.g - llthc)
plot(devc ~ thc.g, type = "l",
     ylab = "Deviance", xlab = expression(th[c]))
abline(v = thc, h = c(0, qchi), lty = 2)
curve(((x - thc)^2/c(vthc)), add = TRUE, col = 2)
abline(v = thc.lim, lty = 2)
abline(v = thc.piv, lty = 2, col = 2)
layout(1)

#-----------------------------------------------------------------------
# Intervalos de confiança baseados na deviance e na aproximação
# quadrática.

n0 <- update(n0, trace = FALSE)

confint(n0)
confint.default(n0)

#-----------------------------------------------------------------------
