#-----------------------------------------------------------------------
#                                            Prof. Dr. Walmes M. Zeviani
#                                leg.ufpr.br/~walmes · github.com/walmes
#                                        walmes@ufpr.br · @walmeszeviani
#                      Laboratory of Statistics and Geoinformation (LEG)
#                Department of Statistics · Federal University of Paraná
#                                       2018-nov-21 · Curitiba/PR/Brazil
#-----------------------------------------------------------------------

#=======================================================================
# Inferência para funções de parâmetros.

#-----------------------------------------------------------------------
# Pacotes.

library(lattice)
library(latticeExtra)
library(rpanel)
library(nlme)
library(car)
library(wzRfun)
# library(plyr)

source("http://leg.ufpr.br/~walmes/cursoR/mrnl2013/fun%c3%a7%c3%b5es/filleddensityplot.R")
# source("../funções/filleddensityplot.R")

#-----------------------------------------------------------------------
# Curvas de retenção de água do solo.

cra <- read.table("http://www.leg.ufpr.br/~walmes/data/cra_manejo.txt",
                  header = TRUE, sep = "\t")

xtabs(~condi + prof, data = cra)

cras <- subset(cra, condi == "LVA3,5" & prof == 80)
cras$tens[cras$tens == 0] <- 0.1

cras <- aggregate(umid ~ posi + tens, data = cras, FUN = mean)
cras$ltens <- log(cras$tens)
str(cras)

p0 <- xyplot(umid ~ log(tens), groups = posi, data = cras, type = "p",
             ylab = "Conteúdo de água do solo",
             xlab = "Tensão matricial", auto.key = list(columns = 2))
p0

model <- umid ~ Ur + (Us - Ur)/(1 + exp(n * (al + ltens)))^(1 - 1/n)

start <- list(Ur = c(init = 0.2, from = 0, to = 0.5),
              Us = c(init = 0.6, from = 0.4, to = 0.8),
              al = c(1, -2, 4),
              n = c(1.5, 1, 4))

cra.fit <- rp.nls(model = model,
                  data = cras,
                  start = start,
                  subset = "posi")

dput(sapply(cra.fit, coef))

cols <- trellis.par.get()$superpose.line$col[1:2]
p0 +
    glayer({
        b <- coef(cra.fit[[group.number]])
        with(as.list(b),
             panel.curve(Ur + (Us - Ur)/(1 + exp(n * (al + x)))^(1 - 1/n),
                         add = TRUE, col = col.line))
    })

# Taxa no ponto de inflexão (S é de slope).
Scalc <- function(obj) {
    if (class(obj) == "nls") {
        with(as.list(coef(obj)), {
            -(Us - Ur) * n * (1 + 1/(1 - 1/n))^(-(1 - 1/n) - 1)
        })
    } else {
        names(obj) <- c("Ur", "Us", "al", "n")
        with(as.list(obj), {
            -(Us - Ur) * n * (1 + 1/(1 - 1/n))^(-(1 - 1/n) - 1)
        })
    }
}

sapply(cra.fit, Scalc)

# Inferência por simulação Monte Carlo. --------------------------------

S_dist <- lapply(cra.fit,
                 FUN = function(fit) {
                     m <- coef(fit)
                     V <- vcov(fit)
                     boot <- mvtnorm::rmvnorm(n = 1000,
                                              mean = m,
                                              sigma = V)
                     S <- apply(boot, MARGIN = 1, FUN = Scalc)
                     return(S)
                 })

str(S_dist)

densityplot(~EL + L, data = S_dist) +
    glayer({
        panel.abline(v = mean(x, na.rm = TRUE), col = col.line, lty = 2)
    })

densityplot(~EL + L, outer = TRUE, data = S_dist) +
    layer({
        qnt <- quantile(x, c(0.025, 0.975), na.rm = TRUE)
        panel.abline(v = mean(x, na.rm = TRUE), lty = 2)
        panel.abline(v = qnt, lty = 3)
    })

qqmath(~EL + L, outer = TRUE, data = S_dist, scales = "free")

#-----------------------------------------------------------------------
# Bootstrap.

n0 <- cra.fit[[1]]
n1 <- cra.fit[[2]]

# Número de reamostragens ou números simulados.
B <- 500

# Resíduos com os dois sinais (+ e -).
res <- c(-residuals(n0), residuals(n0))

# Estimativa da variância dos resíduos.
s <- summary(n0)$sigma

R <- data.frame(
    resid_rea = sample(res, B * length(res)/2, replace = TRUE),
    resid_sim = rnorm(B * length(res)/2, 0, sd = s))
R$b <- rep(1:B, each = length(res)/2)
R$tens <- unique(cras$tens)
R$ltens <- log(R$tens)
R$f <- fitted(n0)

# Produz novas observações ou pseudo observações.
R$y_rea <- R$f + R$resid_rea
R$y_sim <- R$f + R$resid_sim

xyplot(y_rea + y_sim ~ ltens,
       outer = TRUE,
       groups = b,
       data = R,
       jitter.x = TRUE)

str(R)

n0_rea <- nlsList(y_rea ~ Ur + (Us - Ur)/(1 + exp(n * (al + ltens)))^(1 - 1/n) | b,
                  data = R,
                  start = as.list(coef(n0)))

n0_sim <- nlsList(y_sim ~ Ur + (Us - Ur)/(1 + exp(n * (al + ltens)))^(1 - 1/n) | b,
                  data = R,
                  start = as.list(coef(n0)))

boot0_npar <- coef(n0_rea)
boot0_par <- coef(n0_sim)

# Determina a distribuição amostral de S.
S0_npar <- apply(boot0_npar, 1, Scalc)
S0_par <- apply(boot0_par, 1, Scalc)

densityplot(~S0_npar + S0_par) +
    glayer({
        qnt <- quantile(x, c(0.025, 0.975), na.rm = TRUE)
        panel.abline(v = mean(x, na.rm = TRUE), lty = 2, col = col.line)
        panel.abline(v = qnt, lty = 3, col = col.line)
    })

#-----------------------------------------------------------------------
# Método delta.

g <- expression(-(Us - Ur) * n * (1 + 1/(1 - 1/n))^(-(1 - 1/n) - 1))
pars <- c("Ur", "Us", "al", "n")

# Derivadas parciais com relação aos parâmetros.
lapply(pars, D, expr = g)

# Para termos o gradiente da função.
grad <- deriv3(~ -(Us - Ur) * n * (1 + 1/(1 - 1/n))^(-(1 - 1/n) - 1),
               pars,
               function(Ur, Us, al, n) {
                   NULL
               })

# Vetor gradiente avaliação em \hat{\theta}.
gr <- do.call(grad, as.list(coef(n0)))
gr <- attr(gr, "gradient")
gr

Scalc(n0)

# Erro padrão de \hat{\vartheta} = g(\hat{\theta}).
sqrt(gr %*% vcov(n0) %*% t(gr))

# Expressão da função dos parâmetros.
expr_string <- "-(Us - Ur) * n * (1 + 1/(1 - 1/n))^(-(1 - 1/n) - 1)"

# Método delta.
dm0 <- deltaMethod(n0, g = expr_string)
dm0

# Método delta.
dm1 <- deltaMethod(n1, g = expr_string)
dm1

#-----------------------------------------------------------------------
# Uma reparametrização do van Genuchten. Produto da minha tese de
# Doutorado.

# Fórmula do modelo reparametrizado para inferência no ponto de
# inflexão.
model <- umid ~ Ur - Si * (1 + 1/(1 - 1/n))^((1 - 1/n) + 1)/(n * (1 + exp(n * (ltens - psii))/(1 - 1/n))^(1 - 1/n))

# Lista de limites para os parâmetros.
start <- list(Ur = c(init = 0.2, from = 0, to = 0.5),
              Si = c(-0.08, -0.01),
              psii = c(0, 3),
              n = c(1.5, 1, 4))

# Ajuste com GUI.
cra.fit <- rp.nls(model = model,
                  data = cras,
                  start = start,
                  subset = "posi")

# Resultado do ajuste da reparametrização do modelo.
sapply(cra.fit, coef)
lapply(cra.fit, summary)
lapply(cra.fit, confint.default)

# O que foi obtido com o método delta.
dm0
dm1

#-----------------------------------------------------------------------
# Ajuste dos modelos simultaneamente.

theta <- sapply(cra.fit, coef)
theta

start <- theta %*% cbind(c(1, 0), c(-1, 1))
start

# start <- theta %*% cbind(c(1, 0), c(0, 0))
start <- c(t(start))

# IMPORTANT: muita atenção com a ordem dos parâmetros em `params =` e a
# ordem dos valores iniciais. Além disso, esteja ciente da
# parametrização/codificação dos efeitos.

# Ajuste do modelo conjunto.
fit <- gnls(model,
            data = cras,
            params = Ur + Si + n + psii ~ posi,
            start = start,
            verbose = TRUE)

# Resumo do modelo.
summary(fit)

# Testes de hipótese.
anova(fit, type = "sequential")
anova(fit, type = "marginal")

#-----------------------------------------------------------------------
