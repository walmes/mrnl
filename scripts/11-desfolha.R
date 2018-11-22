#-----------------------------------------------------------------------
#                                            Prof. Dr. Walmes M. Zeviani
#                                leg.ufpr.br/~walmes · github.com/walmes
#                                        walmes@ufpr.br · @walmeszeviani
#                      Laboratory of Statistics and Geoinformation (LEG)
#                Department of Statistics · Federal University of Paraná
#                                       2018-nov-21 · Curitiba/PR/Brazil
#-----------------------------------------------------------------------

#=======================================================================
# Ajuste de curvas com estrutura experimental.

#-----------------------------------------------------------------------
# Parâmetros.

library(lattice)
library(latticeExtra)
library(nlme)
library(car)
library(nlstools)
library(rootSolve)

#-----------------------------------------------------------------------
# Experimento dos capulhos.

da <- read.table("http://www.leg.ufpr.br/~walmes/data/algodão.txt",
                 header = TRUE, sep = "\t", encoding = "latin1")
da$desf <- da$desf/100
da$estag <- factor(da$estag, labels = c("Vegetativo", "Botão floral",
                                        "Florescimento", "Maçã",
                                        "Capulho"))
str(da)

# O valor de desfolha não pode ser 0 exato porque dá problema numérico
# na exponenciação com expoentes negativos. Então a dose zero será
# deslocada para direita com adição de um pequeno valor para evitar
# problemas.

0^-1     # Problema.
0.001^-1 # Solução.

# Desloca o zero para evitar problemas ao avaliar potência negativa.
da$desf[da$desf == 0] <- 0.001

xtabs(~estag + desf, da)

p0 <- xyplot(pcapu ~ desf | estag, data = da, layout = c(5, 1),
             xlab = "Nível de desfolha artificial",
             ylab = "Peso de capulhos produzidos (g)")
p0

p0 +
    layer({
        panel.smoother(method = "lm",
                       form = y ~ poly(x, degree = 2),
                       se = TRUE)
    })

# Modelos que serão ajustados.
models <- list("original" = pcapu ~ f0 - f1 * (desf + 0.02)^exp(C),
               "reparam" = pcapu ~ f0 - f1 * (desf + 0.02)^((log(5) - log(f1))/log(xde)))

db <- groupedData(formula = pcapu ~ desf | estag,
                  data = da,
                  order.groups = FALSE, )

n0 <- nlsList(model = models[[1]],
              data = db,
              start = list(f0 = 30, f1 = 10, C = 0))
coef(n0)

n1 <- nlsList(model = models[[2]],
              data = db,
              start = list(f0 = 30, f1 = 10, xde = 0.65))
coef(n1)

# Predição.
pred <- with(da,
             expand.grid(desf = seq(min(desf), max(desf), length.out = 21),
                         estag = levels(estag)))
pred$y <- predict(n1, newdata = pred)

# Sobresição do ajuste.
p0 + as.layer(xyplot(y ~ desf | estag,
                     data = pred,
                     type = "l",
                     col = 1))

#-----------------------------------------------------------------------
# Qual parametrização escolher? Verificar medidas de curvatura.

# Ajuste individual com `nls()` dentro de `by()`.

d3 <- list(deriv3(models[[1]][-2], c("f0", "f1", "C"),
                  function(desf, f0, f1, C) NULL),
           deriv3(models[[2]][-2], c("f0", "f1", "xde"),
                  function(desf, f0, f1, xde) NULL))

fits <- by(data = da,
           INDICES = da$estag,
           FUN = function(s) {
               n0 <- nls(pcapu ~ d3[[1]](desf, f0, f1, C),
                         data = s,
                         start = list(f0 = 30, f1 = 10, C = 0))
               n1 <- nls(pcapu ~ d3[[2]](desf, f0, f1, xde),
                         data = s,
                         start = list(f0 = 30, f1 = 10, xde = 0.65))
               list("orig" = n0, "repa" = n1)
           })

# Cada elemento é um estágio fenológico.
length(fits)

# Cada elemento é um modelo ajustado a cada estágio.
length(fits[[1]])

# Desaninhar a lista (deixar rasa).
fits <- unlist(fits, recursive = FALSE)
length(fits)
names(fits)

# Covariâncias.
lapply(fits, function(.) cov2cor(vcov(.)))

# Medidas de curvatura.
lapply(fits, MASS::rms.curv)

names(fits)

# Região de confiança para os parâmetros via aceitação-rejeição Monte
# Carlo.
i <- levels(da$estag)[4]
s <- subset(da, estag == i)
c0 <- nlsConfRegions(fits[[paste0(i, ".orig")]], exp = 3.5)
c1 <- nlsConfRegions(fits[[paste0(i, ".repa")]], exp = 3.5)

x11()
plot(c0)
x11()
plot(c1)

i <- levels(da$estag)[3]
s <- subset(da, estag == i)

r0 <- nlsContourRSS(fits[[paste0(i, ".orig")]])
r1 <- nlsContourRSS(fits[[paste0(i, ".repa")]])

x11()
plot(r0)
x11()
plot(r1)

# Lista apenas com os modelos reparametrizados.
# fit_repa <- rlist::list.match(fits, "repa")
fit_repa <- fits[seq(fits) %% 2 == 0]
names(fit_repa)

# Preparando os valores iniciais conforme fórmula `. ~ estag`, ou seja,
# primeiro nível é a categoria de referência, os demais são contrastes.
start_0 <- sapply(fit_repa, coef)
start_0

# Funções lineares para passar de `~0 + estag` para `~estag`.
k <- diag(ncol(start_0))
k[1, 2:ncol(k)] <- -1

start <- start_0 %*% k
start

#-----------------------------------------------------------------------
# Ajuste do modelo conjunto considerando estrutura experimental de
# efeitos.

nn0 <- gnls(pcapu ~ f0 - f1 * desf^((log(5) - log(f1))/log(xde)),
            data = da,
            params = f0 + f1 + xde ~ estag,
            start = c(t(start)))
coef(nn0)

summary(nn0)

# Teste de hipótese para o efeito dos termos experimentais nos
# parâmetros do modelo não linear.
anova(nn0, type = "marginal")

# Estimativas de `~estag` para `~0 + estag`
ctr <- cbind(1, contrasts(da$estag))
ctr <- kronecker(diag(3), ctr)
th <- ctr %*% coef(nn0)

matrix(th, ncol = 3)
sapply(fit_repa, coef)

# Erro padrão de cada estimativa.
vcov <- ctr %*% vcov(nn0) %*% t(ctr)
ep <- sqrt(diag(vcov))
ep

ztable <- data.frame(Estimate = th,
                     Std.Err = ep,
                     tvalue = th/ep,
                     pvalue = 2 * pnorm(-abs(th/ep)))
ztable$star <- cut(ztable$pvalue,
                   c(0, 0.01, 0.05, 0.1, 1),
                   labels = c("**", "*", ".", ""))
ztable$ics <- sweep(
    outer(X = ep,
          Y = c(lwr = -1, upr = 1) * qnorm(0.975),
          FUN = "*"), MARGIN = 1, STATS = th, FUN = "+")
ztable[, -(3:4)]

# Ajuste do modelo com estimativas por nível de estag, i.e. `~0 + estag`.
nn1 <- gnls(pcapu ~ f0 - f1 * desf^((log(5) - log(f1))/log(xde)),
            data = da,
            params = f0 + f1 + xde ~ -1 + estag,
            start = c(t(start_0)))
coef(nn1)

cbind(coef(nn1), th, c(t(start_0)))

pred <- with(da,
             expand.grid(desf = seq(min(desf),
                                    max(desf) - 0.001,
                                    length.out = 31),
                         estag = levels(estag),
                         KEEP.OUT.ATTRS = FALSE))
pred$pcapu <- predict(nn1, newdata = pred)
str(pred)

# Cria lista slot por nível de estag contendo vetor de parâmetros.
th <- coef(nn1)
u <- gsub("^.*\\.", "", names(th))
ths <- split(th, f = u)
names(ths)

# ATTENTION: o split não preserva ordem dos fatores. Cuide para
# restabeler a ordem dos itens da lista. Esse BUG me tirou 2 horas.

# Para ter a ordem dos itens da lista conforme ordem dos níveis do
# fator.
ths <- ths[unique(u)]
names(ths)

# Simplifica os nomes.
th_names <- c("f0", "f1", "xde")
ths <- lapply(ths, FUN = "names<-", th_names)
ths

# Para usar derivadas analíticas.
model <- deriv3(expr = formula(nn1)[-2],
                namevec = th_names,
                function.arg = function(desf, f0, f1, xde) { NULL })

# Para usar derivadas numéricas.
model <- function(theta, desf) {
    with(as.list(theta),
         f0 - f1 * desf^((log(5) - log(f1))/log(xde)))
}
# model(ths[[1]], desf = 0.1)
# gradient(f = model, ths[[1]], desf = 0.1)

i <- 0
f1 <- lapply(ths,
             FUN = function(theta) {
                 i <<- i + 1
                 desf <- pred$desf[pred$estag == levels(pred$estag)[i]]
                 if ("theta" %in% names(formals(model))) {
                     gradient(f = model,
                              x = theta,
                              desf = desf)
                 } else {
                     f012 <- do.call(what = model,
                                     args = c(list(desf = desf),
                                              as.list(theta)))
                     attr(f012, "gradient")
                 }
             })
str(f1)

bf1 <- as.matrix(Matrix::bdiag(f1))
str(bf1)
head(bf1)

# DANGER ATTENTION: cuidado com a ordem das colunas. Elas devem estar na
# ordem dos parâmetros ajustados do modelo. Será necessário trocas as
# colunas de lugar.

# Para adequar a ordem das colunas conforme a ordem das estimativas.
names(coef(nn1))
i <- matrix(seq(ncol(bf1)), nrow = length(f1), byrow = TRUE)
c(i)

# Corrige a ordem das colunas.
bf1 <- bf1[, c(i)]
head(bf1)

levelplot(cov2cor(vcov(nn1)), scales = list(x = list(rot = 90)))
# round(cov2cor(vcov(nn1)), digits = 2)

# Cholesky para eficiência.
U <- chol(vcov(nn1))

# Quantil para fazer IC.
tval <- qt(p = c(lwr = 0.025, fit = 0.5, upr = 0.975),
           df = nrow(da) - ncol(bf1))

# Erro padrão e margem de erro.
pred$se <- sqrt(apply(bf1 %*% t(U),
                      MARGIN = 1,
                      FUN = function(x) sum(x^2)))
pred$me <- outer(pred$se, tval, FUN = "*")
pred <- cbind(pred, sweep(pred$me,
                          MARGIN = 1,
                          STATS = pred$pcapu,
                          FUN = "+"))
pred$me <- NULL
str(pred)

xyplot(pcapu ~ desf | estag,
       data = da,
       col = 1,
       layout = c(5, 1),
       xlim = c(-0.2, 1.2),
       xlab = "Desfolha artifical",
       ylab = "Peso de capulhos (g)") +
    as.layer(xyplot(fit + lwr + upr ~ desf | estag,
                    data = pred,
                    type = "l",
                    col = 1,
                    lty = c(1, 2, 2)))

# Medidas de curvatura.
lapply(fit_repa, MASS::rms.curv)

#-----------------------------------------------------------------------
# Ajuste da parametrização original menos interpretável para testar
# hipótese sobre o parâmetro C.

fit_orig <- fits[seq(fits) %% 2 == 1]
names(fit_orig)
start_0 <- sapply(fit_orig, coef)
start_0

nn2 <- gnls(pcapu ~ f0 - f1 * desf^exp(C),
            data = da,
            params = f0 + f1 + C ~ -1 + estag,
            start = c(t(start_0)))
coef(nn2)

# Tabela com testes de hipótese para \theta == 0.
round(summary(nn2)$tTable, digits = 4)

#-----------------------------------------------------------------------
# Testes de hipótese.

saved_coef <- coef(nn1)

# COMMENT: não se rejeita a hipótese nula de que C == 0 para "Botão
# Floral", "Florescimento" e "Capulho".

# ATTENTION: Testar a hipótese de C == 0 é equivalente a hipótese de xde
# == 5/f1 que é o que acontece se o modelo for linear (xde é determinado
# por interpolação linear).

# Teste de hipótese usando Wald.
coef(nn1)["xde.estagCapulho"]
b <- 5/coef(nn1)["f1.estagCapulho"]
b

L <- matrix(0, nrow = 1, ncol = length(coef(nn1)))
L[1, 15] <- 1
L

# Teste para theta == 0.
anova(nn1, L = L)

# DANGER: Subtrai do valor sob hipótese dentro do objeto. Isso só é
# possível porque objetos de classe S3 permitem que você modique o
# conteúdo do objeto. Cuidado para não corromper o objeto.
nn1$coefficients[15] <- nn1$coefficients[15] - b
coef(nn1)

# Agora o teste para theta == 0 é o que se pretendia.
anova(nn1, L = L)

# Restaura os coeficientes.
nn1$coefficients <- saved_coef

# Fazendo para todos os níveis.
for (i in 11:15) {
    b <- 5/coef(nn1)[i - 5]
    cat("H_0: xde == ", b, "\n")
    cat("Estimated xde:", nn1$coefficients[i], "\n\n")
    L <- matrix(0, nrow = 1, ncol = length(coef(nn1)))
    L[1, i] <- 1
    nn1$coefficients[i] <- nn1$coefficients[i] - b
    print(anova(nn1, L = L))
    cat("\n")
}

# Restaura os coeficientes.
nn1$coefficients <- saved_coef

# ATTENTION: está sendo usada a estatítica Chi-quadrado mesmo
# explicitando a estatística F. A suposição é de que a função não esteja
# pronta para objetos de classe gls/gnls. Averiguar.

linearHypothesis(model = nn1,
                 test = "F",
                 hypothesis.matrix = paste0("xde.estagCapulho = ", b))

# Testando em todos os níveis de estag.
hip <- paste0(grep("xde", names(coef(nn1)), value = TRUE),
              " = ",
              5/coef(nn1)[grep("f1", names(coef(nn1)))])
hip

lapply(hip, linearHypothesis, model = nn1, test = "F")

#-----------------------------------------------------------------------
