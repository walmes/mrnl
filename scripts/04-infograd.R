#-----------------------------------------------------------------------
#                                            Prof. Dr. Walmes M. Zeviani
#                                leg.ufpr.br/~walmes · github.com/walmes
#                                        walmes@ufpr.br · @walmeszeviani
#                      Laboratory of Statistics and Geoinformation (LEG)
#                Department of Statistics · Federal University of Paraná
#                                       2018-nov-18 · Curitiba/PR/Brazil
#-----------------------------------------------------------------------

#=======================================================================
# Uso de informação de gradiente.

#-----------------------------------------------------------------------
# Pacotes.

library(lattice)
library(latticeExtra)
library(rootSolve)
ls("package:rootSolve")

#-----------------------------------------------------------------------

data(wtloss, package = "MASS")

p0 <- xyplot(Weight ~ Days,
             data = wtloss,
             col = 1,
             xlab = "Dias em dieta",
             ylab = "Peso (kg)")
p0 +
    layer({
        panel.curve(f0 + f1 * 2^(-x/K), add = TRUE, col = 2)
    }, data = list(f0 = 81, f1 = 102, K = 141))

# Expressão do modelo e as derivadas parciais em relação a \theta.
model_eq <- expression(f0 + f1 * 2^(-x/K))
pars <- c("f0", "f1", "K")
sapply(pars, D, expr = model_eq)

# Cria a função que retorna o gradiente no atributo.
model_fun <- function(x, f0, f1, K, gradient = TRUE) {
    # A função f(x, \theta).
    value <- f0 + f1 * 2^(-x/K)
    # A função f'(x, \theta) como atributo.
    if (gradient) {
        attr(value, "gradient") <-
            cbind(f0 = 1,
                  f1 = 2^(-x/K),
                  K = f1 * (2^(-x/K) * (log(2) * (x/K^2))))
    }
    return(value)
}

# Avaliando a função.
model_fun(x = wtloss$Days[1:5],
          f0 = 81, f1 = 102, K = 141)

# Valores iniciais.
start <- list(f0 = 60, f1 = 200, K = 41)

# Uso da função com gradiente no ajuste.
n0 <- nls(Weight ~ model_fun(Days, f0, f1, K),
          data = wtloss,
          start = start,
          trace = TRUE)
summary(n0)

# Sem o uso do gradiente.
n1 <- nls(Weight ~ f0 + f1 * 2^(-Days/K),
          data = wtloss,
          start = start,
          trace = TRUE)
summary(n1)

#-----------------------------------------------------------------------
# Exemplo com modelos segmentados.

mlp <- deriv3(~f0 + tx * x * (x < xde) + tx * xde * (x >= xde),
              c("f0", "tx", "xde"),
              function(x, f0, tx, xde) { NULL })

# IMPORTANT: A função `deriv3()` trata modelos com funções indicadoras
# ou inequações.

# Declara a função tendo como primeiro argumento o vetor de parâmetros.
model_fun <- function(theta, x_val) {
    if (getOption("my_fun_prints", default = FALSE)) {
        print(theta)
    }
    with(as.list(theta), {
        ind <- (x_val < xde)
        f0 + tx * x_val * ind + tx * xde * (!ind)
    })
}

model_fun(c(f0 = 10, tx = 1, xde = 5),
          x_val = 0:10)

# Matriz com o vetor gradiente avaliado em cada ponto de suporte.
gradient(model_fun,
         c(f0 = 10, tx = 1, xde = 5),
         x_val = 0:10)

# ATTENTION: esse gradiente é numérico.  Sua determinação depende de
# várias avaliações da função.  Use novamente depois de passar TRUE no
# `options()`.
options(my_fun_prints = TRUE)

# Função com gradiente obtido de forma numérica.
fun_num <- function(x, f0, tx, xde) {
    value <- model_fun(c(f0 = f0, tx = tx, xde = xde),
                       x_val = x)
    attr(value, "gradient") <- gradient(model_fun,
                                        c(f0 = f0, tx = tx, xde = xde),
                                        x_val = x)
    return(value)
}
fun_num(x = 0:10, f0 = 10, tx = 1, xde = 5)

# Função com gradiente obtido de forma analítica.
fun_ana <- function(x, f0, tx, xde) {
    value <- model_fun(c(f0 = f0, tx = tx, xde = xde),
                       x_val = x)
    ind <- x < xde
    attr(value, "gradient") <- cbind(f0 = 1,
                                     tx = x * (ind) + xde * (!ind),
                                     xde = 0 * (ind) + tx * (!ind))
    return(value)
}
fun_ana(x = 0:10, f0 = 10, tx = 1, xde = 5)

#-----------------------------------------------------------------------
# Avaliação da peformance do código.

library(microbenchmark)

options(my_fun_prints = FALSE)

microbenchmark(
    "G. numérico" = fun_num(x = 0:10,
                            f0 = 10, tx = 1, xde = 5),
    "G. analítico" = fun_ana(x = 0:10,
                             f0 = 10, tx = 1, xde = 5),
    times = 1000)

# TIP: quando for usar código em simulação, preocupe-se com questões de
# performance. Muitas vezes, detalhes mínimos encurtam drasticamente o
# tempo de execução.

# ATTENTION: O exemplo mostrado aqui não significa que sempre o
# gradiente analítico será mais rápido que o numérico. Depende muito da
# complexidade do modelo, do número de dados e de parâmetros.

#-----------------------------------------------------------------------
