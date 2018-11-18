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

library(rpanel)
library(wzRfun)

#-----------------------------------------------------------------------
# Para compreender visualmente como um modelo não linear se comporta.

# Modelo.
model <- y ~ A/(1 + exp(-(x - M)/S))

# Define os limites para os parâmetros.
limits <- list(A = c(2, 5),
               M = c(0, 5),
               S = c(-1, 0.5, 1))

# Cria conjuntro de controles.
eyefun(model = model,
       start = limits,
       dots = list(from = 0,
                   to = 10,
                   ylim = c(0, 5),
                   lwd = 2,
                   col = "blue",
                   xlab = "Age",
                   ylab = "Weight",
                   main = expression(A/(1 + exp(-(x - M)/S)))))

#-----------------------------------------------------------------------
# Regressão linear simples.

rls <- function(panel) {
    with(panel, {
        curve(b0 + b1 * x,
              xlim[1],
              xlim[2],
              ylim = ylim,
              lwd = 2,
              ylab = ylab,
              xlab = xlab,
              main = main)
        ncur <- seq(xlim[1], xlim[2], length = nc)
        ncur <- ncur[-length(ncur)]
        lapply(ncur,
               FUN = function(x) {
                   xx <- seq(ylim[1], ylim[2], l = 100)
                   mu <- b0 + b1 * x
                   s <- sqrt(s2)
                   y <- dnorm(xx, mu, s)/scale
                   ind <- y > 1e-04
                   if (curvas) {
                       lines(y[ind] + x, xx[ind], col = "green4")
                       segments(x, mu, x + dnorm(mu, mu, s)/scale,
                                mu, col = "green4", lty = 2)
                       segments(x, mu - 4 * s, x, mu + 4 * s,
                                col = "gray70", lty = 2)
                   }
                   if (panel$simula) {
                       y <- rnorm(n, mu, s)
                       points(rep(x, n), y)
                   }
               })
    })
    return(panel)
}

panel <- rp.control(xlim = c(0, 10), ylim = c(0, 10), nc = 7, scale = 1,
                    n = 4, xlab = "x", ylab = "[Y|x]",
                    main = "Regressão linear simples")
rp.slider(panel, b0, 0, 10, title = "beta_0", initval = 5,
          showvalue = TRUE, action = rls, res = 0.05)
rp.slider(panel, b1, -4, 4, title = "beta_1", initval = 0,
          showvalue = TRUE, action = rls, res = 0.05)
rp.checkbox(panel, curvas, labels = "Densidade?", initval = FALSE,
            action = rls)
rp.slider(panel, s2, 0, 1, title = "sigma^2", initval = 0.1,
          showvalue = TRUE, action = rls, res = 0.01)
rp.checkbox(panel, simula, labels = "Simular?", initval = FALSE,
            action = rls)

#-----------------------------------------------------------------------
# Regressão não linear com erros heterocedásticos.

rnlh <- function(panel) {
    with(panel, {
        curve(b0 + A/(1 + exp(-(x - x0)/S)), xlim[1], xlim[2],
              ylim = ylim, lwd = 2, ylab = ylab, xlab = xlab,
              main = main)
        ncur <- seq(xlim[1], xlim[2], length = nc)
        ncur <- ncur[-length(ncur)]
        lapply(ncur,
               FUN = function(x) {
                   xx <- seq(ylim[1], ylim[2], l = 100)
                   mu <- b0 + A/(1 + exp(-(x - x0)/S))
                   s <- sqrt(s2) * exp(d * mu)
                   y <- dnorm(xx, mu, s)/scale
                   ind <- y > 1e-04
                   if (curvas) {
                       lines(y[ind] + x, xx[ind], col = "green4")
                       segments(x, mu, x + dnorm(mu, mu, s)/scale, mu,
                                col = "green4", lty = 2)
                       segments(x, mu - 4 * s, x, mu + 4 * s, col = "gray70",
                                lty = 2)
            }
            if (panel$simula) {
                y <- rnorm(n, mu, s)
                points(rep(x, n), y)
            }
        })
    })
    return(panel)
}

panel <- rp.control(xlim = c(0, 10), ylim = c(0, 10), nc = 6, scale = 1,
                    n = 4, xlab = "x", ylab = "[Y|x]",
                    main = "Regressão não linear - modelo logístico")
rp.slider(panel, b0, 0, 10, title = "Intercepto", initval = 0,
          showvalue = TRUE, action = rnlh, res = 0.05)
rp.slider(panel, A, 5, 10, title = "Assíntota", initval = 5,
          showvalue = TRUE, action = rnlh, res = 0.05)
rp.slider(panel, x0, 0, 10, title = "Inflexão", initval = 5,
          showvalue = TRUE, action = rnlh, res = 0.05)
rp.slider(panel, S, 0, 4, initval = 1, title = "Escala",
          showvalue = TRUE, action = rnlh, res = 0.05)
rp.checkbox(panel, curvas, labels = "Densidades?", initval = FALSE,
            action = rnlh)
rp.slider(panel, s2, 0, 0.5, title = "sigma^2", initval = 0.1,
          showvalue = TRUE, action = rnlh, res = 0.01)
rp.slider(panel, d, -0.5, 0.5, title = "delta", initval = 0,
          showvalue = TRUE, action = rnlh, res = 0.01)
rp.checkbox(panel, simula, labels = "Simular?", initval = FALSE,
            action = rnlh)

#-----------------------------------------------------------------------
# Regressão beta.

rb <- function(panel) {
    dbt <- function(x, mu, phi) {
        gamma(phi)/(gamma(phi * mu) * gamma((1 - mu) * phi)) *
            x^(mu * phi - 1) * (1 - x)^((1 - mu) * phi - 1)
    }
    with(panel, {
        curve(1/(1 + exp(-(x - x0)/S)), xlim[1], xlim[2], ylim = ylim,
              lwd = 2, ylab = ylab, xlab = xlab, main = main)
        ncur <- seq(xlim[1], xlim[2], length = nc)
        ncur <- ncur[-length(ncur)]
        lapply(ncur,
               FUN = function(x) {
                   xx <- seq(0, 1, l = 500)
                   mu <- 1/(1 + exp(-(x - x0)/S))
                   y <- dbt(xx, mu, phi)/scale
                   ind <- y > 1e-04
                   if (curvas) {
                       lines(y[ind] + x, xx[ind], col = "green4")
                       segments(x, mu, x + dbt(mu, mu, phi)/scale, mu,
                                col = "green4")
                   }
                   if (panel$simula) {
                       y <- rbeta(n, mu * phi, (1 - mu) * phi)
                       points(rep(x, n), y)
                   }
               })
    })
    return(panel)
}

panel <- rp.control(xlim = c(0, 10), ylim = c(0, 1), nc = 9, scale = 15,
                    n = 4, xlab = "x", ylab = "[Y|x]",
                    main = "Regressão beta - ligação logística")
rp.slider(panel, x0, 0, 10, title = "Inflexão (x0)", initval = 5,
          showvalue = TRUE, action = rb, res = 0.05)
rp.slider(panel, S, 0, 4, title = "Escala (S)", initval = 1,
          showvalue = TRUE, action = rb, res = 0.05)
rp.checkbox(panel, curvas, labels = "Densidade?", initval = FALSE,
            action = rb)
rp.slider(panel, phi, 0, 170, title = "Dispersão (phi)", initval = 70,
          showvalue = TRUE, action = rb, res = 1)
rp.checkbox(panel, simula, labels = "Simular?", initval = FALSE,
            action = rb)

#-----------------------------------------------------------------------
# Regressão de Poisson.

rpoisson <- function(panel) {
    with(panel, {
        curve(exp(b0 + b1 * x), xlim[1], xlim[2], ylim = ylim, lwd = 2,
              ylab = ylab, xlab = xlab, main = main)
        ncur <- seq(xlim[1], xlim[2], length = nc)
        ncur <- ncur[-length(ncur)]
        lapply(ncur,
               FUN = function(x) {
                   xx <- floor(ylim[1]):ceiling(ylim[2])
                   nu <- exp(b0 + b1 * x)
                   px <- dpois(xx, lambda = nu)
                   ind <- px > 0.001
                   if (curvas) {
                       segments(x, xx[ind], x + px[ind]/scale, xx[ind],
                                col = "green4", lwd = 2)
                   }
                   if (panel$simula) {
                       y <- rpois(n, nu)
                       points(rep(x, n), y)
                   }
               })
    })
    if (panel$fale) {
        system("espeak -v pt -s 110 'Esse é o modelo de regressão Poisson'")
        panel$fale <- FALSE
    }
    return(panel)
}

panel <- rp.control(xlim = c(0, 10), ylim = c(0, 50), nc = 6,
                    scale = 0.3, n = 4, xlab = "x", ylab = "[Y|x]",
                    main = "Regressão Poisson - ligação log")
rp.slider(panel, b0, -2, 8, title = "Intercepto", initval = 3,
          showvalue = TRUE, action = rpoisson, res = 0.1)
rp.slider(panel, b1, -1, 1, title = "Coef. angular", initval = 0,
          showvalue = TRUE, action = rpoisson, res = 0.01)
rp.checkbox(panel, curvas, labels = "Probabilidade?", initval = FALSE,
            action = rpoisson)
rp.checkbox(panel, simula, labels = "Simular?", initval = FALSE,
            action = rpoisson)
rp.checkbox(panel, fale, labels = "Comentar?", initval = FALSE,
            action = rpoisson)

#-----------------------------------------------------------------------
