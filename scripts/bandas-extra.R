#-----------------------------------------------------------------------
# Adicionando bandas de confiança.

data(turk0, package = "alr4")
str(turk0)

plot(Gain ~ A, data = turk0)

n0 <- nls(Gain ~ Int + (Ass - Int)* A/(Half + A),
          data = turk0,
          start = list(Int = 600, Ass = 200, Half = 0.1))
summary(n0)

plot(Gain ~ A, data = turk0, col = "gray70")
with(as.list(coef(n0)), {
    curve(Int + (Ass - Int) * A / (Half + A),
        xname = "A",
        add = TRUE,
        col = "red",
        lwd = 2)
})

tb_pred <- data.frame(A = seq(0, 0.6, length.out = 31))

#-----------------------------------------------------------------------
# {nlraa}: R package for Nonlinear Regression for Agricultural

# browseURL("https://github.com/femiguez/nlraa")
# install.packages("nlraa")
# remotes::install_github("femiguez/nlraa")
library(nlraa)
packageVersion("nlraa")

# browseURL("https://femiguez.github.io/nlraa-docs/index.html")
# browseURL("https://femiguez.github.io/nlraa-docs/Confidence-Bands.html")
# predict_nls: Monte Carlo method and model averaging.
# predict2_nls: Delta method.

aux <-
    nlraa::predict_nls(n0,
                       newdata = tb_pred[, "A", drop = FALSE],
                       interval = "confidence",
                       nsim = 2999) |>
    as.data.frame() |>
    setNames(c("fit", "stderr", "lwr", "upr"))
str(aux)

tb_pred_nlraa_boot <-
    cbind(tb_pred[, "A", drop = FALSE], aux)

plot(Gain ~ A, data = turk0)
with(tb_pred_nlraa_boot,
     matlines(x = A,
              y = cbind(fit, lwr, upr),
              lty = c(1, 2, 2),
              lwd = c(2, 1, 1),
              col = "red"))

# NOTE: Como é baseado em bootstrap, a banda pode ser assimétrica.

# Baseado em abordagem assintótica usando método Delta.
Int <- coef(n0)["Int"]
Ass <- coef(n0)["Ass"]
Half <- coef(n0)["Half"]
aux <-
    nlraa::predict2_nls(n0,
                        newdata = tb_pred[, "A", drop = FALSE],
                        interval = "confidence") |>
    as.data.frame() |>
    setNames(c("fit", "stderr", "lwr", "upr"))
str(aux)

tb_pred_nlraa_dm <-
    cbind(tb_pred[, "A", drop = FALSE], aux)

plot(Gain ~ A, data = turk0)
with(tb_pred_nlraa_dm,
     matlines(x = A,
              y = cbind(fit, lwr, upr),
              lty = c(1, 2, 2),
              lwd = c(2, 1, 1),
              col = "red"))

# DANGER: Que estranho precisar dos valores no `.GlobalEnv`.

#-----------------------------------------------------------------------
# {gslnls}: GSL Multi-Start Nonlinear Least-Squares Fitting in R

# browseURL("https://github.com/JorisChau/gslnls")
# install.packages("gslnls")
# remotes::install_github("JorisChau/gslnls")
library(gslnls)
packageVersion("gslnls")

# GSL (Levenberg-Marquardt algorithm)
# library(gslnls)  ## v1.1.1
n1 <-
    gslnls::gsl_nls(
        fn = Gain ~ Int + (Ass - Int)* A/(Half + A),
        data = transform(turk0, Gain = as.numeric(Gain)),
        start = list(Int = 600, Ass = 200, Half = 0.1)
    )
summary(n1)
class(n1)

aux <-
    predict(n1, newdata = tb_pred, interval = "confidence") |>
    as.data.frame()
tb_pred_gslnls <-
    cbind(tb_pred[, "A", drop = FALSE], aux)

predict(n1, newdata = tb_pred, interval = "prediction") |>
    head()

plot(Gain ~ A, data = turk0)
with(tb_pred_gslnls,
     matlines(x = A,
              y = cbind(fit, lwr, upr),
              lty = c(1, 2, 2),
              lwd = c(2, 1, 1),
              col = "red"))

# DANGER: Tem coisa errada nessa implementação que a banda fica mais
# estreita se a malha for mais fina.

#-----------------------------------------------------------------------
# investr: Inverse Estimation in R

# browseURL("https://github.com/bgreenwell/investr")
# install.packages("investr", dep = TRUE)
# remotes::install_github("bgreenwell/investr")
library(investr)
packageVersion("investr")

aux <-
    investr::predFit(
        object = n0,
        newdata = tb_pred[, "A", drop = FALSE],
        interval = "confidence"
    ) |>
    as.data.frame()
str(aux, max.level = 1)

tb_pred_investr <-
    cbind(tb_pred[, "A", drop = FALSE], aux)

plot(Gain ~ A, data = turk0)
with(tb_pred_investr,
     matlines(x = A,
              y = cbind(fit, lwr, upr),
              lty = c(1, 2, 2),
              lwd = c(2, 1, 1),
              col = "red"))

#-----------------------------------------------------------------------
# {propagate}: Propagation of Uncertainty

# install.packages("propagate")
library(propagate)
packageVersion("propagate")

aux <-
    propagate::predictNLS(
        model = n0,
        newdata = tb_pred[, "A", drop = FALSE],
        interval = "confidence",
        propagate = propagate::propagate(
            second.order = TRUE,
            do.sim = TRUE,
            nsim = 100)
    )
str(aux, max.level = 1)
str(aux$summary, max.level = 1)

tb_pred_propagate_prop <-
    aux$summary |>
    subset(select = c("Prop.Mean.1", "Prop.2.5%", "Prop.97.5%")) |>
    setNames(c("fit", "lwr", "upr")) |>
    cbind(tb_pred[, "A", drop = FALSE])

tb_pred_propagate_sim <-
    aux$summary |>
    subset(select = c("Sim.Mean", "Sim.2.5%", "Sim.97.5%")) |>
    setNames(c("fit", "lwr", "upr")) |>
    cbind(tb_pred[, "A", drop = FALSE])

plot(Gain ~ A, data = turk0)
with(tb_pred_propagate_prop,
     matlines(x = A,
              y = cbind(fit, lwr, upr),
              lty = c(1, 2, 2),
              lwd = c(2, 1, 1),
              col = "red"))

# ATTENTION: Muito demorado se o grid for fino.

#-----------------------------------------------------------------------
# {nlstools}.

# browseURL("https://github.com/lbbe-software/nlstools")
# install.packages("nlstools")
# remotes::install_github("lbbe-software/nlstools")
library(nlstools)
packageVersion("nlstools")

aux <-
    nlstools::nlsBootPredict(
        nlstools::nlsBoot(n0, niter = 1999),
        newdata = tb_pred[, "A", drop = FALSE],
        interval = "confidence"
    ) |>
    as.data.frame() |>
    setNames(c("fit", "lwr", "upr"))
str(aux)

tb_pred_nlstools_boot <-
    cbind(tb_pred[, "A", drop = FALSE], aux)

plot(Gain ~ A, data = turk0)
with(tb_pred_nlstools_boot,
     matlines(x = A,
              y = cbind(fit, lwr, upr),
              lty = c(1, 2, 2),
              lwd = c(2, 1, 1),
              col = "red"))

cr <- nlstools::nlsConfRegions(n0, exp = 2.5, length = 3000)
# plot(cr)

tb_pred_nlstools_confregion <-
    cr$cr |>
    # as.data.frame() |>
    # head(n = 8) |>
    apply(MARGIN = 1,
          # simplify = FALSE,
          FUN = function(x) {
              Int <- x["Int"]
              Ass <- x["Ass"]
              Half <- x["Half"]
              A <- tb_pred$A
              fit <- Int + (Ass - Int) * A/(Half + A)
              fit
          }) |>
    apply(MARGIN = 1, range) |>
    t() |>
    as.data.frame() |>
    setNames(c("lwr", "upr")) |>
    transform(A = tb_pred$A,
              fit = predict(n0, newdata = tb_pred[, "A", drop = FALSE]))

plot(Gain ~ A, data = turk0)
with(tb_pred_nlstools_confregion,
     matlines(x = A,
              y = cbind(fit, lwr, upr),
              lty = c(1, 2, 2),
              lwd = c(2, 1, 1),
              col = "red"))

#-----------------------------------------------------------------------
# Método delta do pacote {car}.

aux <-
    lapply(tb_pred$A,
           function(x) {
               car::deltaMethod(
                   n0,
                   g = glue::glue("Int + (Ass - Int) * {x}/(Half + {x})",
                                  .envir = list(x = x))
               ) |>
                   as.data.frame()
           }) |>
    do.call(what = "rbind") |>
    setNames(c("fit", "stderr", "lwr", "upr"))
rownames(aux) <- NULL

tb_pred_car <- cbind(tb_pred[, "A", drop = FALSE], aux)

#-----------------------------------------------------------------------
# Juntando todas as opções.

library(ggplot2)

rm(tb_pred_all)
tb_pred_all <-
    lapply(ls(pattern = "tb_pred_"),
           FUN = function(obj) {
               # print(str(get(obj)))
               cbind("approach" = obj,
                     get(obj)[, c("A", "fit", "lwr", "upr")])
           }) |>
    do.call(what = "rbind")
str(tb_pred_all)

ggplot(data = tb_pred_all,
       mapping = aes(x = A,
                     y = fit,
                     ymin = lwr,
                     ymax = upr)) +
    facet_wrap(facets = ~approach) +
    geom_ribbon(alpha = 0.5, fill = "orange") +
    geom_line(col = "red") +
    geom_point(data = turk0,
              mapping = aes(x = A, y = Gain),
              inherit.aes = FALSE,
              col = "gray70")

#-----------------------------------------------------------------------
