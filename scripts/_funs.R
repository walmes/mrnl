#-----------------------------------------------------------------------
#                                            Prof. Dr. Walmes M. Zeviani
#                                leg.ufpr.br/~walmes · github.com/walmes
#                                        walmes@ufpr.br · @walmeszeviani
#                      Laboratory of Statistics and Geoinformation (LEG)
#                Department of Statistics · Federal University of Paraná
#                                       2018-nov-18 · Curitiba/PR/Brazil
#-----------------------------------------------------------------------

build_deriv3 <- function(model) {
    #
    # @param model Um objeto retornado pela função \code{nls()}.
    #     Assume-se que a chamada da função \code{nls()} foi feita
    #     usando-se o parâmero \code{data}.
    #
    # @return O retorno da chamada da \code{deriv3()}.
    #
    model_vars <- all.vars(formula(model))
    model_pars <- names(coef(model))
    model_input <- setdiff(model_vars, model_pars)[-1]
    model_input <- setNames(object = replicate(n = length(model_input),
                                               expr = numeric(1),
                                               simplify = FALSE),
                            nm = model_input)
    model_fun  <- function() NULL
    formals(model_fun) <- c(model_input,
                            as.list(coef(model)))
    model_der <- deriv3(expr = formula(model)[-2],
                        namevec = model_pars,
                        function.arg = model_fun)
    return(model_der)
}
cat("Função `build_deriv3()` carregada.\n")

# Retorna as bandas de confiança para um modelo não linear.
confbands <- function(object, newdata, conf = 0.95) {
    #
    # @param object Um objeto retornado pela função \code{nls()}.
    #
    # @param newdata Um \code{data.frame} contendo as variáveis
    #     preditoras do modelo.
    #
    # @param conf Nível de confiança do intervalo.
    #
    # @return Uma matriz de 3 colunas com o limite inferior e superior
    #     do intervalo de confiança e também os valores preditos.
    #
    der <- build_deriv3(object)
    u <- do.call(der,
                 args = c(as.list(newdata),
                          as.list(coef(object))))
    X <- attr(u, "gradient")
    U <- chol(vcov(object))
    se <- sqrt(apply(X %*% t(U),
                     MARGIN = 1,
                     FUN = function(x) sum(x^2)))
    prob_tail <- (1 - conf)/2
    tval <- qt(p = c(lwr = prob_tail, fit = 0.5, upr = 1 - prob_tail),
               df = df.residual(object))
    me <- outer(X = se, Y = tval, FUN = "*")
    sweep(x = me, MARGIN = 1, STATS = u, FUN = "+")
}
cat("Função `confbands()` carregada.\n")
