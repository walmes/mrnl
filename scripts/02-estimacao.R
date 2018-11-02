# ---- der-parcial -----------------------------------------------------

# Dados artificiais.
y <- matrix(c(0.1, 0.4, 0.6, 0.9))
x <- matrix(c(0.2, 0.5, 0.7, 1.8))
# plot(y ~ x, type = "o")

# Derivada de \eta() em relação à \theta.
D(expression(1 - exp(-x^theta)), "theta")

# funções que retornam as avalições das funções
