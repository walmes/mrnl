#-----------------------------------------------------------------------
#                                            Prof. Dr. Walmes M. Zeviani
#                                leg.ufpr.br/~walmes · github.com/walmes
#                                        walmes@ufpr.br · @walmeszeviani
#                      Laboratory of Statistics and Geoinformation (LEG)
#                Department of Statistics · Federal University of Paraná
#                                       2018-nov-17 · Curitiba/PR/Brazil
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
# Instalação dos pacotes.

# Principais pacotes que serão usados no Curso.
pkg <- c("bbmle",
         "car",
         "formatR",
         "lattice",
         "latticeExtra",
         "MASS",
         "nlme",
         "nls2",
         "nlstools",
         "plyr",
         "proto",
         "rootSolve",
         "rpanel",
         "segmented")

# Instalação dos pacotes.
install.packages(pkg,
                 dependencies = TRUE,
                 repos = "http://cran-r.c3sl.ufpr.br/")

# Carregamento dos pacotes em série.
lapply(pkg, FUN = library, character.only = TRUE)

# Definições da sesssão.
sessionInfo()

#-----------------------------------------------------------------------
# Funções que serão usadas.

# No passado, tais funções estavam disponíveis nos endereços.
# browseURL("http://fisher.osu.edu/~schroeder_9/AMIS900/blockdiag.R")
# browseURL("http://nls2.googlecode.com/svn-history/r4/trunk/R/as.lm.R")

# O função `as.lm()` foi incorporada ao pacote wzRfun. Instale o pacote
# ou carrege usando `source()`.
source("https://raw.githubusercontent.com/walmes/wzRfun/master/R/as.lm.R")

# Opções para soma bloco diagonal de matrizes.
#  limma::blockDiag().
#  magic::adiag().
#  Matrix::bdiag().

# NOTE: retorna uma matriz esparsa.
Matrix::bdiag(matrix(1, 3, 3), matrix(10, 2, 2))

#-----------------------------------------------------------------------
