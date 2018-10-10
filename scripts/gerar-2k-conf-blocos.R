#-----------------------------------------------------------------------
#                                            Prof. Dr. Walmes M. Zeviani
#                                leg.ufpr.br/~walmes · github.com/walmes
#                                        walmes@ufpr.br · @walmeszeviani
#                      Laboratory of Statistics and Geoinformation (LEG)
#                Department of Statistics · Federal University of Paraná
#                                       2018-Out-10 · Curitiba/PR/Brazil
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
# Funções.

# Para gerar desenhos do fatorial 2^k completo com duas opções de
# codificação numérica dos níveis. Ver também:
# AlgDesign::gen.factorial().
generate_2k_design <- function(k = 3, coding = c("-1:1", "0:1")) {
    fct_lev <- switch(coding[1],
                      "-1:1" = c(-1L, 1L),
                      "0:1" = c(0L, 1L))
    each_fct <- replicate(k, fct_lev, simplify = FALSE)
    names(each_fct) <- LETTERS[1:k]
    design <- do.call(what = expand.grid,
                      args = c(KEEP.OUT.ATTRS = FALSE, each_fct))
    return(design)
}

# Função que retorna no nome das celas experimentais ((1), a, ab, ...) a
# partir da codificação "0:1" ou "-1:1". A função atua sobre o data
# frame com os k fatores codificados.
cell_names <- function(dataset) {
    nms <- apply(dataset,
                 MARGIN = 1,
                 FUN = function(i) {
                     u <- paste(letters[1:length(i)][i == 1],
                                collapse = "")
                     return(u)
                 })
    nms[nms == ""] <- "(1)"
    return(nms)
}

# Para calcular amplitude: max(x) - min(x).
width <- function(x) diff(range(x))

# Para codificar fatores de dois níveis.
codify_2k <- function(x, coding = c("-1:1", "0:1")) {
    coding <- coding[1]
    to_center <- function(x, coding) {
        if (coding == "-1:1") {
            x <- 2 * x - 1
        }
        return(x)
    }
    u <- unique(x)
    if (length(u) != 2) stop("O fator não tem dois níveis!")
    q <- is.character(x)
    if (q) {
        x <- as.integer(factor(x, labels = u)) - 1
        return(to_center(x, coding))
    }
    q <- is.factor(x)
    if (q) {
        x <- as.integer(x) - 1
        return(to_center(x, coding))
    }
    q <- is.numeric(x)
    if (q) {
        x <- (x - min(x))/diff(range(x))
        return(to_center(x, coding))
    }
}

# codify_2k(c(10, 20))
# codify_2k(c(10, 20), "0:1")
# codify_2k(c("A", "B"), "0:1")
# codify_2k(factor(c("A", "B")), "0:1")
# codify_2k(factor(c("A", "B")))

#-----------------------------------------------------------------------
# Gerar 2^k com confundimento em 2^1 blocos.


# Gerar o confundimento em 2^3 usando a interação tripla A:B:C.
da <- generate_2k_design(k = 3)
da

# Efeito usado para confundir com blocos: ABC.
ef <- with(da, A * B * C)
colnames(X)

da$blc <- factor(ef, labels = as.roman(1:2))
da

#--------------------------------------------
# Usando contrastes de definição.
# Requer codificação começando em 0.

da <- generate_2k_design(k = 3, coding = "0:1")
da

# A interação ABC tem coeficientes \alpha = [1, 1, 1].
alpha <- rbind(1, 1, 1)
rownames(alpha) <- names(da)
alpha

# Só precisa das colunas dos efeitos principais (sem intecepto).
X <- as.matrix(da)
X

# L = (\alpha_1 * x_1 + \alpha_2 * x_2 + ··· + \alpha_k * x_k) mod 2.
l <- (X %*% alpha) %% 2
l

da$blc <- factor(l, labels = as.roman(1:2))
da

# Corridas que e seus respectivos blocos.
da$nom <- cell_names(da)
split(da$nom, da$blc)

#--------------------------------------------
# Determinar os efeitos confundidos com os blocos.

da <- generate_2k_design(k = 3)
X <- model.matrix(~A * B * C, data = da)

# Usando a interação ABC para confundir com blocos.
da$blc <- factor(X[, "A:B:C"], labels = as.roman(1:2))
da

# Quais as colunas que não variam valores separado por bloco?
a <- aggregate(X[, -1] ~ blc,
               data = da,
               FUN = width)
a

# Identifica o termo confundido.
names(a)[sapply(a, FUN = identical, y = rep(0, nrow(a)))]

#-----------------------------------------------------------------------
# Gerar 2^k com confundimento em 2^2 blocos.

da <- generate_2k_design(k = 5)

# Usando as interações ABCD e ABE para gerar o confundimento.
b <- with(da, {
    data.frame(ef1 = factor(A * B * C * D,
                            labels = as.roman(1:2)),
               ef2 = factor(A * B * E,
                            labels = as.roman(1:2)))
})

# Combinando os 2 * 2 níveis para ter os 4 blocos.
da$blc <- with(b,
               factor(interaction(ef1, ef2),
                      labels = as.roman(1:4)))
da

# Quais as colunas que não variam valores separado por bloco?
X <- model.matrix(~A * B * C * D * E, data = da)
a <- aggregate(X[, -1] ~ blc, data = da, FUN = width)
a

# Identifica os termos confundidos.
names(a)[sapply(a, FUN = identical, y = rep(0, nrow(a)))]

#--------------------------------------------
# Usando contrastes de definição.

da <- generate_2k_design(k = 5, coding = "0:1")
X <- as.matrix(da)

# Usando ABCD e ABE.
b <- data.frame(ef1 = (X %*% rbind(1, 1, 1, 1, 0)) %% 2,
                ef2 = (X %*% rbind(1, 1, 0, 0, 1)) %% 2)
da$blc <- with(b,
               factor(interaction(ef1, ef2),
                      labels = as.roman(1:4)))
da

# Corridas e seus blocos.
da$nom <- cell_names(da)
split(da$nom, da$blc)

#-----------------------------------------------------------------------
# Análise de um 2^k com confundimento em 2^1 blocos.

# da <- read.table("clipboard", header = TRUE, sep = "\t", dec = ",")
# dput(da)
da <-
structure(list(Block = structure(c(1L, 2L, 2L, 1L, 2L, 1L, 1L,
2L, 2L, 1L, 1L, 2L, 1L, 2L, 2L, 1L, 2L, 1L, 1L, 2L, 1L, 2L, 2L,
1L, 1L, 2L, 2L, 1L, 2L, 1L, 1L, 2L), .Label = c("Block 1", "Block 2"
), class = "factor"), Aperture = structure(c(2L, 1L, 2L, 1L,
2L, 1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L, 1L,
2L, 1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L, 1L), .Label = c("large",
"small"), class = "factor"), Exposure.Time = c(-20L, -20L, 20L,
20L, -20L, -20L, 20L, 20L, -20L, -20L, 20L, 20L, -20L, -20L,
20L, 20L, -20L, -20L, 20L, 20L, -20L, -20L, 20L, 20L, -20L, -20L,
20L, 20L, -20L, -20L, 20L, 20L), Develop.Time = c(30L, 30L, 30L,
30L, 45L, 45L, 45L, 45L, 30L, 30L, 30L, 30L, 45L, 45L, 45L, 45L,
30L, 30L, 30L, 30L, 45L, 45L, 45L, 45L, 30L, 30L, 30L, 30L, 45L,
45L, 45L, 45L), Mask.Dimension = structure(c(2L, 2L, 2L, 2L,
2L, 2L, 2L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L,
2L, 2L, 2L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L), .Label = c("Large",
"Small"), class = "factor"), Etch.Time = c(14.5, 14.5, 14.5,
14.5, 14.5, 14.5, 14.5, 14.5, 14.5, 14.5, 14.5, 14.5, 14.5, 14.5,
14.5, 14.5, 15.5, 15.5, 15.5, 15.5, 15.5, 15.5, 15.5, 15.5, 15.5,
15.5, 15.5, 15.5, 15.5, 15.5, 15.5, 15.5), Yield = c(7L, 9L,
34L, 55L, 16L, 20L, 40L, 60L, 8L, 10L, 32L, 50L, 18L, 21L, 44L,
61L, 8L, 12L, 35L, 52L, 15L, 22L, 45L, 65L, 6L, 10L, 30L, 53L,
15L, 20L, 41L, 63L)), .Names = c("Block", "Aperture", "Exposure.Time",
"Develop.Time", "Mask.Dimension", "Etch.Time", "Yield"), class = "data.frame", row.names = c(NA,
-32L))

# Nomes curtos.
names(da) <- c("blc", LETTERS[1:5], "y")
str(da)

# Codifica os níveis em -1 e 1.
i <- 2:6
da[, i] <- sapply(da[, i], FUN = codify_2k)
da

# Quais as colunas que não variam valores separado por bloco?
X <- model.matrix(~A * B * C * D * E, data = da)
a <- aggregate(X[, -1] ~ blc, data = da, FUN = width)
a

# Identifica os termos confundidos.
names(a)[sapply(a, FUN = identical, y = rep(0, nrow(a)))]

# Ajuste do modelo aos dados (saturado de propósito).
m0 <- lm(y ~ blc + (A * B * C * D * E), data = da)
anova(m0)

# Note a estimativa com NA.
names(coef(m0))

# Coeficientes dos efeitos dos fatores.
b <- na.omit(coef(m0)[-c(1:2)])

# QQ-plot (inverte os eixos para anotar melhor).
qq <- qqnorm(b, plot.it = FALSE)
plot(qq$y,
     qq$x,
     # col = 2,
     col = as.integer(factor(nchar(names(b)))),
     pch = 19)
abline(v = 0, h = 0, lty = 3)
text(qq$y,
     qq$x,
     labels = names(b),
     cex = 0.6,
     pos = ifelse(qq$y < 0, 4, 2))

# Declara com interações até as duplas.
m0 <- lm(y ~ blc + (A + B + C + D + E)^2, data = da)
anova(m0)

par(mfrow = c(2, 2))
plot(m0)
layout(1)

#-----------------------------------------------------------------------
# Pacotes para geração de delineamentos com confundimento com blocos.

browseURL("https://CRAN.R-project.org/package=AlgDesign")
browseURL("https://CRAN.R-project.org/package=DoE.base")

library(AlgDesign)

dat <- gen.factorial(levels = 2, nVars = 4)
od <- optBlock(frml = ~ .,
               withinData = dat,
               blocksizes = c(8, 8),
               criterion = "OBS")
od

# Como saber o efeito usado para gerar o desenho?
od$Blocks
ncol(od$Blocks[[1]])

X <- model.matrix(~(.)^4, data = od$Blocks[[1]])
apply(X, 2, sd)

# DANGER ATTENTION: estranhamente essa função não produz a mesma
# blocagem que praticamos porque o algorítmo que ela possui gera
# experimentos com confundimentos com blocos mais gerais que os 2^k.

#-----------------------------------------------------------------------
