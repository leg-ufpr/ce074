#-----------------------------------------------------------------------
# Gerar 2^k com confundimento em 2^1 blocos.

# Para gerar desenhos do fatorial 2^k completo com duas opções de
# codificação numérica dos níveis.
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

# Gerar o confundimento em 2^3 usando a interação tripla A:B:C.
da <- generate_2k_design(k = 3)

X <- model.matrix(~A * B * C, data = da)
colnames(X)

da$blc <- factor(X[, "A:B:C"], labels = as.roman(1:2))
da

#--------------------------------------------
# Usando contrastes de definição.
# Requer codificação começando em 0.

da <- generate_2k_design(k = 3, coding = "0:1")
da

# A interação A:B:C tem coeficientes \alpha = [1, 1, 1].
alpha <- rbind(1, 1, 1)
rownames(alpha) <- names(da)
alpha

# Só precisa das colunas dos efeitos principais (sem intecpetop).
X <- as.matrix(da)

# L = (\alpha_1 * x_1 + \alpha_2 * x_2 + ··· + \alpha_k * x_k) mod 2.
l <- (X %*% alpha) %% 2
l

da$blc <- factor(l, labels = as.roman(1:2))
da

da$nom <- apply(da,
                MARGIN = 1,
                FUN = function(i) {
                    u <- paste(letters[1:length(i)][i == 1],
                               collapse = "")
                    if (u == "") {
                        u <- "(1)"
                    }
                    return(u)
                })

#--------------------------------------------
# Determinar os efeitos confundidos com os blocos.

da <- generate_2k_design(k = 3)
X <- model.matrix(~A * B * C, data = da)
da$blc <- factor(X[, "A:B:C"], labels = as.roman(1:2))
da

# Quais as colunas que não variam valores separado por bloco?
a <- aggregate(X[, -1] ~ blc, data = da, FUN = sd)
a

names(a)[sapply(a, FUN = identical, y = rep(0, nrow(a)))]

#-----------------------------------------------------------------------
# Gerar 2^k com confundimento em 2^2 blocos.


#-----------------------------------------------------------------------
# Análise de um 2^k com confundimento em 2^1 blocos.


#-----------------------------------------------------------------------
# Análise de um 2^k com confundimento em 2^2 blocos.

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
apply(X, 2, sum)

# ATTENTION: estranhamente essa função não produz a mesma blocagem que
# praticamos porque o algorítmo que ela possui gera experimentos com
# confundimentos com blocos mais gerais que os 2^k.

#-----------------------------------------------------------------------
