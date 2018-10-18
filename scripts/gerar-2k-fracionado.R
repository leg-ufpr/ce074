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

# Gerar as frações 1/2 de 2^3 usando a interação tripla A:B:C.
da <- generate_2k_design(k = 3)
da$nms <- cell_names(da)
da

# Efeito usado para fracionar: ABC.
ef <- with(da, A * B * C)

da$fracao <- factor(ef, labels = c("principal", "complementar"))
split(da$nms, da$fracao)

#--------------------------------------------
# Usando contrastes de definição.
# Requer codificação começando em 0.

da <- generate_2k_design(k = 3, coding = "0:1")
da

# Só precisa das colunas dos efeitos principais (sem intecepto).
X <- as.matrix(da)

# A interação ABC tem coeficientes \alpha = [1, 1, 1].
alpha <- rbind(1, 1, 1)
rownames(alpha) <- names(da)
alpha

# L = (\alpha_1 * x_1 + \alpha_2 * x_2 + ··· + \alpha_k * x_k) mod 2.
l <- (X %*% alpha) %% 2
l

da$nms <- cell_names(da)
da$fracao <- factor(l, labels = c("princ", "compl"))
da

# Corridas que e seus respectivos blocos.
split(da$nms, da$fracao)

#--------------------------------------------
# Determinar os efeitos confundidos nas frações.

da <- generate_2k_design(k = 4)

# Usando a interação ABC para fracionar.
da$fracao <- with(da, A * B * C * D)
da

# Pegando uma só fração.
db <- subset(da, fracao > 0)
db

X <- model.matrix(~A * B * C * D, data = db)
colnames(X)[1] <- "I"
X

# Atenção para as entradas não nulas fora da diagonal.
xx <- t(X) %*% X
xx

# Efeitos que estão confundidos entre si.
ef_conf <- apply(xx,
                 MARGIN = 2,
                 FUN = function(column) {
                     r <- rownames(xx)[which(column != 0)]
                     r[order(nchar(r))]
                 })
ef_conf <- unique(t(ef_conf))
ef_conf

#-----------------------------------------------------------------------
# Gerar 2^k com 4 frações.

da <- generate_2k_design(k = 5)

# Usando as interações ABCD e ABE para obter frações 1/4.
b <- with(da, {
    data.frame(ef1 = factor(A * B * C * D),
               ef2 = factor(A * B * E))
})

# Combinando os 2 * 2 níveis para ter os 4 blocos.
da$fracao <- with(b,
                  factor(interaction(ef1, ef2),
                         labels = 1:4))
da

split(da, da$fracao)

#--------------------------------------------
# Para obter apenas uma das frações.

# I_1 = ABCD e I_2 = ABE, então pode-se fazer D * I_1 -> D = ABC e E *
# I_2 -> E = AB.

da <- generate_2k_design(k = 3)
da$D <- with(da, A * B * C)
da$E <- with(da, A * B)
da

#--------------------------------------------
# Usando contrastes de definição.

da <- generate_2k_design(k = 5, coding = "0:1")
X <- as.matrix(da)

# Usando ABCD e ABE.
b <- X %*% cbind(L1 = c(1, 1, 1, 1, 0),
                 L2 = c(1, 1, 0, 0, 1)) %% 2

da$frac <- apply(b, MARGIN = 1, FUN = paste0, collapse = "")
da

# Corridas e seus blocos.
da$nom <- cell_names(da)
split(da$nom, da$frac)

#=======================================================================
#-----------------------------------------------------------------------
# Análise de experimentos fatoriais fracionados.

#--------------------------------------------
# Problem 9.
# da <- read.table("clipboard", header = TRUE, sep = "\t", dec = ",")
# dput(da)
da <-
structure(list(Solvent = c(-1L, 1L, -1L, 1L, -1L, 1L, -1L, 1L,
-1L, 1L, -1L, 1L, -1L, 1L, -1L, 1L), Catalyst = c(-1L, -1L, 1L,
1L, -1L, -1L, 1L, 1L, -1L, -1L, 1L, 1L, -1L, -1L, 1L, 1L), Temperature = c(-1L,
-1L, -1L, -1L, 1L, 1L, 1L, 1L, -1L, -1L, -1L, -1L, 1L, 1L, 1L,
1L), Purity = c(-1L, -1L, -1L, -1L, -1L, -1L, -1L, -1L, 1L, 1L,
1L, 1L, 1L, 1L, 1L, 1L), Reactant.pH = c(1L, -1L, -1L, 1L, -1L,
1L, 1L, -1L, -1L, 1L, 1L, -1L, 1L, -1L, -1L, 1L), Color = c(-0.63,
2.51, -2.68, 1.66, 2.06, 1.22, -2.09, 1.93, 6.79, 5.47, 3.45,
5.68, 5.22, 4.38, 4.3, 4.05)), .Names = c("Solvent", "Catalyst",
"Temperature", "Purity", "Reactant.pH", "Color"), class = "data.frame", row.names = c(NA,
-16L))

str(da)
names(da)[1:5] <- LETTERS[1:5]

X <- model.matrix(~A * B * C * D * E, data = da)
colnames(X)[1] <- "I"
dim(X)

# Atenção para as entradas não nulas fora da diagonal.
xx <- t(X) %*% X
dim(xx)

# Efeitos que estão confundidos entre si.
ef_conf <- apply(xx,
                 MARGIN = 2,
                 FUN = function(column) {
                     r <- rownames(xx)[which(column != 0)]
                     r[order(nchar(r))]
                 })
ef_conf <- unique(t(ef_conf))
ef_conf

X[, c("A", "B:C:D:E")]

# Análise.
m0 <- lm(Color ~ A * B * C * D * E, data = da)
anova(m0)

# Estimativas.
coef(m0)

# Para fazer o qq-plot.
b <- na.omit(coef(m0)[-1])
b

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

# FrF2::DanielPlot(m0)
b[order(-abs(b))]

m0 <- lm(Color ~ A + B + C + D + E + A:D + A:B + A:C, data = da)
anova(m0)

#--------------------------------------------
# Problem 32.
# da <- read.table("clipboard", header = TRUE, sep = "\t", dec = ",")
# dput(da)
da <-
structure(list(A = c(-1L, 1L, -1L, 1L, -1L, 1L, -1L, 1L, -1L,
1L, -1L, 1L, -1L, 1L, -1L, 1L), B = c(-1L, -1L, 1L, 1L, -1L,
-1L, 1L, 1L, -1L, -1L, 1L, 1L, -1L, -1L, 1L, 1L), C = c(-1L,
-1L, -1L, -1L, 1L, 1L, 1L, 1L, -1L, -1L, -1L, -1L, 1L, 1L, 1L,
1L), D = c(-1L, -1L, -1L, -1L, -1L, -1L, -1L, -1L, 1L, 1L, 1L,
1L, 1L, 1L, 1L, 1L), E = c(1L, -1L, -1L, 1L, -1L, 1L, 1L, -1L,
-1L, 1L, 1L, -1L, 1L, -1L, -1L, 1L), y = c(7.93, 5.56, 5.77,
12, 9.17, 3.65, 6.4, 5.69, 8.82, 17.55, 8.87, 8.94, 13.06, 11.49,
6.25, 26.05)), .Names = c("A", "B", "C", "D", "E", "y"), class = "data.frame", row.names = c(NA,
-16L))

#--------------------------------------------
# Problem 37.
# da <- read.table("clipboard", header = TRUE, sep = "\t", dec = ",")
# dput(da)
da <-
structure(list(Squeegee.Pressure = c(0.1, 0.3, 0.1, 0.3, 0.1,
0.3, 0.1, 0.3, 0.1, 0.3, 0.1, 0.3, 0.1, 0.3, 0.1, 0.3, 0.1, 0.3,
0.1, 0.3, 0.1, 0.3, 0.1, 0.3, 0.1, 0.3, 0.1, 0.3, 0.1, 0.3, 0.1,
0.3), Printing.Speed = c(24L, 24L, 32L, 32L, 24L, 24L, 32L, 32L,
24L, 24L, 32L, 32L, 24L, 24L, 32L, 32L, 24L, 24L, 32L, 32L, 24L,
24L, 32L, 32L, 24L, 24L, 32L, 32L, 24L, 24L, 32L, 32L), Squeegee.Angle = c(45L,
45L, 45L, 45L, 65L, 65L, 65L, 65L, 45L, 45L, 45L, 45L, 65L, 65L,
65L, 65L, 45L, 45L, 45L, 45L, 65L, 65L, 65L, 65L, 45L, 45L, 45L,
45L, 65L, 65L, 65L, 65L), Temperature = c(20L, 20L, 20L, 20L,
20L, 20L, 20L, 20L, 28L, 28L, 28L, 28L, 28L, 28L, 28L, 28L, 20L,
20L, 20L, 20L, 20L, 20L, 20L, 20L, 28L, 28L, 28L, 28L, 28L, 28L,
28L, 28L), Viscosity = c(1125L, 1125L, 1125L, 1125L, 1125L, 1125L,
1125L, 1125L, 1125L, 1125L, 1125L, 1125L, 1125L, 1125L, 1125L,
1125L, 1275L, 1275L, 1275L, 1275L, 1275L, 1275L, 1275L, 1275L,
1275L, 1275L, 1275L, 1275L, 1275L, 1275L, 1275L, 1275L), Cleaning.Interval = c(8L,
15L, 15L, 8L, 15L, 8L, 8L, 15L, 8L, 15L, 15L, 8L, 15L, 8L, 8L,
15L, 8L, 15L, 15L, 8L, 15L, 8L, 8L, 15L, 8L, 15L, 15L, 8L, 15L,
8L, 8L, 15L), Separation.Speed = c(0.4, 0.8, 0.8, 0.4, 0.4, 0.8,
0.8, 0.4, 0.8, 0.4, 0.4, 0.8, 0.8, 0.4, 0.4, 0.8, 0.4, 0.8, 0.8,
0.4, 0.4, 0.8, 0.8, 0.4, 0.8, 0.4, 0.4, 0.8, 0.8, 0.4, 0.4, 0.8
), Relative.Humidity = c(70L, 70L, 30L, 30L, 30L, 30L, 70L, 70L,
30L, 30L, 70L, 70L, 70L, 70L, 30L, 30L, 30L, 30L, 70L, 70L, 70L,
70L, 30L, 30L, 70L, 70L, 30L, 30L, 30L, 30L, 70L, 70L), PVM = c(1,
1.04, 1.02, 0.99, 1.02, 1.01, 1.01, 1.03, 1.04, 1.14, 1.2, 1.13,
1.14, 1.07, 1.06, 1.13, 1.02, 1.1, 1.09, 0.96, 1.02, 1.07, 0.98,
0.95, 1.1, 1.12, 1.19, 1.13, 1.2, 1.07, 1.12, 1.21), NPU = c(5L,
13L, 16L, 12L, 15L, 9L, 12L, 17L, 21L, 20L, 25L, 21L, 25L, 13L,
20L, 26L, 10L, 13L, 17L, 13L, 14L, 11L, 10L, 14L, 28L, 24L, 22L,
15L, 21L, 19L, 21L, 27L)), .Names = c("Squeegee.Pressure", "Printing.Speed",
"Squeegee.Angle", "Temperature", "Viscosity", "Cleaning.Interval",
"Separation.Speed", "Relative.Humidity", "PVM", "NPU"), class = "data.frame", row.names = c(NA,
-32L))

names(da)[1:8] <- LETTERS[1:8]
str(da)

# Codificando com -1 e 1.
db <- da
db[, 1:8] <- sapply(da[, 1:8], FUN = codify_2k)
str(db)

X <- model.matrix(~(.)^8, data = db[, 1:8])
colnames(X)[1] <- "I"
dim(X)

# Atenção para as entradas não nulas fora da diagonal.
xx <- t(X) %*% X
dim(xx)

# Efeitos que estão confundidos entre si.
ef_conf <- apply(xx,
                 MARGIN = 2,
                 FUN = function(column) {
                     r <- rownames(xx)[which(column != 0)]
                     r[order(nchar(r))]
                 })
dim(ef_conf)
ef_conf <- unique(t(ef_conf))
ef_conf

# Quais são as funções geradoras?
# Qual a resolução do experimento?

m0 <- lm(PVM ~ A * B * C * D * E * F * G * H, data = db)
anova(m0)

# Para fazer o qq-plot.
b <- na.omit(coef(m0)[-1])
b[order(-abs(b))]

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

m1 <- update(m0, . ~ (A + B + C + D + E + F + G + H)^2)
anova(m1)

bb <- b[order(-abs(b))]
length(bb)

bb <- names(head(bb, n = 15))
cat(paste(sort(bb[nchar(bb) == 3]), collapse = " + "), "\n")

m2 <- update(m0, . ~ A + B + C + D + E + F + G + H + A:B + A:E + B:C + B:D + B:E + B:H + C:H + D:F + E:G)
anova(m2)

#-----------------------------------------------------------------------
