#-----------------------------------------------------------------------
# Pacotes e funções.

library(lattice)
library(latticeExtra)
# library(wzRfun)
source(paste0("https://raw.githubusercontent.com/walmes/",
              "wzRfun/master/R/panel.3d.contour.R"))

#-----------------------------------------------------------------------
# Análise de experimento fatorial 3^k.

# y <- scan()
# dput(y)

# Uma máquina é usada para encher recipientes de metal. A variável de
# interesse é o desperdício de líquido devido ao espumamento que
# acontece durante o enchimento. Os fatores estudados foram: tipo de
# bucal (A), velocidade de enchimento (B) e pressão de operação (C). O
# experimento foi corrido com duas repetições dos 3^3 níveis.

y <- c(-35, -25, 110, 75, 4, 5,
       -45, -60, -10, 30, -40, -30,
       -40, 15, 80, 54, 31, 36,
       17, 24, 55, 120, -23, -5,
       -65, -58, -55, -44, -64, -62,
       20, 4, 110, 44, -20, -31,
       -39, -35, 90, 113, -30, -55,
       -55, -67, -28, -26, -61, -52,
       15, -30, 110, 135, 54, 4)

da <- expand.grid(r = 1:2,
                  C = -1:1,
                  B = -1:1,
                  A = -1:1,
                  KEEP.OUT.ATTRS = FALSE)
da$y <- y
str(da)

#-----------------------------------------------------------------------
# Vizualização dos dados.

# A x B.
xyplot(y ~ A, groups = B, data = da,
       type = c("p", "a"), auto.key = TRUE)

# A x C.
xyplot(y ~ A, groups = C, data = da,
       type = c("p", "a"))

# B x C.
xyplot(y ~ B, groups = C, data = da,
       type = c("p", "a"))

# A x B x C.
xyplot(y ~ B | C, groups = A, data = da,
       type = c("p", "a"), layout = c(NA, 1))

ftable(xtabs(~A + B + C, data = da))

#-----------------------------------------------------------------------
# Ajuste do modelo.

# poly(): gera as funções polinômio. Se raw = FALSE, polinômios
# ortogonais são usados.

poly(1:4, degree = 2)
poly(1:4, degree = 2, raw = TRUE)

# Modelo saturado.
m0 <- lm(y ~ poly(A, 2) * poly(B, 2) * poly(C, 2),
         data = da)

par(mfrow = c(2, 2))
plot(m0)
layout(1)

anova(m0)

# Modelo reduzido.
m0 <- lm(y ~ (poly(A, 2) + poly(B, 2) + poly(C, 2))^2,
         data = da)

par(mfrow = c(2, 2))
plot(m0)
layout(1)

anova(m0)
summary(m0)

# NOTE: Polinômios ortogonais foram usados. Veja a matriz do modelo.
X <- model.matrix(m0)
head(X)

levelplot(t(X),
          scales = list(x = list(rot = 90)))

levelplot(crossprod(X),
          scales = list(x = list(rot = 90)))

m1 <- lm(y ~ A + I(A^2) + B + I(B^2) + C + I(C^2) + A:B + A:C + B:C,
         data = da)
anova(m1, m0)

#-----------------------------------------------------------------------
# Predição.

db <- expand.grid(A = -1:1,
                  B = seq(-1.2, 1.2, by = 0.05),
                  C = seq(-1.2, 1.2, by = 0.05))

db$y <- predict(m0, newdata = db)

# Coordenadas do delineamento em B e C.
p <- unique(da[, c("B", "C")])

levelplot(y ~ B + C | factor(A),
          data = db,
          contour = TRUE,
          layout = c(NA, 1),
          aspect = "iso") +
    layer(panel.points(x = B, y = C, pch = 19, col = 1),
          data = p) +
    layer(panel.segments(x0 = c(-1, -1, -1, -1,  0,  1),
                         y0 = c(-1,  0,  1, -1, -1, -1),
                         x1 = c( 1,  1,  1, -1,  0,  1),
                         y1 = c(-1,  0,  1,  1,  1,  1),
                         lty = 2))

#-----------------------------------------------------------------------
# Gráfico em 3D.

# display.brewer.all()
colr <- brewer.pal(11, "Spectral")
colr <- colorRampPalette(colr, space = "rgb")

# Gráfico com contornos de nível.
contourplot(y ~ B + C | factor(A),
            data = db,
            region = TRUE,
            col.regions = colr,
            as.table = TRUE,
            cuts = 20,
            layout = c(NA, 1),
            aspect = "iso")

# Gráfico com superfície e contornos.
wireframe(y ~ B + C | factor(A),
          data = db,
          scales = list(arrows = FALSE),
          col = "gray50",
          as.table = TRUE,
          col.contour = 1,
          panel.3d.wireframe = panel.3d.contour,
          type = "on",
          col.regions = colr(100),
          drape = TRUE,
          alpha.regions = 0.5,
          screen = list(z = 30, x = -60))

#-----------------------------------------------------------------------
# Decomposição das SQ em experimento fatorial 3^k.

library(lattice)
library(latticeExtra)

# O efeito do *developer strengh* (A) e *development time* (B) na
# densidade fotográfica do filme estão sendo estudados. Três níveis de
# cada fator foram utilizados e quatro repetições de um 3^2 foram
# executadas.

da <- expand.grid(rept = 1:4, A = -1:1, B = -1:1,
                  KEEP.OUT.ATTRS = FALSE)

# # da$y <- scan()
# dput(da$y)
da$y <- c(0, 5, 2, 4, 4, 7, 6, 5, 7, 8, 10, 9, 1, 4, 3, 2, 6, 7, 8, 7,
          10, 8, 10, 7, 2, 4, 5, 6, 9, 8, 10, 5, 12, 9, 10, 8)

xyplot(y ~ A,
       groups = B,
       data = da,
       type = c("p", "a"),
       jitter.x = TRUE)

p <- function(x) {
    # Retorna a matriz para os termos linear e quadrático com polinômios
    # crus.
    poly(x, degree = 2, raw = TRUE)
}

m0 <- lm(y ~ p(A) * p(B), data = da)

par(mfrow = c(2, 2))
plot(m0)
layout(1)

anova(m0)

m0 <- lm(y ~ (A + I(A^2)) * (B + I(B^2)), data = da)

# Quadro de anova com unidades de 1 GL.
anova(m0)

# Para obter o fatiamento da SQ em duas unidades de 2 GL cada. Somar 1
# neles é para passar de {-1, 0, 1} para {0, 1, 2}. Dois novos fatores,
# artificiais, serão criados a partir dos níveis e A e B. Cada um deles
# tem um contraste de definição. Cada um deles têm 3 níveis, por isso,
# vão produzir termos com 2 GL cada.

da <- within(da, {
    # Termo AB que significa A¹B¹.
    AB1 <- (1 * (A + 1) + 1 * (B + 1)) %% 3
    # Termo AB que significa A¹B².
    AB2 <- (1 * (A + 1) + 2 * (B + 1)) %% 3
    # Converte para fator.
    AB1 <- factor(AB1)
    AB2 <- factor(AB2)
})

# Eles são completamente cruzados o que implica na ortogonalidade.
xtabs(~AB1 + AB2, da)

# Não é declara interação mas sim usado os termos criados.
m1 <- lm(y ~ A + I(A^2) + B + I(B^2) + AB1 + AB2, data = da)

anova(m1)
anova(m0)

# A importância dessa partição não é para obter o quadro de anova, uma
# vez que não existe interpretação física para as partições. Elas são
# úteis para a construção de delineamentos com confundimento e com
# frações.

#-----------------------------------------------------------------------
# Estudo de caso.

da <- expand.grid(dose = 1:3 - 2,
                  brand = 1:3 - 2,
                  dog = 1:3 - 2,
                  rept = 1:2,
                  KEEP.OUT.ATTRS = FALSE)
names(da)[1:3] <- LETTERS[1:3]

# y <- scan()
# dput(y)

da$y <- c(96, 94, 101, 85, 95, 108, 84, 95, 105, 84, 99, 106, 84, 98,
          114, 83, 97, 100, 85, 98, 98, 86, 97, 109, 81, 93, 106, 84,
          95, 105, 80, 93, 110, 83, 92, 102, 85, 97, 104, 82, 99, 102,
          80, 96, 111, 86, 90, 103, 84, 95, 100, 79, 93, 108)

str(da)
ftable(xtabs(~A + B + C, data = da))

# Um pesquisador está estudando o nível de lidocaina no nível de enzima
# no coração de cachorros Beagle. Três marcas comerciais de locaina (A),
# três dosagens (B) e três cachorros (C) foram usados no experimento,
# com duas repetições de um 3^3.

m0 <- lm(y ~ p(A) * p(B) * p(C), data = da)
anova(m0)

da <- within(da, {
    AB1C1 <- factor(((1 * (A + 1) + 1 * (B + 1)) + 1 * (C + 1)) %% 3)
    AB1C2 <- factor(((1 * (A + 1) + 1 * (B + 1)) + 2 * (C + 1)) %% 3)
    AB2C1 <- factor(((1 * (A + 1) + 2 * (B + 1)) + 1 * (C + 1)) %% 3)
    AB2C2 <- factor(((1 * (A + 1) + 2 * (B + 1)) + 2 * (C + 1)) %% 3)
})

str(da)

m0 <- lm(y ~ (p(A) + p(B) + p(C))^2 + AB1C1 + AB1C2 + AB2C1 + AB2C2,
         data = da)
anova(m0)

# Confundimento com 3 blocos usando o termo AB²C².
split(da, f = da$AB2C2)

#-----------------------------------------------------------------------
# Gerando delineamentos 3^k com confundimento com blocos.

db <- expand.grid(A = 0:2, B = 0:2)
xtabs(~A + B, data = db)

db <- within(db, {
    # Termo AB que significa A¹B¹.
    AB1 <- (1 * (A + 1) + 1 * (B + 1)) %% 3
    # Termo AB que significa A¹B².
    AB2 <- (1 * (A + 1) + 2 * (B + 1)) %% 3
})

# Efeito escolhido para confundir: AB².
db$bloco <- factor(db$AB2, labels = as.roman(1:3))

split(db, f = db$bloco)

# Simulando uma resposta (ruído branco).
db$y <- rnorm(nrow(db))

anova(lm(y ~ bloco + factor(A) * factor(B), data = db))

#-----------------------------------------------------------------------
