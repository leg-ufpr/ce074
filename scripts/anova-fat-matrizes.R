#-----------------------------------------------------------------------
# CE 074 - Controle de Processos Industriais
#
# Análise de experimento fatorial de forma matricial
#
#                                            Prof. Dr. Walmes M. Zeviani
#                                leg.ufpr.br/~walmes · github.com/walmes
#                                        walmes@ufpr.br · @walmeszeviani
#                      Laboratory of Statistics and Geoinformation (LEG)
#                Department of Statistics · Federal University of Paraná
#                                       2018-Ago-28 · Curitiba/PR/Brazil
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
# Obtenção do quadro de anova com operações matriciais.

# Experimento fatorial completo 3 x 2 em DIC com 9 repetições.
xtabs(~tension + wool, warpbreaks)

# Pegando uma fração dos dados para ter matrizes de dimensão menor.
da <- subset(warpbreaks, wool == "A", select = c(1, 3))
da <- da[seq(1, to = nrow(da), by = 3), ]
names(da) <- c("y", "trt")
da

# Para exemplificar será considerado é um experimento de um fator
# apenas.

# Para criar matrizes de projeção.
proj <- function(X) {
    X %*% solve(t(X) %*% X) %*% t(X)
}

# Para obter o traço de matrizes.
tr <- function(X) {
    sum(diag(X))
}

# Soma de quadrados total.
with(da, sum(y^2))

# Soma de quadrados total corrigida pela média.
with(da, sum((y - mean(y))^2))

# Projeção no espaço das observaçoes (n).
I <- diag(nrow(da))
y <- cbind(da$y)

# Modelo: y_{ijk} = \mu + \epsilon_{ijk}.
X_mu <- model.matrix(~1, data = da)
H_mu <- proj(X_mu)

# Modelo: y_{ijk} = \tau_{i} + \epsilon_{ijk}.
X_tau <- model.matrix(~0 + trt, data = da)
H_tau <- proj(X_tau)

# Modelo: y_{ijk} = \mu + \alpha_{i} + \epsilon_{ijk} (\tau_{i} = \mu + \alpha_{i}).
X_mualpha <- model.matrix(~trt, data = da)
H_mualpha <- proj(X_mualpha)

# Embora o modelo esteja paramatrizado diferente, ainda é o mesmo.
H_mualpha == H_tau

# Média geral do experimento.
H_mu %*% y
mean(y)

# Somas de quadrados para o efeito do fator A.
t(y) %*% H_tau %*% y - t(y) %*% H_mu %*% y
t(y) %*% (H_tau -  H_mu) %*% y

# SQ residual.
t(y) %*% (I - H_tau) %*% y

# Graus de liberdade.
tr(H_tau -  H_mu)
tr(I - H_tau)

# Anova para verificação.
anova(lm(y ~ 1, data = da))
anova(lm(y ~ trt, data = da))

#-----------------------------------------------------------------------
# Análise do experimento em blocos.

xtabs(~tension + wool, warpbreaks)

da <- warpbreaks
da <- da[seq(1, to = nrow(da), by = 3), ]
names(da) <- c("y", "blc", "trt")
xtabs(~blc + trt, da)

da

# Quadro de análise de variância considerando que foi um ensaio em
# blocos com um único fator.
anova(lm(y ~ blc + trt, data = da))

# Matriz coluna com a resposta.
y <- cbind(da$y)

# Matrizes do modelo.
X_1 <- model.matrix(y ~ 1, data = da)
X_12 <- model.matrix(y ~ 1 + blc, data = da)
X_123 <- model.matrix(y ~ 1 + blc + trt, data = da)

# As matrizes de projeção.
H_1 <- proj(X_1)
H_12 <- proj(X_12)
H_123 <- proj(X_123)
I <- diag(nrow(da))

# Somas de quadrados.
t(y) %*% (H_12 - H_1) %*% y
t(y) %*% (H_123 - H_12) %*% y
t(y) %*% (I - H_123) %*% y

# Graus de liberdade.
tr(H_12 - H_1)
tr(H_123 - H_12)
tr(I - H_123)

#-----------------------------------------------------------------------

X_mu <- model.matrix(y ~ 1, data = da)
X_trt <- model.matrix(y ~ 1 + trt, data = da)
X_blc <- model.matrix(y ~ 1 + blc, data = da)
X_blctrt <- model.matrix(y ~ 1 + blc + trt, data = da)

H_mu <- proj(X_mu)
H_trt <- proj(X_trt)
H_blc <- proj(X_blc)
H_blctrt <- proj(X_blctrt)

# As projeções são ortogonais!
(H_trt - H_mu) %*% (H_blc - H_mu)

# As projeções são ortogonais!
(H_trt - H_mu) %*% (H_blctrt - H_trt)

# IMPORTANT: no caso acima, as somas de quadrados sequenciais e
# marginais serão iguais porque as projeções são ortogonais.

#-----------------------------------------------------------------------
# Análise do experimento fatorial.

da <- warpbreaks
names(da) <- c("y", "A", "B")

library(ggplot2)

gg1 <-
    ggplot(da, aes(y = y, x = A, color = B)) +
    geom_jitter(width = 0.025) +
    stat_summary(aes(group = B), geom = "line", fun.y = "mean") +
    scale_color_discrete(name = "Fator B") +
    xlab("Fator A") +
    ylab("Resposta") +
    theme(legend.position = c(0.95, 0.95),
          legend.justification = c(1, 1))

gg2 <-
    ggplot(da, aes(y = y, x = B, color = A)) +
    geom_jitter(width = 0.025) +
    stat_summary(aes(group = A), geom = "line", fun.y = "mean") +
    scale_color_discrete(name = "Fator A") +
    xlab("Fator B") +
    ylab("Resposta") +
    theme(legend.position = c(0.95, 0.95),
          legend.justification = c(1, 1))

gridExtra::grid.arrange(gg1, gg2, ncol = 2)

# São todos modelos equivalentes.
# y_{abr} = \mu + \alpha_{a} + \beta_{b} + \gamma_{ab} + \epsilon_{abr}.
m0 <- aov(y ~ A * B, data = da)
m0 <- aov(y ~ B * A, data = da)
m0 <- aov(y ~ A + B + A:B, data = da)

anova(m0)      # SQ sequencial.
car::Anova(m0) # SQ marginal.

# y_{abr} = \mu + \alpha_{a} + \theta_{a(b)} + \epsilon_{abr}.
mAB <- aov(y ~ A/B, data = da)
anova(mAB)
coef(mAB)

# Partição da soma de quadrados em hipóteses mais específicas.
aAB <- summary(mAB, split = list("A:B" = list(
                                     "B|A1" = c(1, 3), # H_0: ef. de B em A1.
                                     "B|A2" = c(2, 4)  # H_0: ef. de B em A2.
                                 )))
aAB

# y_{abr} = \mu + \beta_{b} + \psi_{b(a)} + \epsilon_{abr}.
mBA <- aov(y ~ B/A, data = da)
anova(mBA)
coef(mBA)

aBA <- summary(mBA, split = list("B:A" = list(
                                     "A|B1" = 1, # H_0: ef. de A em B1.
                                     "A|B2" = 2, # H_0: ef. de A em B2.
                                     "A|B3" = 3  # H_0: ef. de A em B3.
                                 )))
aBA

#--------------------------------------------
# No sentido A com B dentro de A.

# Como essas somas de quadrados são obtidas?
y <- cbind(da$y)
X <- model.matrix(mAB)
unique(X)

X_mu <- X[, 1]
X_A <- X[, 1:2]
X_B_A1 <- X[, c(1:2, 3, 5)]
X_B_A2 <- X[, c(1:2, 4, 6)]

H_mu <- proj(X_mu)
H_A <- proj(X_A)
H_B_A1 <- proj(X_B_A1)
H_B_A2 <- proj(X_B_A2)
I <- diag(nrow(y))

L <- list(A = H_A - H_mu,
          B_A1 = H_B_A1 - H_A,
          B_A2 = H_B_A2 - H_A)

t(sapply(L,
         FUN = function(mat_proj) {
             c(DF = tr(mat_proj),
               SQ = t(y) %*% mat_proj %*% y)
         }))
aAB

#--------------------------------------------
# No sentido B com A dentro de B.

# Como essas somas de quadrados são obtidas?
y <- cbind(da$y)
X <- model.matrix(mBA)
unique(X)

X_mu <- X[, 1]
X_B <- X[, 1:3]
X_A_B1 <- X[, c(1:3, 4)]
X_A_B2 <- X[, c(1:3, 5)]
X_A_B3 <- X[, c(1:3, 6)]

H_mu <- proj(X_mu)
H_B <- proj(X_B)
H_A_B1 <- proj(X_A_B1)
H_A_B2 <- proj(X_A_B2)
H_A_B3 <- proj(X_A_B3)
I <- diag(nrow(y))

L <- list(B = H_B - H_mu,
          A_B1 = H_A_B1 - H_B,
          A_B2 = H_A_B2 - H_B,
          A_B3 = H_A_B3 - H_B)

t(sapply(L,
         FUN = function(mat_proj) {
             c(DF = tr(mat_proj),
               SQ = t(y) %*% mat_proj %*% y)
         }))
aBA

#-----------------------------------------------------------------------
# Essas mesmas hipóteses (de forma marginal) podem ser avaliadas pelo
# teste de Wald.

beta <- coef(mAB)
V <- vcov(mAB)

# Testar se o efeito de B dentro de A1 é nulo.
K <- rbind(c(0, 0, 1, 0, 0, 0),
           c(0, 0, 0, 0, 1, 0))

# Estimativa e erro padrão.
pest <- K %*% beta       # Estimativa pontual.
vest <- K %*% V %*% t(K) # Variância.

# Estatística F.
(t(pest) %*% solve(vest) %*% pest)/nrow(K)

# Testar se o efeito de B dentro de A2 é nulo.
K <- rbind(c(0, 0, 0, 1, 0, 0),
           c(0, 0, 0, 0, 0, 1))

# Estimativa e erro padrão.
pest <- K %*% beta       # Estimativa pontual.
vest <- K %*% V %*% t(K) # Variância.

# Estatística F.
(t(pest) %*% solve(vest) %*% pest)/nrow(K)

#-----------------------------------------------------------------------
