#-----------------------------------------------------------------------
# Controle de Processos Industriais
#
# Análise de experimento fatorial duplo
#
#                                            Prof. Dr. Walmes M. Zeviani
#                                leg.ufpr.br/~walmes · github.com/walmes
#                                        walmes@ufpr.br · @walmeszeviani
#                      Laboratory of Statistics and Geoinformation (LEG)
#                Department of Statistics · Federal University of Paraná
#                                       2018-Set-11 · Curitiba/PR/Brazil
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
# Pacotes.

library(tidyverse)

#-----------------------------------------------------------------------
# Importação dos dados.

# Experimento com algodão em delineamento inteiramente casualizado.
# Foram avaliados o efeito de intensidade de desfolha (`desf`) e fase de
# desenvolvimento da planta na qual foi feita a desfolha (`estag`) em um
# arranjo fatorial completo 5 x 5 com 5 repetições.  Várias respostas
# foram observadas no experimento.  A unidade experimental foi um vaso
# cultivado com duas plantas de algodão.  O objetivo do experimento é
# verificar em quais fases de desenvolvimento a planta é mais afetada
# pela desfolha, qual o impacto na produção e quanto de desfolha é
# tolerado sem comprometer a produção.

alg <- read_tsv("http://leg.ufpr.br/~walmes/data/desfolha.txt")
alg

# Experimento com soja no delineamento de blocos casualizados
# completos. Foram avaliados o efeito da quantidade de adubação com
# potássio (`potassio`) e nível de irrigação (`agua`) com um arranjo
# fatorial completo 5 x 3 e 5 repetições.  Várias variáveis respostas
# foram registradas.  A unidade experimental foi um vaso cultivado com
# duas plantas de soja.  O objetivo do experimento foi verificar de
# adubação com potássio conferia a planta maior resistência ao déficit
# hídrico.

soj <- read_tsv("http://leg.ufpr.br/~walmes/data/soja.txt",
                comment = "#",
                locale = locale(decimal_mark = ","))
soj

#=======================================================================
# Análise do experimento do algodão.

#--------------------------------------------
# Exploratória.

# IMPORTANTE: Verificando a estrutura.
xtabs(~estag + desf, data = alg)

# Com função suave na média.
ggplot(data = alg,
       mapping = aes(y = pcapu, x = desf)) +
    geom_point() +
    geom_smooth(se = FALSE) +
    facet_wrap(facets = ~estag, nrow = 1) +
    xlab("Intensidade de desfolha artifical") +
    ylab("Peso de capulhos (g)")

# Com polinômio de grau 2.
ggplot(data = alg,
       mapping = aes(y = pcapu, x = desf)) +
    geom_point() +
    geom_smooth(method = lm, formula = y ~ poly(x, degree = 2)) +
    facet_wrap(facets = ~estag, nrow = 1) +
    xlab("Intensidade de desfolha artifical") +
    ylab("Peso de capulhos (g)")

#--------------------------------------------
# Expecificação do modelo.

# Cria fatores.
alg <- alg %>%
    mutate(Estag = factor(estag),
           Desf = factor(desf))

# Declara o modelo para o experimento fatorial com efeitos cruzados.
m0 <- lm(pcapu ~ Estag * Desf, data = alg)
anova(m0)

# Modelo para o experimento fatorial com efeitos aninhados.
m1 <- aov(pcapu ~ Estag/Desf, data = alg)
anova(m1)

# Mesma deviance.
deviance(m0)
deviance(m1)

# Organização dos coeficientes.
coef(m1)

# Partição das SQs para H0: efeito nulo de desfolha em cada nível de
# estágio.
summary(m1,
        split = list("Estag:Desf" = list(
                         "Desf/1veg" = c(1, 6, 11, 16),
                         "Desf/2bot" = c(1, 6, 11, 16) + 1,
                         "Desf/3flo" = c(1, 6, 11, 16) + 2,
                         "Desf/4mac" = c(1, 6, 11, 16) + 3,
                         "Desf/5cap" = c(1, 6, 11, 16) + 4)))

# Modelo para o experimento fatorial com efeitos aninhados (sentido
# contrário).
m2 <- aov(pcapu ~ Desf/Estag, data = alg)
anova(m2)

# Partição das SQs para H0: efeito nulo estágio em cada nível de
# desfolha.
summary(m2,
        split = list("Desf:Estag" = list(
                         "Estag/Desf0"   = c(1, 6, 11, 16),
                         "Estag/Desf25"  = c(1, 6, 11, 16) + 1,
                         "Estag/Desf50"  = c(1, 6, 11, 16) + 2,
                         "Estag/Desf75"  = c(1, 6, 11, 16) + 3,
                         "Estag/Desf100" = c(1, 6, 11, 16) + 4)))

#=======================================================================
# Análise do experimento da soja.

#--------------------------------------------
# Exploratória.

# IMPORTANTE: Verificando a estrutura. A obs. 74 está "corrompida".
xtabs(~potassio + agua, data = soj[-74, ])

ggplot(data = soj,
       mapping = aes(y = rengrao, x = potassio)) +
    geom_point() +
    geom_smooth(se = FALSE, span = 1) +
    facet_wrap(facets = ~agua, nrow = 1) +
    xlab("Quantidade de potássio") +
    ylab("Rendimento de grãos (g)")

#-----------------------------------------------------------------------

# Transforma para fator.
soj <- soj %>%
    mutate(A = factor(agua),
           K = factor(potassio),
           bloco = factor(bloco))

# Ajuste do modelo fatorial de efeitos cruzados em DBC.
m0 <- lm(rengrao ~ bloco + A * K, data = soj[-74, ])
anova(m0)

# Resíduos.
par(mfrow = c(2, 2))
plot(m0)
layout(1)

# Modelo fatorial de efeitos aninhados.
m0 <- aov(rengrao ~ bloco + A/K, data = soj[-74, ])
anova(m0)

# Organização dos coeficientes.
coef(m0)

# Partição das SQs em H0: efeito nulo de K em cada nível de A.
summary(m0,
        split = list("A:K" = list("K/A1" = c(1, 4, 7, 10),
                                  "K/A2" = c(1, 4, 7, 10) + 1,
                                  "K/A3" = c(1, 4, 7, 10) + 2)))

library(ExpDes)

# Estudo da interação com comparações múltiplas.
with(soj,
     fat2.rbd(factor1 = A,
              factor2 = K,
              block = bloco,
              resp = rengrao,
              fac.names = c("A", "K")))

# Acomodando o efeito de potássio com polinômio.
m0 <- aov(rengrao ~ bloco + A * potassio, data = soj[-74, ])
anova(m0)

# Coeficientes estimados.
coef(m0)

# Resíduos que apontam falta de ajuste.
par(mfrow = c(2, 2))
plot(m0)
layout(1)

# Ajuste de polinômio de 2 grau.
# m0 <- aov(rengrao ~ bloco + A * poly(potassio, raw = TRUE), data = soj[-74, ])
m1 <- aov(rengrao ~ bloco + A * (potassio + I(potassio^2)), data = soj[-74, ])
anova(m1)

# Resíduos.
par(mfrow = c(2, 2))
plot(m1)
layout(1)

# Topo da matriz do modelo.
head(model.matrix(m1))

# Para obter o valor predito de y com potássio = 10 no nível de água 50
# na média do efeito de blocos.
rbind(c(1,                    # Intercepto.
        0.2, 0.2, 0.2, 0.2,   # Efeito de blocos.
        1, 0,                 # Efeito de água.
        10, 100,              # Efeito linear e quadrático de potássio.
        10, 0,                # Interações.
        100, 0)) %*% coef(m1)

# Fazendo para um grid de valores de potássio x níveis de água.
newsoj <- expand.grid(bloco = factor("I", levels = levels(soj$bloco)),
                      potassio = seq(0, 180, 2),
                      A = levels(soj$A),
                      KEEP.OUT.ATTRS = FALSE)

# IMPORTANTE: Obter a matriz do modelo usando a mesma fórmula.
# Xpred <- model.matrix(~bloco + A * (potassio + I(potassio^2)), data = newsoj)
Xpred <- model.matrix(formula(m1)[-2], data = newsoj)
Xpred[, 2:5] <- 0.2 # 1/5 das colunas que são do efeito de bloco.
head(Xpred)

# Predição da resposta.
newsoj$predito <- c(Xpred %*% coef(m1))

# Visualização.
ggplot(data = soj,
       mapping = aes(y = rengrao, x = potassio)) +
    geom_point() +
    # geom_smooth() +
    geom_line(data = newsoj, aes(y = predito, x = potassio)) +
    facet_wrap(facets = ~A, nrow = 1) +
    xlab("Quantidade de potássio") +
    ylab("Rendimento de grãos (g)")

# Adicionar as bandas de confiança.
newsoj$se_fit <- c(sqrt(diag(Xpred %*% vcov(m1) %*% t(Xpred))))

# Quantil da distribuição t para fazer os IC.
qtl <- qt(0.975, df = df.residual(m1))

newsoj <- newsoj %>%
    mutate(lwr = predito - qtl * se_fit,
           upr = predito + qtl * se_fit)

# Visualização.
ggplot(data = soj[-74, ],
       mapping = aes(y = rengrao, x = potassio)) +
    geom_point() +
    geom_line(data = newsoj, aes(y = predito, x = potassio), size = 1) +
    geom_line(data = newsoj, aes(y = lwr, x = potassio), linetype = 2) +
    geom_line(data = newsoj, aes(y = upr, x = potassio), linetype = 2) +
    facet_wrap(facets = ~A, nrow = 1) +
    xlab("Quantidade de potássio") +
    ylab("Rendimento de grãos (g)")

#-----------------------------------------------------------------------
