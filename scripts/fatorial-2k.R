#-----------------------------------------------------------------------
# Controle de Processos Industriais
#
# Análise de experimento fatorial 2^k
#
#                                            Prof. Dr. Walmes M. Zeviani
#                                leg.ufpr.br/~walmes · github.com/walmes
#                                        walmes@ufpr.br · @walmeszeviani
#                      Laboratory of Statistics and Geoinformation (LEG)
#                Department of Statistics · Federal University of Paraná
#                                       2018-Set-11 · Curitiba/PR/Brazil
#-----------------------------------------------------------------------

# da <- read.table("clipboard", header = TRUE, sep = "\t")
# dput(da)
# rm(da)

# Experimento descrito no Montgomery: http://bcs.wiley.com/he-bcs/Books?action=resource&bcsId=7009&itemId=1118146921&resourceId=26715.
da <- structure(list(Cutting.Speed = c(-1L, -1L, -1L, 1L, 1L, 1L, -1L,
-1L, -1L, 1L, 1L, 1L, -1L, -1L, -1L, 1L, 1L, 1L, -1L, -1L, -1L,
1L, 1L, 1L), Tool.Geometry = c(-1L, -1L, -1L, -1L, -1L, -1L,
1L, 1L, 1L, 1L, 1L, 1L, -1L, -1L, -1L, -1L, -1L, -1L, 1L, 1L,
1L, 1L, 1L, 1L), Cutting.Angle = c(-1L, -1L, -1L, -1L, -1L, -1L,
-1L, -1L, -1L, -1L, -1L, -1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
1L, 1L, 1L, 1L), Life.Hours = c(22L, 31L, 25L, 32L, 43L, 29L,
35L, 34L, 50L, 55L, 47L, 46L, 44L, 45L, 38L, 40L, 37L, 36L, 60L,
50L, 54L, 39L, 41L, 47L)), .Names = c("Cutting.Speed", "Tool.Geometry",
"Cutting.Angle", "Life.Hours"), class = "data.frame", row.names = c(NA,
-24L))

str(da)

# Encurta nomes.
names(da) <- c("speed", "geom", "angle", "resp")

#-----------------------------------------------------------------------
# Análise exploratória.

library(tidyverse)

# Distribuição dos pontos de suporte.
ftable(xtabs(~speed + geom + angle, data = da))

ggplot(data = da,
       aes(x = speed, y = resp, color = factor(geom))) +
    geom_point() +
    stat_summary(aes(group = geom), geom = "line", fun.y = mean) +
    facet_wrap(~angle)

ggplot(data = da,
       aes(x = speed, y = resp, color = factor(angle))) +
    geom_point() +
    stat_summary(aes(group = angle), geom = "line", fun.y = mean) +
    facet_wrap(~geom)

ggplot(data = da,
       aes(x = geom, y = resp, color = factor(angle))) +
    geom_point() +
    stat_summary(aes(group = angle), geom = "line", fun.y = mean) +
    facet_wrap(~speed)

#-----------------------------------------------------------------------
# Ajuste do modelo.

# Especificação do fatorial 2^3 completo.
m0 <- lm(resp ~ speed * geom * angle, data = da)

# Resíduos.
par(mfrow = c(2, 2))
plot(m0)
layout(1)

# Quadro de anova.
anova(m0)

# Médias amostrais por nível de geom.
da %>%
    group_by(geom) %>%
    summarise(resp = mean(resp))

# Médias amostrais por níveis de speed:angle.
da %>%
    group_by(speed, angle) %>%
    summarise(resp = mean(resp)) %>%
    spread(key = angle, value = resp)

# Estimativas dos efeitos.
summary(m0)

#-----------------------------------------------------------------------
# Fazendo toda a análise forma matricial.

# Extrai a matriz do modelo.
X <- model.matrix(m0)
unique(X)

# Produtos de matrizes.
XlX <- solve(t(X) %*% X)
Xly <- t(X) %*% cbind(da$resp)

# Estimativas dos efeitos.
XlX %*% Xly

names(coef(m0))

# Médias para os níveis de geom.
cbind(1, 0, -1, 0, 0, 0, 0, 0) %*% coef(m0)
cbind(1, 0,  1, 0, 0, 0, 0, 0) %*% coef(m0)

rbind("geom: {-1}" = c(1, 0, -1, 0, 0, 0, 0, 0),
      "geom: {+1}" = c(1, 0,  1, 0, 0, 0, 0, 0)) %*% coef(m0)

# Médias para os níveis de speed:angle.
rbind("speed x angle: {-1, -1}" = c(1, -1, 0, -1, 0,  1, 0, 0),
      "speed x angle: {-1, +1}" = c(1, -1, 0,  1, 0, -1, 0, 0),
      "speed x angle: {+1, -1}" = c(1,  1, 0, -1, 0, -1, 0, 0),
      "speed x angle: {+1, +1}" = c(1,  1, 0,  1, 0,  1, 0, 0)) %*% coef(m0)

# Joga em objeto para poder reutilizar e calcular a variância.
L <- rbind(c(1, -1, 0, -1, 0,  1, 0, 0),
           c(1, -1, 0,  1, 0, -1, 0, 0),
           c(1,  1, 0, -1, 0, -1, 0, 0),
           c(1,  1, 0,  1, 0,  1, 0, 0))
L %*% coef(m0)

# summary(m0)
round(vcov(m0), 2)
V <- vcov(m0)

summary(m0)

# Erro padrão das médias ajustadas.
L %*% V %*% t(L) %>% diag() %>% sqrt()

#-----------------------------------------------------------------------
