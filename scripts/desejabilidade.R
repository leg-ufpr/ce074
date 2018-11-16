#-----------------------------------------------------------------------
#                                            Prof. Dr. Walmes M. Zeviani
#                                leg.ufpr.br/~walmes · github.com/walmes
#                                        walmes@ufpr.br · @walmeszeviani
#                      Laboratory of Statistics and Geoinformation (LEG)
#                Department of Statistics · Federal University of Paraná
#                                       2018-Nov-13 · Curitiba/PR/Brazil
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
# Pacotes e funções.

library(latticeExtra)
library(gridExtra)
library(rsm)
ls("package:rsm")

source(paste0("https://raw.githubusercontent.com/walmes",
              "/wzRfun/master/R/panel.3d.contour.R"))

trellis.par.set(regions = list(col = rev(heat.colors(100))))

#-----------------------------------------------------------------------
# Tabela 11-6, pág. 422, Montgomety 5th ed.

r <- sqrt(2)
da <- rbind(expand.grid(x2 = c(-1, 1), x1 = c(-1, 1)),
            data.frame(x1 = 0, x2 = 0)[rep(1, 5), ],
            data.frame(x1 = c(r, -r, 0, 0),
                       x2 = c(0, 0, r, -r)))
db <- rbind(expand.grid(B = c(170, 180), A = c(80, 90)),
            data.frame(A = 85, B = 175)[rep(1, 5), ],
            data.frame(A = c(92.07, 77.93, 85, 85),
                       B = c(175, 175, 182.07, 167.93)))
da <- cbind(da, db)
rm("db")
# x <- scan()
# dput(x)
da$y1 <- c(765, 770, 780, 795, 799, 803, 800, 797, 798, 784, 756, 785,
           770)/10
da$y2 <- c(62, 60, 66, 59, 72, 69, 68, 70, 71, 68, 71, 58, 57)
da$y3 <- c(2940, 3470, 3680, 3890, 3480, 3200, 3410, 3290, 3500, 3360,
           3020, 3630, 3150)

str(da)

# A: tempo.
# B: temperatura.
# y1: rendimento.
# y2: viscosidade.
# y3: peso molecular.

#-----------------------------------------------------------------------
# Visualizando o delineamento.

# Visão do delineamento DCC.
xyplot(A ~ B, data = da)

# Adição dos vértices e eixos.
xyplot(x1 ~ x2, data = da, aspect = "iso", pch = 19) +
    layer(panel.segments(c(-r, 0),
                         c(0, -r),
                         c(r, 0),
                         c(0, r),
                         col = 1, lty = 2)) +
    layer(panel.lines(c(-1, 1, 1, -1, -1),
                      c(-1, -1, 1, 1, -1),
                      col = 1, lty = 2))

#-----------------------------------------------------------------------
# A direção de maior inclinação.

library(rpanel)

pred <- expand.grid(x1 = seq(-r, r, length.out = 20),
                    x2 = seq(-r, r, length.out = 20))
X <- model.matrix(~x1 * x2, data = pred)

action <- function(input) {
    pred$y <- X %*% c(1, input$beta_1, input$beta_2, input$beta_12)
    print(
        levelplot(y ~ x1 + x2,
                  data = pred,
                  contour = TRUE,
                  aspect = "iso") +
        layer({
            panel.abline(v = 0, h = 0, lty = 2)
            panel.abline(a = 0, b = delta, col = 2)
            # input$beta_2/input$beta_1
        }, data = list(delta = input$beta_2/input$beta_1)))
    return(input)
}

input <- rp.control()
rp.slider(panel = input,
          variable = "beta_1",
          from = -1,
          to = 1,
          initval = 0.1,
          showvalue = TRUE,
          action = action)
rp.slider(panel = input,
          variable = "beta_2",
          from = -1,
          to = 1,
          initval = 0.1,
          showvalue = TRUE,
          action = action)
rp.slider(panel = input,
          variable = "beta_12",
          from = -1,
          to = 1,
          initval = 0,
          showvalue = TRUE,
          action = action)
rp.do(panel = input, action = action)

#-----------------------------------------------------------------------
# Ajuste para a primeira resposta `y1`.

m0 <- lm(y1 ~ x1 + x2 + x1:x2 + I(x1^2) + I(x2^2),
         data = da)

coef(m0)
summary(m0)

# Ajuste com o pacote rsm.
m0_rsm <- rsm(y1 ~ SO(x1, x2), data = da)
summary(m0_rsm)

# Ponto estacionário.
xs(m0_rsm)

# Malha para a predição.
pred <- expand.grid(x1 = seq(-r, r, length.out = 20),
                    x2 = seq(-r, r, length.out = 20))
pred$y1 <- predict(m0, newdata = pred)

wireframe(y1 ~ x1 + x2, data = pred,
          panel.3d.wireframe = panel.3d.contour,
          type = "on",
          drape = TRUE)

levelplot(y1 ~ x1 + x2, data = pred, contour = TRUE) +
    layer({
        xs <- xs(m0_rsm)
        panel.abline(v = xs[1], h = xs[2], lty = 2)
    })

#-----------------------------------------------------------------------
# Cálculo e análise do ponto estacionário de forma detalhada.

# Coeficientes de primeira ordem.
b <- cbind(coef(m0)[2:3])

# Coeficientes de segunda ordem.
B <- matrix(0, 2, 2)
diag(B) <- coef(m0)[4:5]
B[2, 1] <- coef(m0)[6]/2
B[1, 2] <- coef(m0)[6]/2

# Coordenada do ponto estacionário.
# xs <- -0.5 * solve(B) %*% b
xs <- -0.5 * solve(B, b)

levelplot(y1 ~ x1 + x2, data = pred, contour = TRUE) +
    layer(panel.abline(v = xs[1], h = xs[2], lty = 2)) +
    layer(panel.points(xs[1], xs[2], pch = 19, col = "red"))

# Estudo da concavidade pelo sinal dos autovalores.
eigen(B)

# Ponto estácionário na escala natural.
(xs[1] * 5) + 85
(xs[2] * 5) + 185

# Quando usar a função `rms()`.
xs(m0_rsm)
m0_rsm$b
m0_rsm$B

#-----------------------------------------------------------------------
# Análise das demais respostas.

# m1 <- m0
# m2 <- lm(y2 ~ x1 + x2 + x1:x2 + I(x1^2) + I(x2^2),
#          data = da)
# m3 <- lm(y3 ~ x1 + x2 + x1:x2 + I(x1^2) + I(x2^2),
#          data = da)
# summary(m1)
# summary(m2)
# summary(m3)

m1_rsm <- m0_rsm
m2_rsm <- rsm(y2 ~ SO(x1, x2), data = da)
m3_rsm <- rsm(y3 ~ SO(x1, x2), data = da)

summary(m2_rsm)
summary(m3_rsm)

m3_rsm <- rsm(y3 ~ FO(x1, x2), data = da)

# Predição mais fina para ter contornos mais arredondados.
pred <- expand.grid(x1 = seq(-r, r, length.out = 80),
                    x2 = seq(-r, r, length.out = 80))

# Predição das 3 respostas.
pred$y1 <- predict(m1_rsm, newdata = pred)
pred$y2 <- predict(m2_rsm, newdata = pred)
pred$y3 <- predict(m3_rsm, newdata = pred)

# Argumentos comuns a todas as chamadas da `levelplot()`.
rest <- list(data = pred,
             contour = TRUE,
             labels = TRUE,
             aspect = "iso",
             drape = TRUE)

lv1 <- do.call(levelplot,
               args = c(x = y1 ~ x1 + x2, main = "Produção (y1)", rest))
lv2 <- do.call(levelplot,
               args = c(x = y2 ~ x1 + x2, main = "Viscosidade (y2)", rest))
lv3 <- do.call(levelplot,
               args = c(x = y3 ~ x1 + x2, main = "Peso molecular (y3)", rest))

grid.arrange(lv1, lv2, lv3, nrow = 1)

# Sobreposição dos gráficos de contorno.
contourplot(y1 ~ x1 + x2,
            data = pred,
            col = 1,
            key = list(lines = list(col = c(1, 2, 4)),
                       text = list(c("Produção",
                                     "Viscosidade",
                                     "Peso molecular"))),
            aspect = "iso") +
    as.layer(contourplot(y2 ~ x1 + x2,
                         data = pred,
                         col = 2,
                         aspect = "iso")) +
    as.layer(contourplot(y3 ~ x1 + x2,
                         data = pred,
                         col = 4,
                         aspect = "iso"))

# Gráfico com as restrições.
#   y1 tal que y1 (sem retrições.)
#   y2 tal que 62 <= y2 <= 68.
#   y3 tal que y3 <= 3400.

lv1 <- levelplot(y1 ~ x1 + x2,
                 data = pred,
                 contour = TRUE,
                 labels = TRUE,
                 col.regions = grey.colors,
                 key = list(lines = list(col = c(1, 2, 4)),
                            text = list(c("Produção",
                                          "Viscosidade",
                                          "Peso molecular"))),
                 aspect = "iso")
lv2 <- levelplot((y2 <= 68 & y2 >= 62) ~ x1 + x2,
                 data = pred,
                 at = c(0.9, 1.1),
                 col.regions = rgb(1, 0, 0, 0.25),
                 contour = TRUE,
                 col = 2,
                 aspect = "iso")
lv3 <- levelplot((y3 <= 3400) ~ x1 + x2,
                 data = pred,
                 at = c(0.9, 1.1),
                 col.regions = rgb(0, 0, 1, 0.25),
                 contour = TRUE,
                 col = 4,
                 aspect = "iso")
grid.arrange(lv1, lv2, lv3, nrow = 1)

lv1 + as.layer(lv2) + as.layer(lv3)

optfunc <- function(x1,  x2) {
    # Predito de cada resposta para x1 e x2 fornecidos.
    y1 <- predict(m1_rsm,  newdata = list(x1 = x1, x2 = x2))
    y2 <- predict(m2_rsm,  newdata = list(x1 = x1, x2 = x2))
    y3 <- predict(m3_rsm,  newdata = list(x1 = x1, x2 = x2))
    y1 * (y2 <= 68 & y2 >= 62) * (y3 <= 3400)
}

# Qual o ponto de máximo (na **malha** criada).
pred$z <- optfunc(pred$x1, pred$x2)
i <- which.max(pred$z)
pred[i, ]

# Gráfico de y1 para valores admissíveis de y2 e y3.
levelplot(z ~ x1 + x2,
          data = subset(pred, z > 0),
          contour = TRUE,
          labels = TRUE,
          # col.regions = grey.colors,
          aspect = "iso") +
    layer(panel.points(pred[i, "x1"],
                       pred[i, "x2"],
                       col = "green",
                       cex = 2,
                       pch = 19))

#-----------------------------------------------------------------------
# Funções de desejabilidade (desirability functions).

# Cresce para a direita = desejável à direita.
fr <- function(y, lwr, tgt, r) {
    (0 + ((y - lwr)/(tgt - lwr)) * (y >= lwr & y <= tgt) + (y > tgt))^r
}

# Cresce para a esqueda = desejável à esquerda.
fl <- function(y, upr, tgt, r) {
    (0 + ((upr - y)/(upr - tgt)) * (y >= tgt & y <= upr) + (y < tgt))^r
}

# Desejável no meio.
fb <- function(y, lwr, upr, tgt, r1, r2) {
    (((y - lwr)/(tgt - lwr)) * (y >= lwr & y <= tgt))^r1 +
    (((upr - y)/(upr - tgt)) * (y > tgt & y <= upr))^r2
}

x <- seq(0.5, 2.5, length.out = 100)
par(mfrow = c(1, 3))
curve(fr(x, lwr = 1, tgt = 2, r = 0.7), 0.5, 2.5)
curve(fl(x, upr = 2, tgt = 1, r = 0.7), 0.5, 2.5)
curve(fb(x, lwr = 1, upr = 2, tgt = 1.5, r1 = 1.1, r2 = 0.9),
      0.5, 2.5)
layout(1)

#-----------------------------------------------------------------------

# Desejabilidade (d) de cada resposta.
pred$d1 <- fr(pred$y1, lwr = 70, tgt = 80, r = 1)
pred$d2 <- fb(pred$y2, lwr = 62, upr = 68, tgt = 65,
              r1 = 0.9, r2 = 0.9)
pred$d3 <- with(pred, y3 >= 3200 & y3 <= 3600)

# Argumentos comuns.
rest <- list(data = pred,
             contour = TRUE,
             labels = TRUE,
             aspect = "iso",
             drape = TRUE)

lv1 <- do.call(levelplot,
               args = c(x = d1 ~ x1 + x2, main = "Desejabilidade da produção (d1)", rest))
lv2 <- do.call(levelplot,
               args = c(x = d2 ~ x1 + x2, main = "Desejabilidade da viscosidade (d2)", rest))
lv3 <- do.call(levelplot,
               args = c(x = d3 ~ x1 + x2, main = "Desejabilidade do peso molecular (d3)", rest))

grid.arrange(lv1, lv2, lv3, nrow = 1)

# Desejabilidade média (é uma média geométrica!).
pred$desir <- with(pred, ((d1 * d2 * d3)^(1/3)))

# Superfície da desejabilidade.
lv <- levelplot(desir ~ x1 + x2,
                data = pred,
                contour = TRUE,
                labels = TRUE,
                aspect = "iso",
                drape = TRUE)
lv

# Qual o ponto de máximo da desejabilidade (na região discreta criada).
i <- which.max(pred$desir)
pred[i, ]

# Gráfico de y1 para valores admissíveis de y2 e y3.
lv +
    layer({
    panel.abline(v = pred[i, "x1"],
                 h = pred[i, "x2"],
                 lty = 2)
    panel.points(pred[i, "x1"],
                 pred[i, "x2"],
                 col = "green",
                 cex = 2,
                 pch = 19)
    })

#-----------------------------------------------------------------------
