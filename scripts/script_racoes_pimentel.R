##======================================================================
## Tipos de contrastes em modelos lineares
##======================================================================

##----------------------------------------------------------------------
## Carrega base de dados
## Disponível em:
# url <- "http://leg.ufpr.br/~fernandomayer/data/pimentel_racoes.txt"
da <- read.table("dados/pimentel_racoes.txt",
                 header = TRUE, sep = "\t")
str(da)
xtabs(~ racoes, da)

##----------------------------------------------------------------------
## Visualização
library(lattice)
xyplot(ganhopeso ~ racoes, data = da)
xyplot(ganhopeso ~ racoes, data = da, type = c("p", "a"))
bwplot(ganhopeso ~ racoes, data = da) # pouco útil pelas poucas repetições

##----------------------------------------------------------------------
## Ajuste do modelo
m0 <- lm(ganhopeso ~ racoes, data = da)

## ANOVA e efeitos
anova(m0)   # quadro de análise de variância
summary(m0) # quadro de estimativas dos parâmetros/efeitos
with(da, tapply(ganhopeso, racoes, mean))

##----------------------------------------------------------------------
## E SE... tivessem apenas 2 racoes? O que voce faria? De que forma?

## Faz o subset para ficar apenas com os tratamentos A e B
da2 <- subset(da, racoes %in% c("A", "B"))
## Remove níveis excluidos
da2 <- droplevels(da2)
str(da2)

## Teste t para duas amostras (nao pareado), assumindo que as variancias
## entre os dois grupos sao iguais (suposicao de homcedasticidade da
## ANOVA)
teste <- t.test(ganhopeso ~ racoes, data = da2, var.equal = TRUE)
teste

## Mesmo assim, podemos fazer uma ANOVA so com 2 niveis?
## Por que nao?
m.teste <- lm(ganhopeso ~ racoes, data = da2)
anova(m.teste)
## Qual a diferenca? Nenhuma!
## O teste t é um caso particular de uma ANOVA com um fator de 2
## niveis.
## Lembre também que a distribuição F pode ser obtida através do
## quadrado de uma distribuição t (ou uma t pode ser obtida através da
## raiz quadrada de uma F)
teste$statistic^2 # é o valor de F na tabela de ANOVA
## por isso os p-valores são iguais

## Note ainda que os efeitos estimados
summary(m.teste)
## são a média do tratamento A (intercepto) e a diferença entre a média
## do tratamento A com a média do tratamento B (inclinação)
with(da2, tapply(ganhopeso, racoes, mean))
diff(with(da2, tapply(ganhopeso, racoes, mean)))

## Veja que a parametrização da matriz X padrão usada pela função lm()
## do R é a de assumir o primeiro nível zero (chamada aqui de
## contr.treatment)
model.matrix(m.teste)

##======================================================================
## Ajuste na mao

##----------------------------------------------------------------------
## Sem usar nenhum tipo de restrição
n <- nrow(da2)
nn <- nlevels(da2$racoes)
y <- da2$ganhopeso
X <- matrix(0, nrow = n, ncol = nn)
X[cbind(seq_along(da2$racoes), da2$racoes)] <- 1
(X <- cbind(1, X))
## X'
Xt <- t(X)
## X'X
(Xt %*% X)
## (X'X)^-1
solve(Xt %*% X)
## Matriz de posto incompleto - não admite inversa
Matrix::rankMatrix(X)
Matrix::rankMatrix(Xt %*% X)

##----------------------------------------------------------------------
## Redefinindo o modelo: média de caselas (modelo sem intercepto)

## Guarda a matriz do modelo completo para facilitar
X0 <- X

## Matriz do modelo
(X <- X0[ , -1])
## X'
Xt <- t(X)
## X'X
(XtX <- Xt %*% X)
## (X'X)^-1
solve(XtX)
## X'y
(Xty <- Xt %*% y)
## (X'X)^-1 X'y
solve(XtX) %*% Xty
## Que é o mesmo que
with(da2, tapply(ganhopeso, racoes, mean))
## Usando lm
(m.medias <- lm(ganhopeso ~ racoes - 1, data = da2))
summary(m.medias)
model.matrix(m.medias)

##----------------------------------------------------------------------
## Restrição: Usando restrição soma zero

## Matriz do modelo
(X <- cbind(X0[, 1], X0[, 2] - X0[, 3]))
## X'
Xt <- t(X)
## X'X
(XtX <- Xt %*% X)
## (X'X)^-1
solve(XtX)
## X'y
(Xty <- Xt %*% y)
## (X'X)^-1 X'y
solve(XtX) %*% Xty
## Que é o mesmo que
mean(y) # média geral - intercepto
## O efeito é calculado como:
## média do primeiro menos a média geral, porque:
## \mu_1 = \mu + \alpha_1  =>  \alpha_1 = \mu_1 - \mu
mean(da2$ganhopeso[da2$racoes == "A"]) - mean(da2$ganhopeso)
## OU a média geral menos a média do segundo, porque:
## \mu_2 = \mu - \alpha_1  =>  \alpha_1 = \mu - \mu_2
mean(da2$ganhopeso) - mean(da2$ganhopeso[da2$racoes == "B"])
## Usando lm
(m.soma <- lm(ganhopeso ~ racoes, data = da2,
              contrasts = list(racoes = contr.sum)))
summary(m.soma)
model.matrix(m.soma)
## Novamente, as médias podem então ser obtidas com:
mean(da2$ganhopeso) # media geral = intercepto
mean(da2$ganhopeso) + coef(m.soma)[2] # mu_1 = mu + alpha_1
mean(da2$ganhopeso) - coef(m.soma)[2] # mu_2 = mu - alpha_1
with(da2, tapply(ganhopeso, racoes, mean))

## Usando aov() para acessar alguns métodos específicos
m.soma.aov <- aov(ganhopeso ~ racoes, data = da2)
model.tables(m.soma.aov, type = "means")
model.tables(m.soma.aov, type = "effects")

##----------------------------------------------------------------------
## Restrição: primeiro nível zero

## Matriz do modelo
(X <- X0[, -2])
## X'
Xt <- t(X)
## X'X
(XtX <- Xt %*% X)
## (X'X)^-1
solve(XtX)
## X'y
(Xty <- Xt %*% y)
## (X'X)^-1 X'y
solve(XtX) %*% Xty
## Que é o mesmo que
## Intercepto é a média do primeiro nível
with(da2, tapply(ganhopeso, racoes, mean))
## Inclinação é a diferença entre a média do segundo com a média do
## primeiro
diff(with(da2, tapply(ganhopeso, racoes, mean)))
## mu_2 = mu_1 + alpha_2  =>  alpha_2 = mu_2 - mu_1
mean(da2$ganhopeso[da2$racoes == "B"]) -
    mean(da2$ganhopeso[da2$racoes == "A"])
## Usando lm (já é o padrão do R - não precisa especificar o contraste)
(m.nivel1 <- lm(ganhopeso ~ racoes, data = da2))
summary(m.nivel1)
model.matrix(m.nivel1)

## Efeitos de Yates: diferença entre as médias dos níveis "alto" e
## "baixo" de um fator (usando o objeto criado acima com aov)
## Esse é o alpha_2 do contr.treatment
dae::yates.effects(m.soma.aov, data = da2)
## E esse é o alpha_1 do contr.sum
dae::yates.effects(m.soma.aov, data = da2)/2
## Ou seja: o efeito calculado no contr.sum é a diferença _média_ entre
## as médias dos dois níveis

##----------------------------------------------------------------------
## Restrição: último nível zero

## Matriz do modelo
(X <- X0[, -3])
## X'
Xt <- t(X)
## X'X
(XtX <- Xt %*% X)
## (X'X)^-1
solve(XtX)
## X'y
(Xty <- Xt %*% y)
## (X'X)^-1 X'y
solve(XtX) %*% Xty
## Que é o mesmo que
## Intercepto é a média do último nível
with(da2, tapply(ganhopeso, racoes, mean))
## Inclinação é a diferença entre a média do primeiro com a média do
## segundo
-diff(with(da2, tapply(ganhopeso, racoes, mean)))
## mu_1 = mu_2 + alpha_1  =>  alpha_1 = mu_1 - mu_2
mean(da2$ganhopeso[da2$racoes == "A"]) -
    mean(da2$ganhopeso[da2$racoes == "B"])
## Usando lm (já é o padrão do R - não precisa especificar o contraste)
(m.nivel2 <- lm(ganhopeso ~ racoes, data = da2,
                contrasts = list(racoes = contr.SAS)))
summary(m.nivel2)
model.matrix(m.nivel2)

##----------------------------------------------------------------------
## Forma geral de ver os contrastes
contrasts(C(da2$racoes, treatment)) # (default) primeiro nível tem efeito 0
contrasts(C(da2$racoes, SAS))       # o último nível tem efeito 0
contrasts(C(da2$racoes, sum))       # soma dos efeitos igual zero
options()$contrasts # aqui você verifica o que tá sendo usado por default
help(contrasts) # para criar o seu próprio contraste
