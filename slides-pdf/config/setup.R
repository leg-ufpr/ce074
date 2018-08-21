library(knitr)
opts_chunk$set(dev.args = list(family = "Palatino"),
               fig.align = "center",
               fig.dim = c(5, 4))
# thm <- knit_theme$get("dusk")
thm <- knit_theme$get("solarized-light")
knit_theme$set(thm)

library(xtable)

print_matrix <- function(matrix, ...) {
    x <- print.xtable(xtable(x = matrix, ...),
                      only.contents = TRUE,
                      hline.after = NULL,
                      include.rownames = FALSE,
                      include.colnames = FALSE,
                      comment = FALSE,
                      type = "latex",
                      print.results = FALSE,
                      sanitize.text.function = identity)
    env <- "bmatrix"
    sprintf("\\begin{%s}\n%s\\end{%s}\n", env, x, env)
}
