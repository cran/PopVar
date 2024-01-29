## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  # eval = FALSE,
  collapse = TRUE,
  comment = "#>"
)

## ----example, eval = TRUE-----------------------------------------------------
# Load the package
library(PopVar)

# Load the example data
data("think_barley", package = "PopVar")

## ----example1, message=FALSE, results='hide'----------------------------------
out <- pop.predict(G.in = G.in_ex, y.in = y.in_ex, map.in = map.in_ex,
                   crossing.table = cross.tab_ex,
                   nInd = 1000, nSim = 1, 
                   nCV.iter = 1, models = "rrBLUP")

## ----combine1-----------------------------------------------------------------
predictions1 <- lapply(X = out$predictions, FUN = function(x) {
  x1 <- as.data.frame(apply(X = x, MARGIN = 2, FUN = unlist), stringsAsFactors = FALSE)
  cbind(x1[,c("Par1", "Par2")], sapply(X = x1[,-1:-2], as.numeric)) 
})

# Display the first few lines of the predictions for grain yield
knitr::kable(head(predictions1$Yield_param.df))


## ----example2-----------------------------------------------------------------
out2 <- pop.predict2(G.in = G.in_ex_imputed, y.in = y.in_ex, map.in = map.in_ex,
                     crossing.table = cross.tab_ex, models = "rrBLUP")

knitr::kable(head(subset(out2, trait == "Yield")))


## ----example3-----------------------------------------------------------------
out3 <- pop_predict2(M = G.in_ex_mat, y.in = y.in_ex, map.in = map.in_ex,
                     crossing.table = cross.tab_ex, models = "rrBLUP")

knitr::kable(head(subset(out2, trait == "Yield")))


## ----compare1, message=FALSE--------------------------------------------------
time1 <- system.time({
  capture.output(pop.predict.out <- pop.predict(
    G.in = G.in_ex_imputed, y.in = y.in_ex, map.in = map.in_ex, crossing.table = cross.tab_ex,
    nInd = 1000, nSim = 1, nCV.iter = 1, models = "rrBLUP"))
})

time2 <- system.time({pop.predict2.out <- pop.predict2(
  G.in = G.in_ex_imputed, y.in = y.in_ex, map.in = map.in_ex,
  crossing.table = cross.tab_ex,model = "rrBLUP")})

# Print the time (seconds) required for each function.
c(pop.predict = time1[[3]], pop.predict2 = time2[[3]])


## ----compare2-----------------------------------------------------------------
predictions1 <- lapply(X = pop.predict.out$predictions, FUN = function(x) {
  x1 <- as.data.frame(apply(X = x, MARGIN = 2, FUN = unlist), stringsAsFactors = FALSE)
  cbind(x1[,c("Par1", "Par2")], sapply(X = x1[,-1:-2], as.numeric))
})

pop.predict.out1 <- predictions1$Yield_param.df[,c("Par1", "Par2", "pred.varG")]
pop.predict2.out1 <- subset(pop.predict2.out, trait == "Yield", c(parent1, parent2, pred_varG))

toplot <- merge(pop.predict.out1, pop.predict2.out1, by.x = c("Par1", "Par2"),
                by.y = c("parent1", "parent2"))

plot(pred.varG ~ pred_varG, toplot,
     xlab = "pop.predict2", ylab = "pop.predict",
     main = "Comparsion of predicted genetic variance")


## ----mp.example1--------------------------------------------------------------
# Generate predictions for all possible 4-way crosses of 10 sample parents
sample_parents <- sample(unique(unlist(cross.tab_ex)), 10)

mp_out <- mppop.predict(G.in = G.in_ex_imputed, y.in = y.in_ex, map.in = map.in_ex,
                        parents = sample_parents, n.parents = 4, models = "rrBLUP")

knitr::kable(head(subset(mp_out, trait == "Yield")))


## ----mp.example2--------------------------------------------------------------
# Generate predictions for all possible 4-way crosses of 10 sample parents
mp_out2 <- mppop_predict2(M = G.in_ex_mat, y.in = y.in_ex, map.in = map.in_ex,
                          parents = sample_parents, n.parents = 4, models = "rrBLUP")

knitr::kable(head(subset(mp_out2, trait == "Yield")))


