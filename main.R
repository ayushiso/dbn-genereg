#Load files
stage1 <-as.data.frame(t(readRDS("stage1_df.rds")))
stage2 <-as.data.frame(t(readRDS("stage2_df.rds")))
stage3 <-as.data.frame(t(readRDS("stage3_df.rds")))
stage4 <-as.data.frame(t(readRDS("stage4_df.rds")))
stage5 <-as.data.frame(t(readRDS("stage5_df.rds")))

#Load library
library(bnlearn)
library(gplots)

#Truncate dataset
stage1_trunc = stage1[,c(1:50)]
stage5_trunc = stage5[,c(1:50)]

#Learn DAG initial
dag1 = hc(stage1_trunc)
dag2 = hc(stage2)
dag3 = hc(stage3)
dag4 = hc(stage4)
dag5 = hc(stage5)

#Show DAGs
graphviz.plot(dag1)


#Bootstrap DAGs
str.stage1 = boot.strength(stage1_trunc, algorithm = "hc")
plot(str.stage1)
str.stage5 = boot.strength(stage5_trunc, algorithm = "hc")

#Cross Validating Thresholds
abline(v = 0.65, col = "green", lty = 2, lwd = 2)
abline(v = 0.75, col = "tomato", lty = 2, lwd = 2)
abline(v = 0.85, col = "steelblue", lty = 2, lwd = 2)
nrow(str.stage1[str.stage1$strength > attr(str.stage1, "threshold") &
                  +               str.stage1$direction > 0.5, ])
nrow(str.stage1[str.stage1$strength > 0.65 &
                  +               str.stage1$direction > 0.5, ])
nrow(str.stage1[str.stage1$strength > 0.75 &
                  +               str.stage1$direction > 0.5, ])
nrow(str.stage1[str.stage1$strength > 0.85 &
                  +               str.stage1$direction > 0.5, ])
avg.stage1_50 = averaged.network(str.stage1)
avg.stage1_65 = averaged.network(str.stage1,threshold = 0.65)
avg.stage1_75 = averaged.network(str.stage1,threshold = 0.75)
avg.stage1_85 = averaged.network(str.stage1,threshold = 0.85)
bn.cv(stage1_trunc, bn = avg.stage1_50, loss = "logl-g")
bn.cv(stage1_trunc, bn = avg.stage1_65, loss = "logl-g")
bn.cv(stage1_trunc, bn = avg.stage1_75, loss = "logl-g")
bn.cv(stage1_trunc, bn = avg.stage1_85, loss = "logl-g")

#Plotting graph with threshold of 0.75
avg.stage1 = averaged.network(str.stage1,threshold = 0.75)
graphviz.plot(avg.stage1,main = "Stage 1",sub = "Threshold = 0.75")
avg.stage5 = averaged.network(str.stage5,threshold = 0.75)
graphviz.plot(avg.stage5,main = "Stage 5",sub = "Threshold = 0.75")