#Load files
stage1 <-as.data.frame(t(readRDS("stage1_df.rds")))
stage2 <-as.data.frame(t(readRDS("stage2_df.rds")))
stage3 <-as.data.frame(t(readRDS("stage3_df.rds")))
stage4 <-as.data.frame(t(readRDS("stage4_df.rds")))
stage5 <-as.data.frame(t(readRDS("stage5_df.rds")))

#Load library
library(bnlearn)
library(gplots)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Rgraphviz")
library(viridis)
#Truncate dataset
stage1_trunc = stage1[,c(1:50)]
stage2_trunc = stage2[,c(1:50)]
stage3_trunc = stage3[,c(1:50)]
stage4_trunc = stage4[,c(1:50)]
stage5_trunc = stage5[,c(1:50)]


#Learn DAG initial
dag1 = hc(stage1_trunc)
dag2 = hc(stage2_trunc)
dag3 = hc(stage3)
dag4 = hc(stage4)
dag5 = hc(stage5)

#Show DAGs
graphviz.plot(dag1)
graphviz.plot(dag2)


#Bootstrap DAGs
str.stage1 = boot.strength(stage1_trunc, algorithm = "hc")
plot(str.stage1)
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

avg.stage1 = averaged.network(str.stage1)
graphviz.plot(avg.stage1)

# stage 1
stage1_boot <- readRDS("stage1_bootstrap.rds")
avg.stage1 = averaged.network(stage1_boot, threshold = 0.75)
graphviz.plot(avg.stage1, main = "Stage 1", sub = "Threshold=0.75")

# stage 2
str.stage2 = boot.strength(stage2_trunc, algorithm = "hc")
avg.stage2_075 = averaged.network(str.stage2, threshold = 0.75)
graphviz.plot(avg.stage2_075, main = "Stage 2",sub = "Threshold = 0.75")
saveRDS(str.stage2, file = "stage2_bootstrap_075.rds")

# stage 3 
str.stage3 = boot.strength(stage3_trunc, algorithm = "hc")
avg.stage3 = averaged.network(str.stage3, threshold = 0.75)
graphviz.plot(avg.stage3, main = "Stage 3", sub = "Threshold=0.75")
saveRDS(str.stage3, file = "stage3_bootstrap.rds")

# stage 4
str.stage4 = boot.strength(stage4_trunc, algorithm = "hc")
saveRDS(str.stage4, file = "stage4_bootstrap.rds")
avg.stage4 <- averaged.network(str.stage4, threshold = 0.75)
graphviz.plot(avg.stage4, main = "Stage 4", sub = "Threshold=0.75")

# stage 5 
stage5_boot <- readRDS("stage5_bootstrap.rds")
avg.stage5 <- averaged.network(stage5_boot, threshold = 0.75)
graphviz.plot(avg.stage5, main = "Stage 5", sub= "Threshold = 0.75")


h12 = hamming(avg.stage1, avg.stage2)
h13 = hamming(avg.stage1, avg.stage3)
h14 = hamming(avg.stage1, avg.stage4)
h15 = hamming(avg.stage1, avg.stage5)
h23 = hamming(avg.stage2, avg.stage3)
h24 = hamming(avg.stage2, avg.stage4)
h25 = hamming(avg.stage2, avg.stage5)
h34 = hamming(avg.stage3, avg.stage4)
h35 = hamming(avg.stage3, avg.stage5)
h45 = hamming(avg.stage4, avg.stage5)


r1 = matrix(c(0, h12, h13, h14, h15), ncol=5)
r2 = matrix(c(h12, 0, h23, h24, h25), ncol=5)
r3 = matrix(c(h13, h23, 0, h34, h35), ncol=5)
r4 = matrix(c(h14, h24, h34, 0, h45), ncol=5)
r5 = matrix(c(h15, h25, h35, h45, 0), ncol=5)

dist_matrix = rbind(r1, r2, r3, r4, r5)
rownames(dist_matrix) <- c("Stage 1", "Stage 2", "Stage 3", "Stage 4", "Stage 5")
colnames(dist_matrix) <- c("Stage 1", "Stage 2", "Stage 3", "Stage 4", "Stage 5")
graphviz.compare(avg.stage3, avg.stage4)

graphviz.compare(avg.stage1, avg.stage5)
edges1v2 <- compare(avg.stage1, avg.stage2, arcs = TRUE)
edges2v3 <- compare(avg.stage2, avg.stage3, arcs = TRUE)
edges3v4 <- compare(avg.stage3, avg.stage4, arcs = TRUE)
edges4v5 <- compare(avg.stage4, avg.stage5, arcs = TRUE)
compare(avg.stage2, avg.stage3)
saveRDS(edges1v2, "edges1v2.rds")
saveRDS(edges2v3, "edges2v3.rds")
saveRDS(edges3v4, "edges3v4.rds")
saveRDS(edges4v5, "edges4v5.rds")

# visualize distances 
library(pheatmap)
pheatmap(dist_matrix, main="Heatmap of stage-wise hamming distances", color=viridis(10), display_numbers = T, cluster_rows = F, cluster_cols = F, fontsize_number = 15)


# source: >=1 children, no parents
# sink: no children, >= 1 parent

stage_boot <- readRDS("stage5_bootstrap.rds")
avg.stage = averaged.network(stage_boot, threshold = 0.75)

stage <- data.frame(nodes(avg.stage))
rownames(stage) <- stage$nodes.avg.stage
namevector <- c("in_deg", "out_deg")
stage[ , namevector] <- NA
names(stage)[1] <- "nodes"


for (node in nodes(avg.stage))
{
  stage[node, "in_deg"] <- in.degree(avg.stage, node)
  stage[node, "out_deg"] <- out.degree(avg.stage, node)
}


stage_srce <- stage[stage$in_deg==0, ]
stage_srce <- stage_srce[stage_srce$out_deg>=1,]

stage_sink <- stage[stage$out_deg==0, ]
stage_sink <- stage_sink[stage_sink$in_deg>=1,]

write.csv(stage_srce,"stage5_srce.csv", row.names = FALSE)
write.csv(stage_sink, "stage5_sink.csv", row.names = FALSE)


stage_boot <- readRDS("stage1_bootstrap.rds")
avg.stage1 = averaged.network(stage_boot, threshold = 0.75)
stage_boot <- readRDS("stage2_bootstrap.rds")
avg.stage2 = averaged.network(stage_boot, threshold = 0.75)
stage_boot <- readRDS("stage3_bootstrap.rds")
avg.stage3 = averaged.network(stage_boot, threshold = 0.75)
stage_boot <- readRDS("stage4_bootstrap.rds")
avg.stage4 = averaged.network(stage_boot, threshold = 0.75)
stage_boot <- readRDS("stage5_bootstrap.rds")
avg.stage5 = averaged.network(stage_boot, threshold = 0.75)
mat1 <- amat(avg.stage1)
mat2 <- amat(avg.stage2)
mat3 <- amat(avg.stage3)
mat4 <- amat(avg.stage4)
mat5 <- amat(avg.stage5)

mats <- mat1 + mat2 + mat3 + mat4 + mat5
uniques <- data.frame(which(mats==1, arr.ind=TRUE))
mat4[3,1]
dummy <- mat4[uniques$row, uniques$col]


uniques$stage1 <- NA
uniques$stage2 <- NA
uniques$stage3 <- NA
uniques$stage4 <- NA
uniques$stage5 <- NA

for (row in rownames(uniques)) {
  uniques[row, "stage1"] = mat1[uniques[row, "row"], uniques[row, "col"]]
  uniques[row, "stage2"] = mat2[uniques[row, "row"], uniques[row, "col"]]
  uniques[row, "stage3"] = mat3[uniques[row, "row"], uniques[row, "col"]]
  uniques[row, "stage4"] = mat4[uniques[row, "row"], uniques[row, "col"]]
  uniques[row, "stage5"] = mat5[uniques[row, "row"], uniques[row, "col"]]
}

length(mat5)
uniques$sanity = uniques$stage1 + uniques$stage2 + uniques$stage3 + uniques$stage4 + uniques$stage5
max(uniques$sanity)
min(uniques$sanity)
stage_uniques <- uniques[uniques$stage5==1,]
stage_uniques <- stage_uniques[, -c(3:8)]
stage_uniques$from <- rownames(mats)[stage_uniques$row]
stage_uniques$to <- colnames(mats)[stage_uniques$col]
write.csv(stage_uniques, "stage5_uniquegenes.csv")
write.csv(rownames(mats), "bg.csv")
