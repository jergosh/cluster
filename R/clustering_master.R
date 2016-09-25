library(plyr) 

cluster_master <- data.frame(stable_id=unique(slr_all$stable_id), stringsAsFactors=F)

cluster_master <- join(cluster_master, runs.all_0.05[, c("stable_id", "p")], match="first")
cluster_master$waw_0.05 <- cluster_master$p
cluster_master$p <- NULL

cluster_master <- join(cluster_master, runs.all_0.1[, c("stable_id", "p")], match="first")
cluster_master$waw_0.1 <- cluster_master$p
cluster_master$p <- NULL

cluster_master <- join(cluster_master, runs.all_0.2[, c("stable_id", "p")], match="first")
cluster_master$waw_0.2 <- cluster_master$p
cluster_master$p <- NULL

cluster_master <- join(cluster_master, runs.all_0.5[, c("stable_id", "p")], match="first")
cluster_master$waw_0.5 <- cluster_master$p
cluster_master$p <- NULL

head(cluster_master)
hist(subset(cluster_master, !is.nan(waw_0.2))$waw_0.2, breaks=50)
cluster_master$waw_0.2_adj[!is.nan(cluster_master$waw_0.2)] <-  p.adjust(cluster_master$waw_0.2[!is.na(cluster_master$waw_0.2)], method="BH")
sum(cluster_master$waw_0.2_adj < 0.05, na.rm=T)

sum(subset(cluster_master, !is.nan(waw_0.2)$ )

## Load other results