library(plyr)
library(ggplot2)

setwd("~/Documents/projects/cluster/")
source("R/qqplots.R")

ens_uniprot <- unique(slr_structure_all[, c("stable_id", "uniprot_id")])
ec_numbers$uniprot_id <- ec_numbers$accession

logunif <- function (p, min = 0, max = 1, lower.tail = TRUE, log.p = FALSE) {
  log10(qunif(p, min, max, lower.tail, log.p))
}


## For the bar plots
make.fractions.ec <- function(master.df, variables, thr) {
  ddply(master.df, .variables=variables, function(df) {
    fraction <- sum(df$adj.pval < thr) / nrow(df)
    data.frame(number=nrow(df), Fraction=fraction)
  })
}

make.fractions.str <- function(master.df, variables, thr) {
  ddply(master.df, .variables=variables, function(df) {
    fraction <- sum(df$Adj.Pval < thr & df$omega > 1) / nrow(df)
    data.frame(number=nrow(df), Fraction=fraction)
  })
}

# slr_structure_all$Adj.Pval <- p.adjust(slr_structure_all$Pval, method="BH")
slr_structure_all <- join(slr_structure_all, ec_numbers)
slr_structure_all$ec_1 <- substr(slr_structure_all$ec_number, 1, 1)


# Cluster results for all thrs 
thrs <- c("0.05", "0.1", "0.2", "0.5")
# thrs <- c("0.2")
cluster_results_collated <- c()
for (thr in thrs) {
  print(thr)
  cluster_results_tmp <- read.table(paste0("data/clusters_", thr, ".tab"), sep="\t", header=F)
  colnames(cluster_results_tmp) <- c("stable_id", "pdb_id", "pdb_chain", "N", "n_sign", "pval")
  cluster_results_tmp$adj.pval <- p.adjust(cluster_results_tmp$pval, method="BH")
  cluster_results_tmp$thr <- thr
  cluster_results_collated <- rbind(cluster_results_collated, cluster_results_tmp)
}

cluster_results_collated_uniprot <- join(cluster_results_collated, ens_uniprot, type="inner", match="all")
which(!(cluster_results_collated_uniprot$stable_id %in% cluster_results_collated$stable_id))
cluster_results_stats <- ddply(cluster_results_collated_uniprot, "stable_id", function(df) {
  data.frame(N=length(unique(df$uniprot_id)))
})
subset(cluster_results_stats, N > 1)
subset(cluster_results_collated_uniprot, stable_id == "ENSP00000353099")


cluster_results_collated_ec_all <- join(cluster_results_collated_uniprot, ec_numbers, type="left")
cluster_results_collated_ec_all$ec_1 <- as.numeric(substr(cluster_results_collated_ec_all$ec_number, 1, 1))

ec_names <- c("Oxidoreductase", "Transferase", "Hydrolase", "Lyase", "Isomerase", "Ligase", "Mixed")
ec_cols <- c("red", "springgreen3", "blue", "brown", "darkgreen", "lightblue")

# cluster_results_stats <- ddply(cluster_results_collated_ec_all, "stable_id", function(df) {
#   data.frame(N=length(unique(df$ec_number)))
# })

# Alternative: how many different enzyme classes per protein?
cluster_results_stats <- ddply(cluster_results_collated_ec_all, "stable_id", function(df) {
  data.frame(N=length(unique(df$ec_1)))
})

# TODO could further distinguish between enzymes in multiple different classes 
# or with multiple annotations within a class
multi_ec <- unique(subset(cluster_results_stats, N > 1)$stable_id)
cluster_results_collated_multi_ec <- subset(cluster_results_collated_ec_all, stable_id %in% multi_ec)

cluster_results_collated_ec_filtered <- ddply(cluster_results_collated_ec_all, c("stable_id", "pdb_id", "pdb_chain", "thr"), function(df) {
  if(length(unique(df$ec_1)) > 1) {
    df$ec_1 <- 7
  }
  
  df[1, ]
})

cluster_results_collated_ec_filtered$type <- factor(ec_names[cluster_results_collated_ec_filtered$ec_1],
                                                       levels=c("Oxidoreductase", "Transferase", "Hydrolase", "Lyase", "Isomerase", "Ligase", "Mixed"))

nrow(cluster_results_collated_ec_filtered)
head(cluster_results_collated_ec_filtered)
subset(cluster_results_collated_ec_filtered, ec_1 == 7)

write.table(cluster_results_collated_ec_filtered, file="data/cluster_results_filtered.tab", sep="\t", quote=F, row.names=F)

cluster_results_collated_single_ec <- subset(cluster_results_collated_ec_all, !is.na(ec_number) & !(stable_id %in% multi_ec))
cluster_results_collated_single_no_ec <- subset(cluster_results_collated_ec_all, is.na(ec_number))
sum(c(nrow(cluster_results_collated_multi_ec), nrow(cluster_results_collated_single_ec), nrow(cluster_results_collated_single_no_ec)))
nrow(cluster_results_collated_ec_all)
# Could process and re-join the multi-EC results to the single_ec table



###   ^^^^ Cleaned up up to here ^^^^

cath_ann <- read.table("data/cath_ann.tab", sep="\t", stringsAsFactors=F, header=T)
cath_ann_short <- cath_ann[, c("pdb_id", "pdb_chain", "domain_id", "CATHCODE", "CLASS", "ARCH")]

cluster_results_collated_cath <- join(cluster_results_collated_ec_filtered, cath_ann_short, match = "first")
write.table(cluster_results_ann, file="data/cluster_results_ann.tab", sep="\t", quote=F, row.names=F)

cluster_results_collated_cath <- adply(cluster_results_collated_cath, 1, function(r) {
  f <- strsplit(r$ec_number[1], "\\.")[[1]]
  
  if (!is.na(f[1])) {
    data.frame(ec_1=as.numeric(f[1]), ec_2=as.numeric(f[2]), ec_3=as.numeric(f[3]), ec_4=as.numeric(f[4]))
  } else {
    data.frame(ec_1=NA, ec_2=NA, ec_3=NA, ec_4=NA)
  }
})

curr_thr <- 0.2
cluster_results <- subset(cluster_results_collated_cath, thr == curr_thr)
cluster_results_all <- subset(cluster_results_collated_ec_filtered, thr == curr_thr)
sum(cluster_results$adj.pval < 0.05)

cluster_results[order(cluster_results[, "adj.pval"]), ][1:30, ]

ggd.qqplot(pmax(0.000000001, cluster_results$pval))

plot.vertical.hist <- function(data, breaks=500) {
  hs <- hist(data, breaks=breaks, plot=FALSE)
  
  old.par <- par(no.readonly=TRUE)
  mar.default <- par('mar')
  mar.left <- mar.default
  mar.right <- mar.default
  mar.left[4] <- 0
  mar.right[2] <- 0
  
  # par (fig=c(0.8,1.0,0.0,1.0), mar=mar.right, new=TRUE)
  plot (NA, axes=FALSE, xaxt='n', yaxt='n', # type='n', 
        xlab=NA, ylab=NA, main=NA,
        xlim=c(0,max(hs$counts)),
        ylim=c(1,length(hs$counts)))
  # axis (1)
  arrows(rep(0,length(hs$counts)),1:length(hs$counts),
         hs$counts,1:length(hs$counts),
         length=0,angle=0)
  
  par(old.par)
  # invisible ()
}

plot.vertical.hist.rev <- function(data, breaks=500) {
  hs <- hist(data, breaks=breaks, plot=FALSE)
  
  old.par <- par(no.readonly=TRUE)
  mar.default <- par('mar')
  mar.left <- mar.default
  mar.right <- mar.default
  mar.left[4] <- 0
  mar.right[2] <- 0
  
  # par (fig=c(0.8,1.0,0.0,1.0), mar=mar.right, new=TRUE)
  plot (NA, axes=FALSE, xaxt='n', yaxt='n', # type='n', 
        xlab=NA, ylab=NA, main=NA,
        xlim=c(-max(hs$counts), 0),
        ylim=c(1,length(hs$counts)))
  # axis (1)
  arrows(rep(0,length(hs$counts)),1:length(hs$counts),
         -hs$counts,1:length(hs$counts),
         length=0,angle=0)
  
  par(old.par)
  # invisible ()
}

pdf("marginal_0.05.pdf", width=2, height=5)
plot.vertical.hist.rev(-log10(pmax(0.000000001, subset(cluster_results_collated, thr == 0.05)$pval)))
dev.off()

pdf("marginal_0.2.pdf", width=2, height=5)
plot.vertical.hist.rev(-log10(pmax(0.000000001, subset(cluster_results_collated, thr == 0.2)$pval)))
dev.off()

pdf("marginal_0.1.pdf", width=2, height=5)
plot.vertical.hist(-log10(pmax(0.000000001, subset(cluster_results_collated, thr == 0.1)$pval)))
dev.off()

pdf("marginal_0.5.pdf", width=2, height=5)
plot.vertical.hist(-log10(pmax(0.000000001, subset(cluster_results_collated, thr == 0.5)$pval)))
dev.off()


ggplot(cluster_results, aes(sample=pval)) +
  stat_qq(distribution=qunif)

hist(cluster_results$pval, breaks=1000)

site_sums <- ddply(slr_structure_all_strict, c("stable_id", "pdb_id", "pdb_chain"), function(df) {s
  data.frame(N=nrow(df), n=nrow(subset(df, Adj.Pval < curr_thr & omega > 1)))
})

site_sums_cluster <- subset(site_sums, n >= 2)
nrow(site_sums_cluster)
length(unique(site_sums_cluster$pdb_id))

ddply(cluster_results, "ec_1", function(df) {
  data.frame(frac=sum(df$adj.pval < 0.05)/nrow(df), mean_n=mean(df$n), N=nrow(df))
})

cluster_results_ec[order(cluster_results_ec[, "adj.pval"]), ][1:60, ]

csa <- read.table("~/Documents/projects/slr_pipeline/data/CSA_2_0_121113.txt", sep=",", stringsAsFactors=FALSE, header=TRUE)
nrow(csa)
head(csa)

csa <- subset(csa, PDB.ID %in% unique(slr_structure_all$pdb_id))
length(unique(csa$PDB.ID))
sum(unique(csa$PDB.ID) %in% subset(cluster_results_ann, n_sign >= 2)$pdb_id)

plot(cluster_results$n, -log10(cluster_results_ec$adj.pval), col=cluster_results_ec$ec_1)

pdf("results/qqplot.pdf", height=7, width=8)
ggd.qqplot(pmax(0.00000001, cluster_results$pval))
hist(-log10(cluster_results$pval), breaks=50, col="black")
dev.off()

pdf("results/EC_qqplots.pdf", height=7, width=8)
old_par <- par(mfrow=c(3, 2))
for (i in seq(1, 6))
  ggd.qqplot(pmax(0.000000001, subset(cluster_results, ec_1 == i)$pval), main=ec_names[i])
par(old_par)
dev.off()

## Plotting starts here
## Bar plots
# Numbers and sizes
ddply(slr_structure_all, c("ec_1"), function(df){
  stats <- ddply(df, c("stable_id", "pdb_id", "pdb_chain"), function(df2) {
    if (nrow(df2) > 50)
      data.frame(N=nrow(df2))
    else
      data.frame()
  })
  
  print(nrow(stats))
  print(summary(stats$N))
})


ggplot(subset(cluster_results_collated, thr==0.2), aes(x=factor(ec_1), y=N, fill=thr)) +
  geom_boxplot()

ggplot(cluster_results, aes(x=factor(ec_1), y=n_sign, colour=thr)) +
  theme_bw() +
  geom_boxplot() + 
  geom_jitter()


# Fraction of sites under pos sel.
str_fractions <- c()
for (thr in thrs) {
  str_fractions_tmp <- make.fractions.str(slr_structure_all, c("ec_1"), thr)
  str_fractions_tmp$thr <- thr
  str_fractions <- rbind(str_fractions, str_fractions_tmp)
}

ggplot(str_fractions, aes(x=factor(ec_1), y=Fraction, fill=thr)) +
  geom_bar(stat="identity", position=position_dodge()) +
  ggtitle("Fraction of positive selection") +
  theme_bw(base_size=18)

## Group by enzyme class
# Sample sizes
make.fractions.ec(slr_structure_all, c("ec_1"), 0.05)

ec_fractions <- make.fractions.ec(cluster_results_collated, c("ec_1", "thr"), 0.05)

ggplot(ec_fractions, aes(x=factor(ec_1), y=Fraction, fill=thr)) +
  geom_bar(stat="identity", position=position_dodge()) +
  ggtitle("Fraction of positive selection") +
  theme_bw(base_size=18)



ggplot(cluster_results_collated, aes(factor(ec_1), n_sign, fill=factor(thr))) +
  geom_boxplot()

# Are p-values correlated with n_sign?
ggplot(cluster_results_ann, aes(n_sign, -log10(pval), colour=factor(ec_1))) +
  geom_point()

cor.test(cluster_results_ann$n_sign, -log10(pmax(0.00001, cluster_results_ann$pval)))


## QQ plots 
clust_results_enzymes <- subset(cluster_results, ec_1 %in% c(1, 2, 3, 4, 5, 6))
pdf("results/enzyme_qqplot.pdf", height=7, width=8)
# png("results/enzyme_qqplot.png", height=500, width=600)
ggd.qqplot.mult(pmax(0.0000001, clust_results_enzymes$pval), as.factor(clust_results_enzymes$ec_1), cols=ec_cols, all=T)
legend("topleft", c(ec_names, "All"), fill=c(ec_cols, "black"), cex=0.75)
dev.off()

cluster_results_2 <- subset(clust_results_enzymes, n_sign <= 20)
ggd.qqplot.mult(pmax(0.000001, cluster_results_2$pval), cluster_results_2$ec_1, cols=ec_cols, all=T)
legend("topleft", c(ec_names, "All"), fill=c(ec_cols, "black"), cex=0.75)


## Why are lyases not clustered?
cluster_results_lyases <- subset(clust_results_enzymes, ec_1 == 4)
cluster_results_lyases <- cluster_results_lyases[order(cluster_results_lyases$pval), ]
nrow(cluster_results_lyases)
head(cluster_results_lyases)

## Non-enzymes
cluster_results_rest <- subset(cluster_results_ann, is.na(ec_number))
nrow(cluster_results_rest)
cluster_results_rest <- cluster_results_rest[order(cluster_results_rest$adj.pval, decreasing=F), ]
head(cluster_results_rest)

## Misc
cluster_subset <- subset(cluster_results_ann, ARCH == "3-Layer(aba) Sandwich")
ggd.qqplot(pmax(0.000000001, cluster_subset$pval))


## cluster_results_collated_multi_ec
multi_ec_unique <- data.frame(stable_id=multi_ec)
cluster_results_collated_multi <- subset(cluster_results_collated_multi_ec, thr == curr_thr)
cluster_results_collated_multi <- join(cluster_results_collated_multi, multi_ec_unique, type="right", match="first")
ggd.qqplot(pmax(0.000000001, cluster_results_collated_multi$pval))
head(cluster_results_collated_multi)
cluster_results_collated_multi[order(cluster_results_collated_multi[, "adj.pval"]), ][1:20, ]
