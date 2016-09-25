## Sort out duplicate datasets


# First work out a common limit
slr_structure_all_strict_sums <- ddply(slr_structure_all_strict, c("stable_id", "pdb_id", "pdb_chain"), function(df) {
  data.frame(nsign=sum(df$Adj.Pval < 0.5 & df$omega > 1))
})
x_lim <- max(slr_structure_all_strict_sums$nsign)
hists_thr <- list()

## Counts
# 0.05  118
# 0.10  137
# 0.20  177
# 0.50  253

for (thr in c(0.05, 0.1, 0.2, 0.5)) {
  slr_structure_all_strict_sums <- ddply(slr_structure_all_strict, c("stable_id", "pdb_id", "pdb_chain"), function(df) {
    data.frame(nsign=sum(df$Adj.Pval < thr & df$omega > 1))
  })
  
  # sum(slr_structure_all_strict_sums$nsign > 1)
  slr_structure_all_strict_sums_subset <- subset(slr_structure_all_strict_sums, nsign > 1)
  
  # This ends up being done with cowplot
  p <- ggplot(slr_structure_all_strict_sums_subset, aes(x=nsign)) +
    scale_x_continuous(limits=c(1, x_lim+1), breaks=seq(2, x_lim, 4)) +
    scale_y_continuous(limits=c(0, 90)) +
    theme_cowplot() +
    theme(axis.title=element_text(size=8), axis.text=element_text(size=8),
          plot.title=element_text(size=10)) +
    ggtitle(paste0("Positively-selected sites, (p < ", thr, ")")) +
    xlab("Number of positively selected sites") +
    ylab("Count") +
    geom_histogram(binwidth=1)
  
  hists_thr[[as.character(thr)]] <- ggplotGrob(p)
  print(p)
  ggsave(file=paste0("results/3D/hist_nsites_", thr, ".pdf"), width=9, height=7, unit="cm")
}

plot_grid(plotlist=hists_thr, nrow=2, ncol=2, align="hv")
ggsave(filename="results/3D/hists_nsites_all.pdf", width=16, height=14, unit="cm")


for (thr in c(0.05, 0.1, 0.2, 0.5)) {
  slr_structure_all_strict_sums <- ddply(slr_structure_all_strict, c("stable_id", "pdb_id", "pdb_chain"), function(df) {
    data.frame(nsign=sum(df$Adj.Pval < thr & df$omega > 1))
  })
  
  print(sum(slr_structure_all_strict_sums$nsign > 1))
}

read.results <- function(indir, method, thr) {
  filename <- paste0(indir, "/", "cluster_", thr, "_", method, ".tab")
  print(filename)
  
  results <- read.table(filename, sep="\t", header=F, stringsAsFactors=F, quote=NULL)
  
  if (method == "gr") { # Need to split the two submethods
    colnames(results) <- c("stable_id", "pdb_id", "pdb_chain", "n_raw", "n_mapped", "method", "cluster", "pval")
  } else if (method == "clumps") { # Nothing special to do
    colnames(results) <- c("stable_id", "pdb_id", "pdb_chain", "n_raw", "n_mapped", "cluster", "pval")
  } else if (method == "cucala") { # There may be multiple clusters
    colnames(results) <- c("stable_id", "pdb_id", "chain_id", "n_raw", "n_mapped", "cluster", "pval")
    results <- ddply(results, c("stable_id", "pdb_id", "chain_id", "n_raw", "n_mapped"), function(df) {
      data.frame(cluster=df$cluster, pval=df$pval, cluster_id=1:nrow(df))
    })
  } else {
    stop(paste("Unknown method", method))
  }
  
  results
}

results <- read.results("results", method="clumps", thr=0.5)
clusters_all <- data.frame(stable_id=unique(results$stable_id), stringsAsFactors=F)

for(method in c("clumps", "cucala", "gr")) {
  for (thr in c(0.05, 0.1, 0.2, 0.5)) {
    if (method == "gr") {
      results <- read.results("results", method=method, thr=thr)

      for (subm in c("gr_n", "gr_max")) {
        res <- subset(results, method == subm)
        clusters_curr <- join(clusters_all, res, type="left", match="first")
        clusters_all[[paste(subm, thr, sep="_")]] <- clusters_curr$pval  
      }
    } else {
      results <- read.results("results", method=method, thr=thr)
      clusters_curr <- join(clusters_all, results, type="left", match="first")
      clusters_all[[paste(method, thr, sep="_")]] <- clusters_curr$pval
    }
  }
}

clusters_dataset <- join(clusters_all, slr_all, match="first")
datasets <- data.frame(dataset=unique(clusters_dataset$dataset), stringsAsFactors=F)
clusters_all <- join(datasets, clusters_dataset, match="first")


plot(-log10(clusters_all$gr_n_0.2), -log10(clusters_all$gr_max_0.2))

plot(-log10(clusters_all$cucala_0.2), -log10(clusters_all$gr_max_0.2))

clusters_0.2 <- subset(clusters_all, !is.na(clumps_0.2))
clusters_0.2$clumps_0.2_adj <- p.adjust(clusters_0.2$clumps_0.2, method="BH")
write.table(subset(clusters_0.2, clumps_0.2_adj < 0.05)$stable_id, quote=F, row.names=F)

clusters_0.2 <- subset(clusters_all, !is.na(cucala_0.2))
clusters_0.2$cucala_0.2_adj <- p.adjust(clusters_0.2$cucala_0.2, method="BH")
write.table(subset(clusters_0.2, cucala_0.2_adj < 0.05)$stable_id, quote=F, row.names=F)

clusters_0.2 <- subset(clusters_all, !is.na(gr_n_0.2))
clusters_0.2$gr_n_0.2_adj <- p.adjust(clusters_0.2$gr_n_0.2, method="BH")
write.table(subset(clusters_0.2, gr_n_0.2_adj < 0.05)$stable_id, quote=F, row.names=F)

library(VennDiagram)

venn.diagram(filename="results/3D/venn.png", imagetype="png",
             list(cucala=subset(clusters_0.2, cucala_0.2_adj < 0.05)$stable_id,
             clumps=subset(clusters_0.2, clumps_0.2_adj < 0.05)$stable_id,
             graph=subset(clusters_0.2, gr_n_0.2_adj < 0.05)$stable_id))

# All three ENSP00000360372 ENSP00000260682 ENSP00000360317 ENSP00000353099
## What to do about multiple clusters from Cucala

## Some speculation as to which method is best
# Probably the easiest if we just quickly compare 

# Venn diagram for the methods

# Zoom in on those found by one method but not the other
# Likely the restriction of Cucala for the clusters to be spheres
# will make it screw up in some circumstances 
# Perhaps it is to be expected for CLUMPS to have the highest power 
# As it takes into account the most information: all residues are considered, the distances are taken into account

## Are their accessory sites? I.e. less significant sites 