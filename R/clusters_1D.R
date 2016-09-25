library(plyr)
library(ggplot2)
library(grid)
library(gtable)

library(cowplot)

setwd("~/Documents/projects/cluster")

nrow(slr_structure_all_strict)


cluster_1d_continuous <- read.table("data/cluster_1d_multi_results.tab", header=F, sep="\t", stringsAsFactors=F)
head(cluster_1d_continuous)

colnames(cluster_1d_continuous) <- c("stable_id", "dataset", "maxI", "centre", "x", "start", "end", "pval", "cluster_id")
# Python has 0-based indexing
cluster_1d_continuous$start <- cluster_1d_continuous$start + 1
cluster_1d_continuous$end <- cluster_1d_continuous$end + 1

# The multiple testing correction can be done two different ways: correct only the p-value from the first cluster, that is for the evidence of clustering
# Use all the clusters, that is for a p-value associated with each cluster.
# OR by the multiplication of values

cluster_1d_continuous_first <- subset(cluster_1d_continuous, cluster_id == 1)
cluster_1d_continuous_first$adj.pval <- p.adjust(cluster_1d_continuous_first$pval, method="BH")
hist(cluster_1d_continuous_first$adj.pval)

cluster_1d_continuous$adj.pval <- p.adjust(cluster_1d_continuous$pval, method="BH")
hist(cluster_1d_continuous$adj.pval)
sum(cluster_1d_continuous$adj.pval < 0.05)

sign_clusters_1d_continuous <- ddply(cluster_1d_continuous, c("stable_id"), function(df) { data.frame(nsign=sum(df$adj.pval < 0.05)) })
hist(sign_clusters_1d_continuous$nsign)
sum(sign_clusters_1d_continuous$nsign >= 4)

ggplot(sign_clusters_1d_continuous, aes(x=nsign)) +
  xlab("Number of clusters") + 
  ylab("Count") +
  scale_x_continuous(breaks=0:6) +
  geom_bar() +
  theme_minimal(base_size = 12, base_family = "Helvetica") +
  theme(panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank())

ggsave(file="results/1D/clust_1D_continuous.pdf", width=12, height=10, unit="cm")

## Interesting case studies
feature_data <- read.table("~/Documents/projects/slr_pipeline/data/protein_features.tab", header=F, fill=T, sep="\t", quote=c())
colnames(feature_data) <- c("stable_id", "start", "end", "type", "intepro_id", "domain_id", "domain_desc")
head(feature_data)

plot.gene()
subset(sign_clusters_1d_continuous, nsign >= 4)

fix_coords <- function(clusters) {
  # The first cluster is correct.
  ret <- data.frame()
  
  offsets <- data.frame()
  
  for (i in 1:nrow(clusters)) {
    start <- clusters[i, "start"]
    end <- clusters[i, "end"]
    
    fixed_start <- start
    fixed_end <- end
    
    if (nrow(offsets))
    for (j in seq(1, nrow(offsets))) {
      offset_start <- offsets[j, "start"]
      len <- offsets[j, "len"]
      
      if (fixed_start >= offset_start) {
        fixed_start <- fixed_start + len
      }
      
      if (fixed_end >= offset_start) {
        fixed_end <- fixed_end + len
      }
    }
    
    
    offsets <- rbind(offsets, data.frame(start=fixed_start, len=(fixed_end-fixed_start+1)))
    # print(offsets)
    
    prev_start <- start
    prev_end <- end
    
    fixed_row <- clusters[i, ]
    fixed_row[1, "start"] <- fixed_start
    fixed_row[1, "end"] <- fixed_end
    ret <- rbind(ret, fixed_row)
  }
  
  ret
}

plot_example <- function(ex, st_id, plot_clusters=T) {
  ex_clusters <- subset(cluster_1d_continuous, stable_id == st_id & adj.pval < 0.05)
  print(ex_clusters)
  ex_clusters_fixed <- fix_coords(ex_clusters)
  print(ex_clusters_fixed)
  
  ex$xmin <- 1:(nrow(ex))
  ex$xmax <- ex$xmin + 1
  ex$ymin <- 0
  ex$ymax <- ex$Omega
  
  ex_clusters_fixed$ymin <- 0
  ex_clusters_fixed$ymax <- max(ex$Omega) + 0.2
  ex_clusters_fixed$end <- ex_clusters_fixed$end +1
  
  p <- ggplot(ex, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)) +
    theme_minimal(base_size = 8, base_family = "Helvetica") +
    theme(plot.margin=margin(l=0.0, unit="cm"),
          panel.margin=margin(l=0.0, r=0, unit="cm"),
          panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(),
          panel.border=element_blank(), panel.background = element_blank(),
          axis.text.y=element_text(margin=margin(l=0.0, r=0.0, unit="cm")),
          axis.ticks=element_line(size = 0.25),
          axis.ticks.x=element_line(), axis.ticks.y=element_line(),
          axis.ticks.length=unit(0.1, units="cm")) +
    scale_x_continuous(expand = c(0, 0)) + # limits=c(1, nrow(ex)+1), 
    scale_y_continuous(limits=c(0, max(ex$ymax+0.2)), expand = c(0, 0)) +
    ggtitle(st_id) +
    geom_rect()
  
  if (plot_clusters) {
    p <- p + geom_rect(data=ex_clusters_fixed, aes(xmin=start, xmax=end, ymin=ymin, ymax=ymax), fill="lightpink")
  }
  
  p + geom_rect()
}

## Actually
## Non-overlapping annotations shouldn't be so hard:
## Within each annotation type, we can check for overlaps
## and move the current annotation to the next row.

# But for now, let's assume a single, non-overlapping ann, type
plot.ann <- function(anns, len) {
  anns$ymin <- 0
  anns$ymax <- 1
  anns$textstart <- anns$start + (anns$end-anns$start)/2
  
  p <- ggplot(anns, aes(xmin=start, xmax=end, ymin=ymin, ymax=ymax, fill=domain_id)) +
      theme_minimal(base_size = 8, base_family = "Helvetica") +
      theme(plot.margin=margin(l=0.0, t=0.1, unit="cm"),
            panel.margin=margin(l=0.0, unit="cm"),
            panel.background=element_rect(fill="lightgrey", colour="white"),
            panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(),
            panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
            panel.border = element_blank(), panel.background = element_blank(),
            # axis.text.y=element_text(margin=margin(r=0, l=0.5)), panel.margin=unit(0, units="cm"),
            axis.title=element_blank(),
            axis.ticks.x=element_line(), axis.ticks.y=element_blank(),
            axis.text=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.ticks.length=unit(0.0, units="cm"),
            legend.position="none"
            ) +
      scale_x_continuous(limits=c(1, len+1), expand = c(0, 0)) +
      scale_y_continuous(limits=c(0, 1), expand = c(0, 0)) +
      geom_rect() +
      geom_text(aes(x=textstart, y=0.5, label=domain_desc), hjust=0.5, vjust=0.5, fontface=2, size=8*0.35)
  
  p
}

plot_all <- function(st_id, plot_clusters=T) {
  ex <- subset(slr_all, stable_id == st_id)
  ex_ann <- subset(feature_data, stable_id == st_id & type == "pfam")

  p_ex <- plot_example(ex, st_id, plot_cluster=plot_clusters)
  if (nrow(ex_ann)) {
    p_ann <- plot.ann(ex_ann, len=nrow(ex))
    plot_grid(p_ex, p_ann, ncol=1, align="v", rel_heights=c(1, 0.2))  
  } else {
    print(p_ex)
  }
}

save_plots <- function(st_id) {
  plot_all(st_id, plot_clusters=T)
  ggsave(file=paste("results/1D/bulk", paste0(st_id, ".pdf"), sep="/"), width=14, height=6, unit="cm")
  
  plot_all(st_id, plot_clusters=F)
  ggsave(file=paste("results/1D/bulk", paste0(st_id, "_nocluster", ".pdf"), sep="/"), width=14, height=6, unit="cm")  
}

ddply(subset(sign_clusters_1d_continuous, nsign >= 4), "stable_id", function(df) { save_plots(df$stable_id[1]) })

plot_all("ENSP00000263063", plot_clusters=T)
ggsave(file="results/1D/ex2_ENSP00000263063.pdf", width=14, height=6, unit="cm")

plot_all("ENSP00000263063", plot_clusters=F)
ggsave(file="results/1D/ex1_ENSP00000263063_nocluster.pdf", width=14, height=6, unit="cm")


plot_all("ENSP00000314348", plot_clusters=T)
ggsave(file="results/1D/ex2_ENSP00000314348.pdf", width=14, height=6, unit="cm")

plot_all("ENSP00000314348", plot_clusters=F)
ggsave(file="results/1D/ex2_ENSP00000314348_nocluster.pdf", width=14, height=6, unit="cm")



plot_all("ENSP00000280346", plot_clusters=T)
ggsave(file="results/1D/ex3_ENSP00000280346.pdf", width=14, height=6, unit="cm")


# ex_2 <- subset(slr_all, stable_id == "ENSP00000268919")
# ex_2_ann <- subset(feature_data, stable_id == "ENSP00000268919" & type == "pfam")
# 
# p_ex <- plot_example(ex_2, "ENSP00000268919")
# p_ann <- plot.ann(ex_2_ann, len=nrow(ex_2))
# 
# plot_grid(p_ex, p_ann, ncol=1, align="v", rel_heights=c(1, 0.2))
# ggsave(file="results/1D/ex2_ENSP00000268919.pdf", width=14, height=6, unit="cm")
# 
# p_ex <- plot_example(ex_2, "ENSP00000268919", plot_clusters=F)
# plot_grid(p_ex, p_ann, ncol=1, align="v", rel_heights=c(1, 0.2))
# ggsave(file="results/1D/ex1_ENSP00000268919_noclusters.pdf", width=14, height=6, unit="cm")

