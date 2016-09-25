library(seqinr)

kinase_fa <- read.fasta("data/uniprot_kinases.fa")
kinase_list <- substr(names(kinase_fa), 4, 9)
length(kinase_list)
sum(kinase_list %in% unique(slr_structure_all$uniprot_id))

slr_kinases <- subset(slr_structure_all, uniprot_id %in% kinase_list)
length(unique(slr_kinases$stable_id))
slr_kinases_pos <- subset(slr_kinases, Adj.Pval < 0.1 & omega > 1)
length(unique(slr_kinases_pos$stable_id))

subset(cluster_results, uniprot_id %in% kinase_list)

TF_fa <- read.fasta("data/Homo_sapiens_transcription_factors.fasta")
TF_list <- names(TF_fa)
length(TF_list)
sum(TF_list %in% unique(slr_structure_all$stable))

slr_TFs <- subset(slr_structure_all, stable_id %in% TF_list)
slr_TFs_pos <- subset(slr_TFs, Adj.Pval < 0.1 & omega > 1)
nrow(slr_TFs_pos)

subset(cluster_results, stable_id %in% TF_list)

# Try to isolate the immune stuff
# This is based on BP -- perhaps we would expect more from 'binding' MF
immune_list <- unique(read.table("data/immune_process.txt", header=T, sep="\t", stringsAsFactors=F)$Ensembl.Protein.ID)
innate_list <- unique(read.table("data/innate_immune_response.txt", header=T, sep="\t", stringsAsFactors=F)$Ensembl.Protein.ID)
adaptive_list <- unique(read.table("data/adaptive_immune_response.txt", header=T, sep="\t", stringsAsFactors=F)$Ensembl.Protein.ID)

list_subset <- function(id_list) {
  sum(id_list %in% unique(slr_structure_all$stable))
  slr_immunity <- subset(slr_structure_all, stable_id %in% id_list)
  print(nrow(slr_immunity))
  slr_immune_pos <- subset(slr_immunity, Adj.Pval < 0.1 & omega > 1)
  print(length(unique(slr_immune_pos$stable_id)))
  
  immune_clusters <- subset(cluster_results_all, stable_id %in% id_list)
  ggd.qqplot(pmax(0.000000001, immune_clusters$pval))
  print(immune_clusters[order(immune_clusters$pval, decreasing=F)[1:20], ])
}

list_subset(immune_list)
list_subset(innate_list)
list_subset(adaptive_list)

receptor_list <- unique(read.table("data/receptor_binding.txt", header=T, sep="\t", stringsAsFactors=F)$Ensembl.Protein.ID)
list_subset(receptor_list)
