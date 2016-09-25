

slr_structure_subset <- subset(slr_structure_all_strict, omega > 1 & Adj.Pval < 0.5)

colour.df <- ddply(slr_structure_subset, c("stable_id", "pdb_id"), function(df) {
  stable_id <- df$stable_id[1]
  pdb_id <- df$pdb_id[1]
  chain_id <- df$pdb_chain[1]
  out.df <- df[, c("pdb_id", "pdb_chain",	"pdb_pos", "Adj.Pval")]
  out.df$color <- "salmon"
  out.df$color[out.df$Adj.Pval < 0.2] <- "warmpink"
  out.df$color[out.df$Adj.Pval < 0.1] <- "tv_red"
  out.df$color[out.df$Adj.Pval < 0.05] <- "red"
  out.df <- out.df[, c("pdb_id", "pdb_chain",  "pdb_pos", "color")]
  
  out.df$group <- "positive"

  csa.residues <- unique(csa[csa$PDB.ID == pdb_id & csa$CHAIN.ID == chain_id, "RESIDUE.NUMBER"])
  if (length(csa.residues)) {
    csa.df <- data.frame(pdb_id=pdb_id, pdb_chain=chain_id, pdb_pos=csa.residues, color="orange", group="CSA")
    
    both.df <- join(out.df[, c("pdb_chain", "pdb_pos")], csa.df[, c("pdb_chain", "pdb_pos")], type="inner")
    print(paste(pdb_id, chain_id))
    
    out.df <- rbind(out.df, csa.df)  
  }
  write.table(out.df, file=paste0("data/color_strict/", pdb_id, '.tab'), row.names=F, quote=F, sep="\t")
  out.df
})

nrow(colour.df)
length(unique(colour.df$pdb_id))
