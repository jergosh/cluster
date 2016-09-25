slr_ref <- subset(slr_structure, stable_id == "ENSP00000316809")
meme_test <- read.table("159_1.txt", sep=",", header=T)[slr_ref$Site, ]
meme_test$human_idx <- slr_ref$human_idx
meme_test$pdb_pos <- slr_ref$pdb_pos

subset(meme_test, p.value < 0.1)
