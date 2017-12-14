context("AUC per region")

# pick an example gene
g = "SAMD11"

# pick example
s_names = c("A549_BR1_Baseline","A549_BR2_Baseline")
s_bigwigs = 


l = "C004GD_mature neutrophil"
coords = "chr1  1012923 1013923"
cmd = "WiggleTools/wiggletools apply AUC test.bed ~/links/CB.HTS.Analysis/CTTV020/data/BLUEPRINT/C004GDH1.ERX207032.H3K27ac.bwa.GRCh38.20150528.bw"

test_chip = make_auc_matrix("data/data_blueprint_parsed_test.csv", roi, "H3K27ac", "tmp/", quantile_norm=FALSE)

blueprint_chip[[1]]$res[which(blueprint_chip[[1]]$annot$Label==l), which(roi$gene==g)]
blueprint_chip_cut[[1]]$res[which(blueprint_chip_cut[[1]]$annot$Label==l), which(roi$gene==g)]
blueprint_chip_cut_final[[1]]$res[which(blueprint_chip_cut_final[[1]]$annot$Label==l), which(roi$gene==g)]

i = which(paste(dat$donor, dat$cell_type, sep="_")==l)
y = dat$data[i]
filter(y[[1]], gene==g) %>% dplyr::select(chip_k27ac) %>% exp()

# chip data looks good, so mapping from res_1 to res_2 should be tight...
# start from all_dat
old_dat = group_by(all_dat, donor, cell_type, group) %>% nest() # group by donor / cell type
j = which(paste(old_dat$donor, old_dat$cell_type, sep="_")==l)
old_dat_1 = old_dat$data[[j]] # pick first donor / cell type
x1 = arrange(old_dat_1, gene)
y1 = arrange(y[[1]], gene)
plot(x1$chip_k27ac, y1$chip_k27ac)