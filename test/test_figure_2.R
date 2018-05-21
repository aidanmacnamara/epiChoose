# comparison of basal cell state - monocyte (4 conditions), thp1, u937

# 1. get the necessary files from blueprint

bp <- read_tsv("inst/extdata/blueprint_files.tsv")
bp[is.na(bp)] = "" # removes nas that affect dplyr selections

# what are the monocyte/macrophage treatment groups in bp?
table(grep("monocyte|macrophage", bp$`Sub-group`, value=TRUE))
trmt_groups = unique(grep("monocyte|macrophage", bp$`Sub-group`, value=TRUE))
bp_filt = bp %>% filter(`Sub-group` %in% trmt_groups, Format=="bigWig")

# which of these has the 3 epi marks + rna?
chip_assays = c("H3K27ac","H3K4me3","H3K27me3")
all_assays = c("H3K27ac","H3K4me3","H3K27me3","RNA-Seq")
group_vars = lapply(c("Group","Sub-group","Name","Sex","Cell type","Tissue","Cell line","Disease","Donor"), as.symbol)

bp_with_assays = bp_filt %>% group_by_(.dots=group_vars) %>% summarise(
  chip_assays = all(chip_assays %in% Experiment),
  all_assays = all(all_assays %in% Experiment),
  files = n(),
  H3K27ac = length(Experiment[Experiment=="H3K27ac"]),
  H3K4me3 = length(Experiment[Experiment=="H3K4me3"]),
  H3K27me3 = length(Experiment[Experiment=="H3K27me3"]),
  rna_seq = length(Experiment[Experiment=="RNA-Seq"])
)
table(bp_with_assays$`Sub-group`)
table(bp_with_assays$all_assays)
table(bp_with_assays$chip_assays)

bp_with_assays_filt = filter(bp_with_assays, chip_assays) # 23 conditions with all epi + rna
table(bp_with_assays_filt$`Sub-group`)

