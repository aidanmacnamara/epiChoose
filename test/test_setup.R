
require(tidyverse)
require(GenomicRanges)
require(biomaRt)
require(devtools)
require(SummarizedExperiment)
require(ggrepel)
require(readxl)
require(TxDb.Hsapiens.UCSC.hg38.knownGene) # for peak to gene
require(ChIPseeker) # for peak to gene
require(org.Hs.eg.db) # for peak to gene
load_all()


# GENERATE ENSEMBL DATA ---------------------------------------------------

mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")

# genes
gene_list_all = getBM(attributes=c("ensembl_gene_id","hgnc_symbol","chromosome_name","start_position","end_position","strand","transcription_start_site"), mart=mart, filters=list(chromosome_name=c(as.character(1:22), "X", "Y"), with_protein_id=TRUE))

# pick out 1 tss per gene (min position)
pick_tss <- function(x) {
  if(x$strand[1]==1) {
    tss_ix = which.min(x$transcription_start_site)
  }
  if (x$strand[1]==-1) {
    tss_ix = which.max(x$transcription_start_site)
  } 
  return(x[tss_ix,])
}

gene_list_all = group_by(gene_list_all, hgnc_symbol) %>% do(pick_tss(.))
gene_list_all$strand[gene_list_all$strand==1] = "+"
gene_list_all$strand[gene_list_all$strand==-1] = "-"
gene_list_all = arrange(gene_list_all, chromosome_name, start_position, end_position)
gene_list_all = makeGRangesFromDataFrame(gene_list_all, keep.extra.columns=TRUE, start.field="start_position", end.field="end_position")
newNames = paste("chr", levels(seqnames(gene_list_all)), sep="")
names(newNames) = levels(seqnames(gene_list_all))
gene_list_all = renameSeqlevels(gene_list_all, newNames)

save(gene_list_all, file="data/gene_list_all.RData") # savepoint

# regulatory build
mart = useMart("ENSEMBL_MART_FUNCGEN", dataset="hsapiens_regulatory_feature")
roi_reg = getBM(attributes=c("chromosome_name","chromosome_start","chromosome_end","feature_type_name"), filters=list(chromosome_name=c(as.character(1:22), "X", "Y")), mart=mart)
roi_reg$chromosome_name = paste("chr", roi_reg$chromosome_name, sep="")
roi_reg = arrange(roi_reg, chromosome_name, chromosome_start, chromosome_end)

roi_reg = makeGRangesFromDataFrame(roi_reg, keep.extra.columns=TRUE, start.field="chromosome_start", end.field="chromosome_end", seqnames.field="chromosome_name")
roi_reg_annot = annotatePeak(roi_reg, TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb="org.Hs.eg.db")
roi_reg = roi_reg_annot@anno
save(roi_reg, file="data/roi_reg.RData") # savepoint

mart = useMart("ENSEMBL_MART_FUNCGEN", dataset="hsapiens_external_feature")
roi_reg_other = getBM(attributes=c("chromosome_name","chromosome_start","chromosome_end","feature_type","feature_type_class"), filters=list(chromosome_name=c(as.character(1:22), "X", "Y")), mart=mart)
roi_reg_other$chromosome_name = paste("chr", roi_reg_other$chromosome_name, sep="")
roi_reg_other = arrange(roi_reg_other, chromosome_name, chromosome_start, chromosome_end)

roi_reg_other = makeGRangesFromDataFrame(roi_reg_other, keep.extra.columns=TRUE, start.field="chromosome_start", end.field="chromosome_end", seqnames.field="chromosome_name")
save(roi_reg_other, file="data/roi_reg_other.RData") # savepoint

roi_tss = gene_list_all
tss_window = 2e3
start(roi_tss) = roi_tss$transcription_start_site - tss_window
end(roi_tss) = roi_tss$transcription_start_site + tss_window
roi_tss_annot = annotatePeak(roi_tss, TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb="org.Hs.eg.db")
roi_tss = roi_tss_annot@anno

# roi_tss not ordered (despite 'sort')
roi_tss = as.data.frame(roi_tss)
roi_tss = arrange(roi_tss, as.character(seqnames), start, end)
roi_tss = makeGRangesFromDataFrame(roi_tss, keep.extra.columns=TRUE)

save(roi_tss, file="data/roi_tss.RData") # savepoint


# ROI ---------------------------------------------------------------------

# look across all regulatory regions
load("data/roi_reg.RData")
load("data/gene_list_all.RData")
marks = c("H3K27ac","H3K4me3","H3K27me3","ATAC","CTCF")


# BLUEPRINT DATA ----------------------------------------------------------

# this returns the metadata for samples that have all of H3K27ac, H3K4me3, H3K27me3 and rna available
rna_annot = read_tsv("inst/extdata/rna/E-MTAB-3827.sdrf.txt") # get rna annotation to check against
data_blueprint = prep_blueprint_chip(blueprint_data="inst/extdata/blueprint_files.tsv", root="/GWD/bioinfo/projects/RD-Epigenetics-NetworkData/otar_020/BLUEPRINT/", out_file="inst/extdata/data_blueprint.csv", rna_annot=rna_annot)

blueprint_input = "inst/extdata/data_blueprint.csv"

require(BiocParallel)

# blueprint_chip = bplapply(seq(along=marks), function(x) make_auc_matrix(blueprint_input, roi, marks[x], "tmp/", quantile_norm=FALSE), BPPARAM=MulticoreParam(workers=3))

blueprint_chip = vector("list", length=length(marks))
for(i in 1:length(blueprint_chip)) {
  blueprint_chip[[i]] = make_auc_matrix(blueprint_input, roi_reg, marks[i], "tmp/", quantile_norm=FALSE)
}
blueprint_chip_filtered = prep_across_datatypes(blueprint_chip) # make sure rows are in the same order

blueprint_tsss = vector("list", length=length(marks))
for(i in 1:length(blueprint_tsss)) {
  blueprint_tsss[[i]] = make_auc_matrix(blueprint_input, roi_tss, marks[i], "tmp/", quantile_norm=FALSE)
}
blueprint_tsss_filtered = prep_across_datatypes(blueprint_tsss)

save(blueprint_chip_filtered, file="data/blueprint_chip_filtered.RData") # savepoint
save(blueprint_tsss_filtered, file="data/blueprint_tsss_filtered.RData") # savepoint


# GSK DATA ----------------------------------------------------------------

gsk_input = "inst/extdata/data_gsk.xlsx"

# check for missing data
data_gsk = read_excel(gsk_input)

# bw_missing = data_gsk[!file.exists(str_replace(data_gsk$Bigwig, "/GWD/bioinfo/projects", "z:/links")),]
bw_missing = data_gsk[!file.exists(data_gsk$Bigwig),]
bw_missing$Label

# bam_missing = data_gsk[!file.exists(str_replace(data_gsk$Bam, "/GWD/bioinfo/projects", "z:/links")),]
bam_missing = data_gsk[!file.exists(data_gsk$Bam),]
bam_missing$Label

# check bw == bam filenames
table(
  str_replace(data_gsk$Bigwig, "^.*/bigwig/(.*)\\.bw$", "\\1") ==
    str_replace(data_gsk$Bam, "^.*/bam/(.*)\\.bam$", "\\1")
)
# 5 files: project 1 ctcf where the bigwig file was not updated

gsk_chip = bplapply(seq(along=marks), function(x) make_auc_matrix(gsk_input, roi_reg, marks[x], "tmp/", quantile_norm=FALSE), BPPARAM=MulticoreParam(workers=5))
gsk_chip_filtered = prep_across_datatypes(gsk_chip)

gsk_tsss = bplapply(seq(along=marks), function(x) make_auc_matrix(gsk_input, roi_tss, marks[x], "tmp/", quantile_norm=FALSE), BPPARAM=MulticoreParam(workers=5))
gsk_tsss_filtered = prep_across_datatypes(gsk_tsss)

save(gsk_chip_filtered, file="data/gsk_chip_filtered.RData") # savepoint
save(gsk_tsss_filtered, file="data/gsk_tsss_filtered.RData") # savepoint


# ENCODE DATA -------------------------------------------------------------

encode_input = "inst/extdata/data_encode.csv"

encode_chip = bplapply(seq(along=marks), function(x) make_auc_matrix(encode_input, roi_reg, marks[x], "tmp/", quantile_norm=FALSE), BPPARAM=MulticoreParam(workers=3))
encode_chip_filtered = prep_across_datatypes(encode_chip)

encode_tsss = bplapply(seq(along=marks), function(x) make_auc_matrix(encode_input, roi_tss, marks[x], "tmp/", quantile_norm=FALSE), BPPARAM=MulticoreParam(workers=3))
encode_tsss_filtered = prep_across_datatypes(encode_tsss)

save(encode_chip_filtered, file="data/encode_chip_filtered.RData") # savepoint
save(encode_tsss_filtered, file="data/encode_tsss_filtered.RData") # savepoint


# DEEP DATA ---------------------------------------------------------------

deep_input = "inst/extdata/data_deep.xlsx"

deep_chip = bplapply(seq(along=marks), function(x) make_auc_matrix(deep_input, roi_reg, marks[x], "tmp/", quantile_norm=FALSE), BPPARAM=MulticoreParam(workers=3))
deep_chip_filtered = prep_across_datatypes(deep_chip)

deep_tsss = bplapply(seq(along=marks), function(x) make_auc_matrix(deep_input, roi_tss, marks[x], "tmp/", quantile_norm=FALSE), BPPARAM=MulticoreParam(workers=3))
deep_tsss_filtered = prep_across_datatypes(deep_tsss)

save(deep_chip_filtered, file="data/deep_chip_filtered.RData") # savepoint
save(deep_tsss_filtered, file="data/deep_tsss_filtered.RData") # savepoint


# COMBINE -----------------------------------------------------------------

total_reg = vector("list", length(marks)) # number of data types
for(i in 1:length(total_reg)) {
  total_reg[[i]]$res = rbind(blueprint_chip_filtered[[i]]$res, gsk_chip_filtered[[i]]$res, encode_chip_filtered[[i]]$res, deep_chip_filtered[[i]]$res)
  # renormalize as we are multiple sources
  total_reg[[i]]$res = quantile_norm(total_reg[[i]]$res)
  if(length(deep_chip_filtered[[i]]$annot$Size)) { # for error with data type conversion
    deep_chip_filtered[[i]]$annot$Size = as.character(deep_chip_filtered[[i]]$annot$Size)
  }
  total_reg[[i]]$annot = bind_rows(blueprint_chip_filtered[[i]]$annot, gsk_chip_filtered[[i]]$annot, encode_chip_filtered[[i]]$annot, deep_chip_filtered[[i]]$annot)
}
names(total_reg) = marks

total_tss = vector("list", length(marks)) # number of data types
for(i in 1:length(total_tss)) {
  total_tss[[i]]$res = rbind(blueprint_tsss_filtered[[i]]$res, gsk_tsss_filtered[[i]]$res, encode_tsss_filtered[[i]]$res, deep_tsss_filtered[[i]]$res)
  # renormalize as we are multiple sources
  total_tss[[i]]$res = quantile_norm(total_tss[[i]]$res)
  if(length(deep_tsss_filtered[[i]]$annot$Size)) { # for error with data type conversion
    deep_tsss_filtered[[i]]$annot$Size = as.character(deep_tsss_filtered[[i]]$annot$Size)
  }
  total_tss[[i]]$annot = bind_rows(blueprint_tsss_filtered[[i]]$annot, gsk_tsss_filtered[[i]]$annot, encode_tsss_filtered[[i]]$annot, deep_tsss_filtered[[i]]$annot)
}
names(total_tss) = marks


# GET LABELS --------------------------------------------------------------

single_labels = rownames(total_reg[[1]]$res)
group_labels = c(
  rep("BLUEPRINT", dim(blueprint_chip_filtered[[1]]$res)[1]),
  rep("GSK", dim(gsk_chip_filtered[[1]]$res)[1]),
  rep("ENCODE", dim(encode_chip_filtered[[1]]$res)[1]),
  rep("DEEP", dim(deep_chip_filtered[[1]]$res)[1])
  # rep("BLUEPRINT", 82),
  # rep("GSK", 82),
  # rep("ENCODE", 31),
  # rep("DEEP", 5)
)

# total_reg_edit = total_tss
# gsk_ix = which(group_labels=="DEEP" | total_reg[[1]]$annot$Project %in% c(4))
# for(i in 1:length(total_reg_edit)) {
#   total_reg_edit[[i]]$res = total_reg_edit[[i]]$res[gsk_ix,]  
# }

# plot_data = prep_for_plot(total_reg_edit, annot_1=group_labels[gsk_ix], annot_2=single_labels[gsk_ix], marks=marks, plot_type="mds")
# pdf(file="out_4_mds.pdf", height=24, width=96)
# ggplot(plot_data, aes(x=x, y=y, color=annot_1)) + geom_point(size=8, shape=17) + theme_thesis(50) + geom_text_repel(aes(label=annot_2), fontface="bold", size=5, force=0.5) + facet_wrap(~mark, nrow=1)
# dev.off()


# ADD RNA -----------------------------------------------------------------

rna_dat = sapply(grep("fpkm", list.files("inst/extdata/rna/"), value=TRUE), function(x) read_tsv(paste0("inst/extdata/rna/",x)))

rna_annot = sapply(grep("sdrf", list.files("inst/extdata/rna/"), value=TRUE), function(x) read_tsv(paste0("inst/extdata/rna/",x)))

# merge data
rna_dat_merge = rna_dat %>% Reduce(function(x, y) full_join(x, y, by="Gene ID"), .)
rna_annot_merge = bind_rows(rna_annot)


# LABEL DATA - NEEDS MANUAL EDITING ---------------------------------------

# match data to annotation by run id
sample_match = match(names(rna_dat_merge)[-1], rna_annot_merge$`Comment[ENA_RUN]`)

# create blueprint labels first
sample_labels = paste(rna_annot_merge$`Comment[donor ID]`[sample_match], rna_annot_merge$`Characteristics[cell type]`[sample_match], sep="_")

# now gsk labels
gsk_label_mapping = read_excel("inst/extdata/rna/rna_chip_mapping.xlsx")
sample_labels[grep("^NA", sample_labels)] = rna_annot_merge$`Source Name`[sample_match[grep("^NA", sample_labels)]]
sample_labels[sample_labels %in% gsk_label_mapping$rna_label] = gsk_label_mapping$chip_label[match(sample_labels[sample_labels %in% gsk_label_mapping$rna_label], gsk_label_mapping$rna_label)]
names(rna_dat_merge)[-1] = sample_labels
names(rna_dat_merge)[-1][is.na(names(rna_dat_merge)[-1])] = paste("Unknown", 1:length(names(rna_dat_merge)[-1][is.na(names(rna_dat_merge)[-1])]), sep="_")
rna_dat_merge = rna_dat_merge[,!duplicated(names(rna_dat_merge))] # IS THIS CORRECT? some biological replicates still appearing twice in the data, removed here ...

# / LABEL DATA / ----------------------------------------------------------


mart_1 = useMart("ensembl", dataset="hsapiens_gene_ensembl")
mapping <- getBM(attributes=c("ensembl_gene_id","hgnc_symbol"), mart=mart_1)
rownames_symbol = mapping$hgnc_symbol[match(rna_dat[[1]]$`Gene ID`, mapping$ensembl_gene_id)]
rna_dat_merge = tbl_df(cbind(rownames_symbol, rna_dat_merge))
names(rna_dat_merge)[1] = "Gene Name"

rna_add = prep_rna(fpkm_table=rna_dat_merge, gene_list=gene_list_all$hgnc_symbol, chip_labels=single_labels, rna_labels=names(rna_dat_merge)[-c(1,2)], quantile_norm=TRUE)


# EDIT NAMES --------------------------------------------------------------

total_reg[[1]]$annot$Project[which(group_labels=="BLUEPRINT")] = "BLUEPRINT"
total_reg[[1]]$annot$Project[which(group_labels=="ENCODE")] = "ENCODE"
total_reg[[1]]$annot$Project[which(group_labels=="DEEP")] = "DEEP"

total_reg[[6]] = rna_add
names(total_reg)[6] = "RNA"

save(total_reg, file="data/total_reg.RData") # savepoint

total_tss[[1]]$annot$Project[which(group_labels=="BLUEPRINT")] = "BLUEPRINT"
total_tss[[1]]$annot$Project[which(group_labels=="ENCODE")] = "ENCODE"
total_tss[[1]]$annot$Project[which(group_labels=="DEEP")] = "DEEP"

save(total_tss, file="data/total_tss.RData") # savepoint


# REG TO GENE COLLAPSES ---------------------------------------------------

dat_max_gb = total_reg # max across gene body
loc_max_gb = vector("list", length(total_reg))
names(loc_max_gb) = names(total_reg)

for(i in 1:length(dat_max_gb[1:5])) {
  print(paste("Processing data type", names(dat_max_gb)[i]))
  my_data = convert_reg_matrix(dat_max_gb[[i]]$res, roi_reg, gene_list_all, reg_window=2e3, summ_method="max")
  dat_max_gb[[i]]$res = my_data$data
  loc_max_gb[[i]] = my_data$locations
}

dat_tss = total_reg # tss sites only
# loc_max_gb = vector("list", length(total_reg))
# names(loc_max_gb) = names(total_reg)
for(i in 1:length(dat_tss[1:5])) {
  print(paste("Processing data type", names(dat_tss)[i]))
  dat_tss[[i]]$res = total_tss[[i]]$res
}

dat_sum_gb = total_reg # normalised sum of aucs across gene body (h3k27me3 relevant)
loc_sum_gb = vector("list", length(total_reg))
names(loc_sum_gb) = names(total_reg)

for(i in 1:length(dat_sum_gb[1:5])) {
  print(paste("Processing data type", names(dat_sum_gb)[i]))
  my_data = convert_reg_matrix(dat_sum_gb[[i]]$res, roi_reg, gene_list_all, reg_window=2e3, summ_method="sum")
  dat_sum_gb[[i]]$res = my_data$data
  loc_sum_gb[[i]] = my_data$locations 
}

dat_max_10 = total_reg # max from 10 closest peaks (ctcf relevant)
loc_max_10 = vector("list", length(total_reg))
names(loc_max_10) = names(total_reg)

for(i in 1:length(dat_max_10[1:5])) {
  print(paste("Processing data type", names(dat_max_10)[i]))
  my_data = convert_reg_matrix(dat_max_10[[i]]$res, roi_reg, gene_list_all, reg_window=2e3, summ_method="closest")
  dat_max_10[[i]]$res = my_data$data
  loc_max_10[[i]] = my_data$locations
}

dat_all = list(dat_max_gb, dat_tss, dat_sum_gb, dat_max_10)
names(dat_all) = c("max", "tss", "sum", "closest")

save(dat_all, file="data/dat_all.RData") # savepoint


# SAVE DATA TO EPIVIEW ----------------------------------------------------

save(dat_all, file="../epiView/data/dat_all.RData") # savepoint
save(gene_list_all, file="../epiView/data/gene_list_all.RData") # savepoint
roi_reg_df = as.data.frame(roi_reg)[,1:3]
save(roi_reg_df, file="../epiView/data/roi_reg_df.RData") # savepoint


