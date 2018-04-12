
require(tidyverse)
require(GenomicRanges)
require(biomaRt)
require(devtools)
require(SummarizedExperiment)
require(ggrepel)
require(VennDiagram)
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
gene_list_all = makeGRangesFromDataFrame(gene_list_all, keep.extra.columns=TRUE, start.field="start_position", end.field="end_position")
newNames = paste("chr", levels(seqnames(gene_list_all)), sep="")
names(newNames) = levels(seqnames(gene_list_all))
gene_list_all = renameSeqlevels(gene_list_all, newNames)
gene_list_all = gene_list_all[-1] # no hgnc symbol
gene_list_all = sort(gene_list_all)
save(gene_list_all, file="data/gene_list_all.RData")

# regulatory build
mart = useMart("ENSEMBL_MART_FUNCGEN", dataset="hsapiens_regulatory_feature")
roi_reg = getBM(attributes=c("chromosome_name","chromosome_start","chromosome_end","feature_type_name"), filters=list(chromosome_name=c(as.character(1:22), "X", "Y")), mart=mart)
roi_reg$chromosome_name = paste("chr", roi_reg$chromosome_name, sep="")
roi_reg = arrange(roi_reg, chromosome_name, chromosome_start, chromosome_end)

roi_reg = makeGRangesFromDataFrame(roi_reg, keep.extra.columns=TRUE, start.field="chromosome_start", end.field="chromosome_end", seqnames.field="chromosome_name")
roi_reg_annot = annotatePeak(roi_reg, TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb="org.Hs.eg.db")
roi_reg = roi_reg_annot@anno
save(roi_reg, file="data/roi_reg.RData")

mart = useMart("ENSEMBL_MART_FUNCGEN", dataset="hsapiens_external_feature")
roi_reg_other = getBM(attributes=c("chromosome_name","chromosome_start","chromosome_end","feature_type","feature_type_class"), filters=list(chromosome_name=c(as.character(1:22), "X", "Y")), mart=mart)
roi_reg_other$chromosome_name = paste("chr", roi_reg_other$chromosome_name, sep="")
roi_reg_other = arrange(roi_reg_other, chromosome_name, chromosome_start, chromosome_end)

roi_reg_other = makeGRangesFromDataFrame(roi_reg_other, keep.extra.columns=TRUE, start.field="chromosome_start", end.field="chromosome_end", seqnames.field="chromosome_name")
save(roi_reg_other, file="data/roi_reg_other.RData")


# ROI ---------------------------------------------------------------------

# look across all regulatory regions
load("data/roi_reg.RData")
load("data/gene_list_all.RData")
marks = c("H3K27ac","H3K4me3","H3K27me3","ATAC","CTCF")


# BLUEPRINT DATA ----------------------------------------------------------

# this returns the metadata for samples that have all of H3K27ac, H3K4me3, H3K27me3 and rna available
rna_annot = read_tsv("inst/extdata/rna/E-MTAB-3827.sdrf.txt") # get rna annotation to check against
blueprint_parsed = prep_blueprint_chip(blueprint_data="inst/extdata/blueprint_files.tsv", root="z:/links/RD-Epigenetics-NetworkData/otar_020/BLUEPRINT/", out_file="inst/extdata/blueprint_parsed.csv", rna_annot=rna_annot)

blueprint_input = "inst/extdata/blueprint_parsed.csv"

require(BiocParallel)
# blueprint_chip = bplapply(seq(along=marks), function(x) make_auc_matrix(blueprint_input, roi, marks[x], "tmp/", quantile_norm=FALSE), BPPARAM=MulticoreParam(workers=3))
blueprint_chip = vector("list", length=length(marks))
for(i in 1:length(blueprint_chip)) {
  blueprint_chip[[i]] = make_auc_matrix(blueprint_input, roi_reg, marks[i], "tmp/", quantile_norm=FALSE)
}

# make sure rows are in the same order
blueprint_chip_filtered = prep_across_datatypes(blueprint_chip)


# GSK DATA ----------------------------------------------------------------

gsk_input = "inst/extdata/data_gsk.xlsx"

# check for missing data
data_gsk = read_excel(gsk_input)

bw_missing = data_gsk[!file.exists(str_replace(data_gsk$Bigwig, "/GWD/bioinfo/projects", "z:/links")),]
# bw_missing = data_gsk[!file.exists(data_gsk$Bigwig),]
bw_missing$Label

bam_missing = data_gsk[!file.exists(str_replace(data_gsk$Bam, "/GWD/bioinfo/projects", "z:/links")),]
# bam_missing = data_gsk[!file.exists(data_gsk$Bam),]
bam_missing$Label

gsk_chip = bplapply(seq(along=marks), function(x) make_auc_matrix(gsk_input, roi_reg, marks[x], "tmp/", quantile_norm=FALSE), BPPARAM=MulticoreParam(workers=5))
gsk_chip_filtered = prep_across_datatypes(gsk_chip)


# ENCODE DATA -------------------------------------------------------------

encode_input = "inst/extdata/data_encode.csv"
encode_chip = bplapply(seq(along=marks), function(x) make_auc_matrix(encode_input, roi_reg, marks[x], "tmp/", quantile_norm=FALSE), BPPARAM=MulticoreParam(workers=3))
encode_chip_filtered = prep_across_datatypes(encode_chip)


# COMBINE -----------------------------------------------------------------

total_data = vector("list", length(marks)) # number of data types
for(i in 1:length(total_data)) {
  total_data[[i]]$res = rbind(blueprint_chip_filtered[[i]]$res, gsk_chip_filtered[[i]]$res, encode_chip_filtered[[i]]$res)
  # renormalize as we are multiple sources
  total_data[[i]]$res = quantile_norm(total_data[[i]]$res)
  total_data[[i]]$annot = bind_rows(blueprint_chip_filtered[[i]]$annot, gsk_chip_filtered[[i]]$annot, encode_chip_filtered[[i]]$annot)
}

names(total_data) = marks


# GET LABELS --------------------------------------------------------------

single_labels = rownames(total_data[[1]]$res)
group_labels = c(
  rep("BLUEPRINT", dim(blueprint_chip_filtered[[1]]$res)[1]), 
  rep("GSK", dim(gsk_chip_filtered[[1]]$res)[1]), 
  rep("ENCODE", dim(encode_chip_filtered[[1]]$res)[1])
)


# ADD RNA -----------------------------------------------------------------

rna_dat = sapply(grep("fpkm", list.files("inst/extdata/rna/"), value=TRUE), function(x) read_tsv(paste0("inst/extdata/rna/",x)))

rna_annot = sapply(grep("sdrf", list.files("inst/extdata/rna/"), value=TRUE), function(x) read_tsv(paste0("inst/extdata/rna/",x)))

mart_1 = useMart("ensembl", dataset="hsapiens_gene_ensembl")
mapping <- getBM(attributes=c("ensembl_gene_id","hgnc_symbol"), mart=mart_1)
rownames_symbol = mapping$hgnc_symbol[match(rna_dat[[1]]$`Gene ID`, mapping$ensembl_gene_id)]



rna = tbl_df(merge(merge(p1_rna, p3_rna, all=TRUE), p2_rna, all=TRUE, by="Gene Name"))
names(rna)[2] = "Gene ID"
rna = rna[,-which(names(rna)=="Gene ID.y")]
names(rna)[3:11] = paste0(names(rna)[3:11], "_BR1_Baseline")
names(rna) = str_replace(names(rna), "-", "")
names(rna) = str_replace(names(rna), "P1", "P-1")

# rna_bp = prep_blueprint_rna() # leave this for the moment

rna_add = prep_rna(rna, gene_list_all$hgnc_symbol, single_labels, names(rna)[-c(1,2)])


# EDIT NAMES --------------------------------------------------------------

total_data[[1]]$annot$Project[which(group_labels=="BLUEPRINT")] = "BLUEPRINT"
total_data[[1]]$annot$Project[which(group_labels=="ENCODE")] = "ENCODE"

total_data[[6]] = rna_add
names(total_data)[6] = "RNA"


# REG TO GENE COLLAPSES ---------------------------------------------------

dat_max_gb = total_data # max across gene body
for(i in 1:length(dat_max_gb[1:5])) {
  print(paste("Processing data type", names(dat_max_gb)[i]))
  dat_max_gb[[i]]$res = convert_reg_matrix(dat_max_gb[[i]]$res, roi_reg, gene_list_all, reg_window=2e3, summ_method="max")
}

dat_tss = total_data # tss sites only
for(i in 1:length(dat_tss[1:5])) {
  print(paste("Processing data type", names(dat_tss)[i]))
  dat_tss[[i]]$res = convert_reg_matrix(dat_tss[[i]]$res, roi_reg, gene_list_all, reg_window=2e3, summ_method="tss")
}

dat_sum_gb = total_data # normalised sum of aucs across gene body (h3k27me3 relevant)
for(i in 1:length(dat_sum_gb[1:5])) {
  print(paste("Processing data type", names(dat_sum_gb)[i]))
  dat_sum_gb[[i]]$res = convert_reg_matrix(dat_sum_gb[[i]]$res, roi_reg, gene_list_all, reg_window=2e3, summ_method="sum")
}

# dat_max_10_closest = total_data # max from 10 closest peaks (ctcf relevant)

dat_all = list(dat_max_gb, dat_tss, dat_sum_gb, NA)
names(dat_all) = c("max", "tss", "sum", "closest")


