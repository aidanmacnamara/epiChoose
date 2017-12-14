
require(tidyverse)
require(GenomicRanges)
require(biomaRt)
require(devtools)
require(SummarizedExperiment)
require(ggrepel)
require(VennDiagram)
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

roi_reg = makeGRangesFromDataFrame(roi_reg, keep.extra.columns=TRUE, start.field="chromosome_start", end.field="chromosome_end", seqnames.field="chromosome_name")
roi_reg = sort(roi_reg)
save(roi_reg, file="data/roi_reg.RData")

mart = useMart("ENSEMBL_MART_FUNCGEN", dataset="hsapiens_external_feature")
roi_reg_other = getBM(attributes=c("chromosome_name","chromosome_start","chromosome_end","feature_type","feature_type_class"), filters=list(chromosome_name=c(as.character(1:22), "X", "Y")), mart=mart)
roi_reg_other$chromosome_name = paste("chr", roi_reg_other$chromosome_name, sep="")

roi_reg_other = makeGRangesFromDataFrame(roi_reg_other, keep.extra.columns=TRUE, start.field="chromosome_start", end.field="chromosome_end", seqnames.field="chromosome_name")
roi_reg_other = sort(roi_reg_other)
save(roi_reg_other, file="data/roi_reg_other.RData")


# ROI ---------------------------------------------------------------------

# look across all regulatory regions
load("data/roi_reg.RData")
load("data/gene_list_all.RData")
marks = c("H3K27ac","H3K4me3","H3K27me3","ATAC","CTCF")


# BLUEPRINT DATA ----------------------------------------------------------

blueprint_input = "inst/extdata/blueprint_parsed.csv"

require(BiocParallel)
# blueprint_chip = bplapply(seq(along=marks), function(x) make_auc_matrix(blueprint_input, roi, marks[x], "tmp/", quantile_norm=FALSE), BPPARAM=MulticoreParam(workers=3))
blueprint_chip = vector("list", length=length(marks))
for(i in 1:length(blueprint_chip)) {
  blueprint_chip[[i]] = make_auc_matrix(blueprint_input, roi_reg, marks[i], "tmp/", quantile_norm=TRUE)
}

# make sure rows are in the same order
blueprint_chip_filtered = prep_across_datatypes(blueprint_chip)

# take out samples with missing data for 1 or more data type
na_df = data.frame(lapply(blueprint_chip_filtered, function(x) apply(x$res, 1, function(y) !all(is.na(y)))))
names(na_df) = NULL
na_ix = which(apply(na_df, 1, all))

for(i in 1:length(blueprint_chip_filtered)) {
  blueprint_chip_filtered[[i]]$res = blueprint_chip_filtered[[i]]$res[na_ix,]
  blueprint_chip_filtered[[i]]$annot = blueprint_chip_filtered[[i]]$annot[na_ix,]
}


# GSK DATA ----------------------------------------------------------------

gsk_input = "inst/extdata/data_gsk.csv"
gsk_chip = bplapply(seq(along=marks), function(x) make_auc_matrix(gsk_input, roi, marks[x], "tmp/", quantile_norm=TRUE), BPPARAM=MulticoreParam(workers=5))
gsk_chip_filtered = prep_across_datatypes(gsk_chip)


# ENCODE DATA -------------------------------------------------------------

encode_input = "inst/extdata/data_encode.csv"
encode_chip = bplapply(seq(along=marks), function(x) make_auc_matrix(encode_input, roi, marks[x], "tmp/", quantile_norm=FALSE), BPPARAM=MulticoreParam(workers=3))
encode_chip_filtered = prep_across_datatypes(encode_chip)


# COMBINE -----------------------------------------------------------------

mask_data = vector("list", 5)
for(i in 1:length(mask_data)) {
  mask_data[[i]]$res = rbind(blueprint_chip_filtered[[i]]$res, gsk_chip_filtered[[i]]$res, encode_chip_filtered[[i]]$res)
  # renormalize as we are multiple sources
  mask_data[[i]]$res = quantile_norm(mask_data[[i]]$res)
  mask_data[[i]]$annot = bind_rows(blueprint_chip_filtered[[i]]$annot, gsk_chip_filtered[[i]]$annot, encode_chip_filtered[[i]]$annot)
}

names(mask_data) = marks


# GET LABELS --------------------------------------------------------------

single_labels = rownames(mask_data[[1]]$res)
group_labels = c(rep("BLUEPRINT",178), rep("GSK",43), rep("ENCODE",31))

# ADD RNA -----------------------------------------------------------------

u937 = list.files("z:/links/bix-analysis-stv/2016/CTTV/U937/data/Outputs/star/bam_files/genecounts/")
thp1 = list.files("z:/links/bix-analysis-stv/2016/CTTV/THP1/data/genecounts/")

u937_df = lapply(as.list(u937), function(x) read_tsv(paste0("z:/links/bix-analysis-stv/2016/CTTV/U937/data/Outputs/star/bam_files/genecounts/", x), col_names=FALSE))
thp1_df = lapply(as.list(thp1), function(x) read_tsv(paste0("z:/links/bix-analysis-stv/2016/CTTV/THP1/data/genecounts/", x), col_names=FALSE))

# merge into sample per column
my_names = u937_df[[1]]$X1
p2_rna = data.frame(do.call("cbind", lapply(u937_df, function(x) x$X2)))
p2_rna = cbind(NA, my_names, p2_rna)
p2_rna = cbind(p2_rna, data.frame(do.call("cbind", lapply(thp1_df, function(x) x$X2[match(my_names, x$X1)]))))

# get names
sample_info = read_csv("z:/links/bix-analysis-stv/2016/CTTV/U937/data/sampleInfo.csv", col_names=FALSE)

all_names = sample_info$X1[c(match(str_extract(u937, "^[[:alnum:]]+"), sample_info$X32), match(str_extract(thp1, "^[[:alnum:]]+"), sample_info$X32))]

names(p2_rna) = c("Gene ID", "Gene Name", all_names)

# change names so they match with epigenetic labels
names(p2_rna) = str_replace(names(p2_rna), "THP1", "THP-1")
names(p2_rna) = str_replace(names(p2_rna), "CTR", "Baseline")
names(p2_rna) = str_replace(names(p2_rna), "_RNA", "")
names(p2_rna) = str_replace(names(p2_rna), "Baseline\\+LPS", "LPS")

p2_rna = tbl_df(p2_rna)

p1_rna = read_tsv("z:/sandbox/epiChoose/data/rna/E-MTAB-4101-query-results.fpkms.tsv", skip=4)
p3_rna = read_tsv("z:/sandbox/epiChoose/data/rna/E-MTAB-4729-query-results.fpkms.tsv", skip=4)
# bp_rna = read_tsv("z:/sandbox/epiChoose/data/rna/E-MTAB-3827-query-results.fpkms.tsv", skip=4)

rna = tbl_df(merge(merge(p1_rna, p3_rna, all=TRUE), p2_rna, all=TRUE, by="Gene Name"))
names(rna)[2] = "Gene ID"
rna = rna[,-which(names(rna)=="Gene ID.y")]
names(rna)[3:11] = paste0(names(rna)[3:11], "_BR1_Baseline")
names(rna) = str_replace(names(rna), "-", "")
names(rna) = str_replace(names(rna), "P1", "P-1")

# rna_bp = prep_blueprint_rna() # leave this for the moment

rna_add = prep_rna(rna, gene_list_all$hgnc_symbol, single_labels, names(rna)[-c(1,2)])


# EDIT NAMES --------------------------------------------------------------

mask_data[[1]]$annot$Project[which(group_labels=="BLUEPRINT")] = "BLUEPRINT"
mask_data[[1]]$annot$Project[which(group_labels=="ENCODE")] = "ENCODE"

mask_data[[6]] = rna_add
names(mask_data)[6] = "RNA"


# REG TO GENE COLLAPSES ---------------------------------------------------

dat_max = mask_data # max across gene body

for(i in 1:length(dat_max[1:5])) {
  print(paste("Processing data type", names(dat_max)[i]))
  dat_max[[i]]$res = convert_reg_matrix(dat_max[[i]]$res, roi, gene_list_all, reg_window=2e3, summ_method="max")
}

