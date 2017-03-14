
require(tidyverse)
require(GenomicRanges)
require(biomaRt)
require(devtools)
load_all()


# GENES -------------------------------------------------------------------

# link to ensembl
mart_1 = useMart("ensembl", dataset="hsapiens_gene_ensembl")

gene_list_all = getBM(attributes=c("ensembl_gene_id","hgnc_symbol","chromosome_name","start_position","end_position","strand","transcription_start_site"), mart=mart_1, filters=list(chromosome_name=c(as.character(1:22), "X", "Y"), with_protein_id=TRUE))

# pick out 1 tss per gene
pick_tss <- function(x, p_window=100) {
  
  if(x$strand[1]==1) {
    promoter_end = min(x$transcription_start_site) + p_window
    promoter_start = min(x$transcription_start_site) - p_window
  } else if (x$strand[1]==-1) {
    promoter_start = max(x$transcription_start_site) - p_window
    promoter_end = max(x$transcription_start_site) + p_window
  } else {
    promoter_start = NA
    promoter_end = NA
  }
  
  return(data.frame(cbind(x[1,], promoter_start, promoter_end)))
}

gene_list_all = group_by(gene_list_all, hgnc_symbol) %>% do(pick_tss(., 500))
gene_list_all$strand[gene_list_all$strand==1] = "+"
gene_list_all$strand[gene_list_all$strand==-1] = "-"
gene_list_all = makeGRangesFromDataFrame(gene_list_all, keep.extra.columns=TRUE, start.field="start_position", end.field="end_position")
newNames = paste("chr", levels(seqnames(gene_list_all)), sep="")
names(newNames) = levels(seqnames(gene_list_all))
gene_list_all = renameSeqlevels(gene_list_all, newNames)
gene_list_all = gene_list_all[-1]
gene_list_all = sort(gene_list_all)

# ROI ---------------------------------------------------------------------

my_promoters_gene_all = data.frame(
  seqnames = as.character(seqnames(gene_list_all)),
  start = mcols(gene_list_all)$promoter_start,
  end = mcols(gene_list_all)$promoter_end,
  gene = mcols(gene_list_all)$hgnc_symbol
)
my_promoters_gene_all = dplyr::arrange(my_promoters_gene_all, seqnames, start, end) # make sure list is sorted
roi = makeGRangesFromDataFrame(my_promoters_gene_all, keep.extra.columns=TRUE)


# FIND AUC ----------------------------------------------------------------

input_data = "data/data.csv"
dat_auc = make_auc_matrix(input_data, roi, "H3K27ac", "tmp/")
dat_count = make_count_matrix(input_data, roi, "H3K27ac")


# PLOT --------------------------------------------------------------------

pca_res = plot_pca(dat_auc, "GSK", dat$Label, out_file="out.png")

