context("AUC per region")

# get regulatory regions
data("roi_reg.RData", package="epiChoose")
data("gene_list_all", package="epiChoose")
data("dat_max_db", package="epiChoose")
data("dat_tss", package="epiChoose")
data_gsk = read_excel(system.file("extdata", "data_gsk.xlsx", package="epiChoose"))

# generate data for test
g = "SAMD11" # pick an example gene
g_loc = gene_list_all[gene_list_all$hgnc_symbol==g]
win = 2e3
s_name = c("A549_BR1_Baseline")
start(g_loc) = start(g_loc) - win
end(g_loc) = end(g_loc) + win

# roi_test = roi_reg[subjectHits(findOverlaps(g_loc, roi_reg))]
# roi_test = as.data.frame(roi_test)[,1:3]
# roi_test = arrange(roi_test, seqnames, start, end)
# write_tsv(roi_test, "inst/extdata/roi_test.bed", col_names=FALSE)
# s_bigwig = filter(data_gsk, Label==s_name, Mark=="H3K27ac") %>% dplyr::select(Bigwig) %>% unlist()
# cmd = paste0("WiggleTools/wiggletools apply AUC inst/extdata/roi_test.bed ", s_bigwig, " > inst/extdata/roi_test.out")

roi_test = read_tsv(system.file("extdata", "roi_test.out", package="epiChoose"), col_names=FALSE)

test_that("make_auc_matrix.R is returning the correct summary gene metrics (max):",
          {
            expect_equal(
              dat_max_gb$H3K27ac$res[
                which(rownames(dat_max_gb$H3K27ac$res)==s_name),
                which(colnames(dat_max_gb$H3K27ac$res)==g)
                ], 
              max(total_data$H3K27ac$res[
                which(rownames(total_data$H3K27ac$res)==s_name),
                subjectHits(findOverlaps(g_loc, roi_reg))
                ])
              )
          }
)

g_loc = gene_list_all[gene_list_all$hgnc_symbol==g]
start(g_loc) = g_loc$transcription_start_site - win
end(g_loc) = g_loc$transcription_start_site + win

test_that("make_auc_matrix.R is returning the correct summary gene metrics (tss):",
          {
            expect_equal(
              dat_tss$H3K27ac$res[
                which(rownames(dat_tss$H3K27ac$res)==s_name),
                which(colnames(dat_tss$H3K27ac$res)==g)
                ], 
              max(total_data$H3K27ac$res[
                which(rownames(total_data$H3K27ac$res)==s_name),
                subjectHits(findOverlaps(g_loc, roi_reg))
                ])
            )
          }
)


