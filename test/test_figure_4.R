
require(tidyverse)
require(cProfile)
require(epiChoose)

data(cProfileDemo)
load("data/roi_tss.RData")


# CELL LINES VS PRIMARY ---------------------------------------------------

# i.e. cell line differentiation (with pma) is very different to primary monocyte - macrophage differentiation

a_1 = compareTests("c.prime", "c.pma.line", p=c(0.05,0.01,0.0001), bundle=results_bundle, scale=2)
png("out.png", height=800, width=800)
a_1$graph + theme_thesis(20) + ylab("Cell Line Stimulation") + xlab("Primary Stimulation") + ggtitle(label="")
dev.off()

a_2 = compareTests("c.prime", "c.pma.line", v1="-log10(pvalue)", bundle=results_bundle, scale=NA)

# none of these are significant
# i.e. there is no common mechanism across thp-1 and u937 that is regulated in the same way upon pma stimulation as in primary monocyte - macrophage

# what behaves in common between the cell lines?

b_1 = compareTests("c.U937.pma.prime", "c.THP1.pma.prime", p=c(0.05,0.01,0.0001), bundle=results_bundle, scale=2)
png("out.png", height=800, width=800)
b_1$graph + theme_thesis(20) + xlab("U937 Stimulation") + ylab("THP-1 Stimulation") + ggtitle(label="")
dev.off()

b_2 = compareTests("c.U937.pma.prime", "c.THP1.pma.prime", v1="-log10(pvalue)", bundle=results_bundle, scale=NA)

dat = tbl_df(b_1$data)
dat$gene = roi_tss$hgnc_symbol
q = qnorm(1e-4/2)

dat = dat %>% filter(x1 > -q, x2 > -q) # upregulated
# dat = dat %>% filter(x1 < q, x2 < q) # downregulated

dat = dat %>% dplyr::mutate(comb = rowMeans(dplyr::select(dat,x1,x2)))
dat = arrange(dat, desc(comb))

my_p = list(
  go_bp = gmtPathways("tmp/c5.bp.v6.2.symbols.gmt"), # go biological processes,
  reactome = gmtPathways("tmp/c2.cp.reactome.v6.2.symbols.gmt") # reactome
)

gsea_in = dat$comb
names(gsea_in) = dat$gene
res_gsea <- lapply(my_p, function(x) fgsea(x, stats=gsea_in, nperm=1000))
res_gsea = lapply(res_gsea, function(x) x %>% as_tibble() %>% dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% arrange(padj))

# look at cell lines separately
# david to add thp1 pma vs. baseline to results_bundle?

c_1 = compareTests("c.prime", "c.pma.line", p=c(0.05,0.01,0.0001), bundle=results_bundle, scale=2)
a_1$graph + theme_thesis(15) + ylab("Cell Line Stimulation") + xlab("Primary Stimulation") + ggtitle(label="")


# DATA --------------------------------------------------------------------

require(plyr)
require(ggplot2)
require(reshape)
require(cProfile)
require(DESeq2)

my_n = 500
n_clusters = 10
data(cProfileDemo)
ls()

source("R/joint_regulation_plot.R")

output_off_diagonals <- function(u) {
  
  sets = subset(u$allocation, x!=y)
  sets$cluster=with(sets, paste0("(",x,",",y,")"))
  res = as.data.frame(sets)
  res$gene_symbol = roi_tss$hgnc_symbol[as.numeric(res$gene)]
  res = tbl_df(res)
  
  return(res)
  
}


# THP-1 PMA VS. VD3 -------------------------------------------------------

n = jointGeneOrdering(bundle=results_bundle, contrasts=c("c.THP1.pma.prime","c.THP1.vd3.prime"), n=my_n, metric=min)

a_pma = filterProfiles(recodeLevels(selectData("c.THP1.pma.prime"),"PMA"), results_bundle$results[["c.THP1.pma.prime"]],n=-n, average=T, trans=log10p1)

a_vd3 = filterProfiles(recodeLevels(selectData("c.THP1.vd3.prime"),"VD3"), results_bundle$results[["c.THP1.vd3.prime"]], n=-n, average=T, trans=log10p1)

a_2 = makeCrossCluster(U=list(pma=a_pma, vd3=a_vd3), nClusters=n_clusters, nstart=10)

png("out.png", height=800, width=1200)
plotCrossClusterGrid(a_2, yscale=NULL) + theme_thesis(15)
dev.off()

my_plot = joint_regulation_plot(a_2, contrasts=list(list(contrast="c.THP1.pma.prime", stimulus="PMA"), list(contrast="c.THP1.vd3.prime", stimulus="VD3")), onlyGeneSets=T, n_sets=5, theme_size=15)
png("out.png", height=500, width=1000)
my_plot$g_12
dev.off()

res = output_off_diagonals(a_2)

write_csv(res, "~/Downloads/thp1.csv")


# U937 PMA VS. VD3  -------------------------------------------------------

n = jointGeneOrdering(bundle=results_bundle, contrasts=c("c.U937.pma.prime","c.U937.vd3.prime"), n=my_n, metric=min)

b_pma = filterProfiles(recodeLevels(selectData("c.U937.pma.prime"),"PMA"), results_bundle$results[["c.U937.pma.prime"]],n=-n, average=T, trans=log10p1)

b_vd3 = filterProfiles(recodeLevels(selectData("c.U937.vd3.prime"),"VD3"), results_bundle$results[["c.U937.vd3.prime"]], n=-n, average=T, trans=log10p1)

b_2 = makeCrossCluster(U=list(pma=b_pma, vd3=b_vd3), nClusters=n_clusters, nstart=10)

png("out.png", height=800, width=1200)
plotCrossClusterGrid(b_2, yscale=NULL) + theme_thesis(15)
dev.off()

my_plot = joint_regulation_plot(b_2, contrasts=list(list(contrast="c.U937.pma.prime", stimulus="PMA"), list(contrast="c.U937.vd3.prime", stimulus="VD3")), onlyGeneSets=T, n_sets=5, theme_size=15)
png("out.png", height=500, width=1000)
my_plot$g_12
dev.off()

res = output_off_diagonals(b_2)

write_csv(res, "~/Downloads/u937.csv")


