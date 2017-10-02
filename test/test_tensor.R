load("r_data/tss_data.RData")
load("r_data/column_annotation/tss_list.RData")
load("r_data/column_annotation/roi_ensembl_multicell.RData")
load("r_data/column_annotation/gene_list_all.RData")
load_all("../CTTV-020-11/CMEA")

# choose max across region
names(tss_data)
dat = tss_data[[5]]
data_shape(dat)
summary(dat)

single_labels = rownames(dat[[1]]$res)
group_labels = c(rep("GSK",19), rep("ENCODE",2))

# slice matrices if necessary
sample_ix = c(seq(from=1,to=15,by=2),18)
sample_ix = c(1,3,15)
for(i in 1:length(dat)) {
  dat[[i]]$res = dat[[i]]$res[sample_ix,]
  dat[[i]]$annot = dat[[i]]$annot[sample_ix,]
}

dat = dat[-4] # remove atac (no data yet)

# a = melt(ldply(dat, function(z) {apply(z$res, 1, function(w) {ifelse(all(is.na(w)), 0, 1)})}))
# ggplot(a, aes(.id,variable, colour=as.factor(value))) + theme_thesis() +theme(axis.text.x=element_text(angle=45, hjust=1)) + geom_point(size=5) + labs(x="", y="", colour="present")

data_shape(dat)

# sub na for 0
for(i in 1:length(dat)) {
  dat[[i]]$res[is.na(dat[[i]]$res)] = 0
}

d = newProject(
  prepareData(
    dat,
    gene_list_all,
    chr=NULL,
    markerNames=names(dat),
    trace=F
  )
)

# set.seed(1)
# d_1 = subsetData(d) # withtout options, no data thinning is performed
# dim(d_1)

m = as.tensor(foldSubsetData(d$data))
cp_1 = cp(m, 5) # slider in app?

u_1 = cp_1$U[[1]]
rownames(u_1)=dimnames(m@data)[[1]]

u_3=cp_1$U[[3]]

marker_names = dimnames(m@data)[[3]]

col_table = list(
  "A549"="A549",
  "BEAS2B"="BEAS2B",
  "NHBE"="NHBE"
)

pcs=c(1,2)
plot_grid(
  simpleScoresPlot(u_1, pcs=pcs, colours=col_table, axes=F),
  cpMarkerProfile(u_3, pcs, marker_names),
  nrow=2, ncol=1,
  rel_heights=c(3,1)
)


