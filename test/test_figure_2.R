# comparison of monocytes/macrophages with u937/thp-1

# 1. establish the amount of paired samples in blueprint
# 2. groups: monocytes, macrophages, u937, thp-1
# 3. how to define dynamic regions?
# 4. figure 4d for dynamic regions i.e. regions that differ between monocytes and macrophages and what they look like in u937/thp-1
# 5. tf inference from deregulated sites (look at luz's review for methods)


# PAIRED SAMPLES IN BLUEPRINT ---------------------------------------------

# get the necessary files from blueprint
bp <- read_tsv("inst/extdata/blueprint_files.tsv")
bp[is.na(bp)] = "" # removes nas that affect dplyr selections

# which donors have monocyte and macrophage samples?
bp_tbl = bp %>% group_by(Donor) %>% summarise("Both"=(any(grepl("monocyte",`Sub-group`,ignore.case=TRUE))) & (any(grepl("macrophage",`Sub-group`,ignore.case=TRUE))))
bp_filt = filter(bp, Donor %in% unlist(bp_tbl[bp_tbl$Both==TRUE, 'Donor']), Format=="bigWig")
bp_filt %>% group_by(Donor) %>% summarise(Marks=all(c("H3K4me3","H3K27ac") %in% Experiment))


# DYNAMIC H3K27AC REGIONS -------------------------------------------------

load("data/total_data.RData")
load("data/dat_all.RData")

dat = dat_all$tss$H3K27ac$res[grep("monocyte|macrophage", rownames(dat_all$tss$H3K27ac$res), ignore.case=TRUE),]
rownames(dat) = str_extract(rownames(dat), "[Mm]onocyte|[Mm]acrophage")
rownames(dat) = tolower(rownames(dat))
rownames(dat) = str_replace(rownames(dat), "_BR[12]", "")

# remove na columns
dim(dat)
na_ix = which(apply(is.na(dat), 2, all))
dat = dat[,-na_ix] # remove regions with no data
# dat = dat[-which(apply(dat, 1, function(x) all(is.na(x)))),]  # remove samples with no data
dim(dat)

# pick high-variance regions
head(order(apply(dat, 2, var, na.rm=TRUE), decreasing=TRUE))
ggplot(data.frame(Type=names(dat[,order(apply(dat, 2, var, na.rm=TRUE), decreasing=TRUE)[2]]), AUC=dat[,order(apply(dat, 2, var, na.rm=TRUE), decreasing=TRUE)[2]]), aes(Type, AUC)) + geom_boxplot() + theme_thesis(20)
col_ix = head(order(apply(dat,2,var,na.rm=TRUE), decreasing=TRUE), 1000)
dat_filt = dat[,col_ix]
pheatmap(log(dat_filt), show_colnames=FALSE)

# kmeans
dat_filt_t = t(dat_filt)
wss <- (nrow(dat_filt_t)-1) * sum(apply(dat_filt_t,2,var))
for(i in 2:15) {
  wss[i] <- sum(kmeans(dat_filt_t, centers=i)$withinss)
}

plot(1:15, wss, type="b", xlab="Number of Clusters", ylab="Within-group sum-of-squares")
fit <- kmeans(dat_filt_t, 7)
col_clust = hclust(dist(dat_filt_t))

pheatmap(log(dat_filt),
         annotation_col = data.frame(Cluster=factor(fit$cluster)),
         cluster_cols = col_clust,
         show_colnames = FALSE
)

which_clust = 7
to_plot = rbind(
  data.frame(Type="Macrophage", AUC=as.numeric(dat_filt[rownames(dat_filt)=="macrophage",fit$cluster==which_clust])),
  data.frame(Type="Monocyte", AUC=as.numeric(dat_filt[rownames(dat_filt)=="monocyte",fit$cluster==which_clust]))
)
ggplot(to_plot, aes(Type, log(AUC))) + geom_boxplot() + theme_thesis(20)

# add thp-1 and u937
dat_add = dat_all$tss$H3K27ac$res[grep("thp-1|u937", rownames(dat_all$tss$H3K27ac$res), ignore.case=TRUE),]
dim(dat_add)
dat_add = dat_add[,-na_ix]
dim(dat_add)[2] == dim(dat)[2]
dat_add_filt = dat_add[,col_ix]
dim(dat_add_filt)
all(colnames(dat_filt) == colnames(dat_add_filt))
dat_add_filt = dat_add_filt[-which(apply(dat_add_filt, 1, function(x) all(is.na(x)))),]
dim(dat_add_filt)

pheatmap(log(rbind(dat_filt, dat_add_filt)),
         annotation_col = data.frame(Cluster=factor(fit$cluster)),
         cluster_cols = col_clust,
         show_colnames = FALSE
)

# what are the clusters?

write_tsv(data.frame(names(which(fit$cluster==7))), "out_1.txt", col_names=FALSE) # monocyte high
# http://amp.pharm.mssm.edu/Enrichr/enrich?dataset=42vf0

write_tsv(data.frame(names(which(fit$cluster==5))), "out_2.txt", col_names=FALSE) # macrophage high
# http://amp.pharm.mssm.edu/Enrichr/enrich?dataset=42vf4


# APPLY T-TEST ------------------------------------------------------------

p_res_g = rep(NA, dim(dat_filt)[2])
for(i in 1:length(p_res_g)) {
  p_res_g[i] = t.test(dat_filt[which(rownames(dat_filt)=="monocyte"),i], dat[which(rownames(dat_filt)=="macrophage"),i], alternative="g")$p.value
}

check_ix = 1
to_plot = rbind(
  data.frame(Type="Macrophage", AUC=as.numeric(dat_filt[rownames(dat_filt)=="macrophage", order(p_res_g)[check_ix]])),
  data.frame(Type="Monocyte", AUC=as.numeric(dat_filt[rownames(dat_filt)=="monocyte", order(p_res_g)[check_ix]]))
)
ggplot(to_plot, aes(Type, log(AUC))) + geom_boxplot() + theme_thesis(20) + ggtitle(colnames(dat_filt)[order(p_res_g)[check_ix]])


