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
dat = dat[,!apply(is.na(dat), 2, all)] # remove regions with no data
# dat = dat[-which(apply(dat, 1, function(x) all(is.na(x)))),]  # remove samples with no data
dim(dat)

# pick high-variance regions
head(order(apply(dat, 2, var, na.rm=TRUE), decreasing=TRUE))
ggplot(data.frame(Type=names(dat[,order(apply(dat, 2, var, na.rm=TRUE), decreasing=TRUE)[2]]), AUC=dat[,order(apply(dat, 2, var, na.rm=TRUE), decreasing=TRUE)[2]]), aes(Type, AUC)) + geom_boxplot() + theme_thesis(20)
dat_filt = dat[,head(order(apply(dat,2,var,na.rm=TRUE), decreasing=TRUE), 1000)]
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

which_clust = 6
to_plot = rbind(
  data.frame(Type="Macrophage", AUC=as.numeric(dat_filt[rownames(dat_filt)=="macrophage",fit$cluster==which_clust])),
  data.frame(Type="Monocyte", AUC=as.numeric(dat_filt[rownames(dat_filt)=="monocyte",fit$cluster==which_clust]))
)
ggplot(to_plot, aes(Type, log(AUC))) + geom_boxplot() + theme_thesis(20)




