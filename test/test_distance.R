# start off with mask_data

load("r_data/mask_data.RData")
load("r_data/all_data.RData")

# look at the project 1 result again

start_data = mask_data

sample_ix = 114:137 # c(grep("monocyte|macrophage", rownames(mask_data[[1]]$res)), 114:137)
sample_ix = c(77:88,46,43)

single_labels = rownames(start_data[[1]]$res)[sample_ix]
# single_labels = str_extract(single_labels, "[[:alnum:]\\+]+$")
# single_labels = rownames(mask_data[[1]]$res)

# single_labels[1:20] = str_replace(single_labels[1:20], "[[:alnum:]]+_", "")
# group_labels = c(rep("THP-1",12), rep("U937", 12)) # c(rep("Primary",20), rep("Cell Line",24))
group_labels = c(rep("GSK",12), rep("Blueprint",2))
# group_labels = c(rep("BLUEPRINT",76), rep("GSK",12), rep("ENCODE",18), rep("GSK",31))

for(i in 1:length(start_data)) {
  # start_data[[i]]$res = log(start_data[[i]]$res[sample_ix,])
  start_data[[i]]$res = start_data[[i]]$res[sample_ix,]
  start_data[[i]]$annot = start_data[[i]]$annot[sample_ix,]
}

pca_data = prep_for_plot(start_data, annot_1=group_labels, annot_2=single_labels, marks=names(start_data), plot_type="mds")

png(filename="out.png", height=1200, width=3600)
ggplot(pca_data, aes(x=x, y=y, color=annot_1)) + geom_point(size=5, shape=17) + theme_thesis() + geom_text_repel(aes(label=annot_2), fontface="bold", size=5, force=0.5) + facet_wrap(~mark, nrow=1, scales="free")
dev.off()

res = dist_mat(start_data, comp_ix=list(1:12, 13:14))

