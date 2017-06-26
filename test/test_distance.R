# start off with mask_data

load("r_data/mask_data.RData")
load("r_data/all_data.RData")

# look at the project 1 result again
my_labels = c("F36P_BR1","F36P_BR2","HEL9217_BR1","HEL9217_BR2","K562_BR1","K562_BR2","KU812_BR1","KU812_BR2","MEG01_BR1","MEG01_BR2","UT7_BR1","UT7_BR2","S004BT_CD34-negative, CD41-positive, CD42-positive megakaryocyte cell","S002R5_erythroblast")
my_labels[1:12] = paste0(my_labels[1:12], "_Baseline")

start_data = all_data

group_labels = c(rep("GSK",12), rep("BLUEPRINT",2))
single_labels = my_labels
sample_ix = match(my_labels, rownames(start_data[[2]]$res))

for(i in 1:length(start_data)) {
  start_data[[i]]$res = start_data[[i]]$res[sample_ix,]
  start_data[[i]]$annot = start_data[[i]]$annot[sample_ix,]
}

pca_data = prep_for_plot(start_data, annot_1=group_labels, annot_2=my_labels, marks=names(start_data), plot_type="mds")

png(filename="out.png", height=800, width=3200)
ggplot(pca_data, aes(x=x, y=y, color=annot_1)) + geom_point(size=5, shape=17) + theme_thesis() + geom_text_repel(aes(label=annot_2), fontface="bold", size=5, force=0.5) + facet_wrap(~mark, nrow=1, scales="free")
dev.off()

res = dist_mat(start_data[2:4], comp_ix=list(1:12, 13:14))

