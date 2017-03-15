cmd_1 = "WiggleTools/wiggletools apply AUC tmp/1a466d277ff3.bed ~/links/CB.HTS.Analysis/CTTV020/data/BLUEPRINT/C004GDH1.ERX207032.H3K27ac.bwa.GRCh38.20150528.bw"
system.time(system(cmd_1, intern=TRUE))
cmd_2 = "WiggleTools/wiggletools apply AUC tmp/1a466d277ff3.bed ~/links/RD-Epigenetics-NetworkData/otar_020/BLUEPRINT/C004GDH1.ERX207032.H3K27ac.bwa.GRCh38.20150528.bw"
system.time(system(cmd_2, intern=TRUE))


