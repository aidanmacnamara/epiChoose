
require(tidyverse)
require(GenomicRanges)
require(biomaRt)
require(devtools)
require(SummarizedExperiment)
require(ggrepel)
require(VennDiagram)
load_all()

# idea: use sushi plots to compare regulatory links to relevant genes between primary and cell line ATAC signals

# 1. don't work directly with the chicago data