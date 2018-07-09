
load("data/dat_all.RData")


# DATA PREP ---------------------------------------------------------------

# list of matrices: features as rows, samples as columns

mofa_input = vector("list", length(dat_all$max))
names(mofa_input) = names(dat_all$max)
mofa_input$H3K27ac = t(as.matrix(dat_all$max$H3K27ac$res))
mofa_input$H3K4me3 = t(as.matrix(dat_all$max$H3K4me3$res))
mofa_input$H3K27me3 = t(as.matrix(dat_all$max$H3K27me3$res))
mofa_input$ATAC = t(as.matrix(dat_all$closest$ATAC$res))
mofa_input$CTCF = t(as.matrix(dat_all$closest$CTCF$res))
mofa_input$RNA = t(as.matrix(dat_all$max$RNA$res))

mofa = createMOFAobject(mofa_input)
mofa

plotTilesData(mofa)

mofa_train = prepareMOFA(mofa)
mofa_model = runMOFA(mofa_train, outfile=tempfile())

