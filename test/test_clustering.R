# start_data

D = newProject(
  prepareData(start_data,
              roi,
              chr=NULL,
              markerNames=names(start_data),
              trace=T)
)

D = addDist(D) # add a simple distance matrix

# too large to cluster, work off pcs?

D = addCellLines(D, list(
  macro=c(8:10,12,13,15,17,18,28:31),
  t.cell=c(11,14,20:21,23,25,27,32)
)) # cell line groups