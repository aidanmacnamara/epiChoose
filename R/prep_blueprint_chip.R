#' @title TO ADD
#' @description TO ADD
#' @description This is a new line ...
#' @details What's this?
#' @return TO ADD
#' @export

prep_blueprint_chip <- function(blueprint_data, assays=c("H3K27ac","H3K4me3","H3K27me3"), out_file, root, rna_annot) {
  
  require(tidyverse)
  require(stringr)
  
  bp <- read_tsv(blueprint_data)
  bp[is.na(bp)] = "" # removes nas that affect dplyr selections
  
  # what are the variables i want to group by
  group_vars = lapply(c("Group","Sub-group","Name","Sex","Cell type","Tissue","Cell line","Disease","Donor"), as.symbol)
  
  chip_with_assays = bp %>% group_by_(.dots=group_vars) %>% summarise (
    Data = all(assays %in% Experiment)
  )
  
  # filter those that have all assays available
  chip_with_assays_filtered = filter(chip_with_assays, Data)
  
  # check if rna data is available
  bp_rna_labels = paste(rna_annot$`Comment[donor ID]`, rna_annot$`Characteristics[cell type]`, sep="_")
  chip_with_assays_filtered$has_rna = paste(chip_with_assays_filtered$Donor, chip_with_assays_filtered$`Cell type`, sep="_") %in% bp_rna_labels
  chip_with_assays_filtered = filter(chip_with_assays_filtered, has_rna)
  
  # add back in novakovic data
  bp_novakovic = filter(bp, grepl("SANQUIN", Name), Format=="bigWig", Experiment %in% assays)
  bp_novakovic_summ = bp_novakovic %>% group_by_(.dots=group_vars) %>% summarise (
    Data = all(assays %in% Experiment)
  )
  
  # add back saeed data
  # saeed donors
  saeed_donors = c(
    "N00031318490130",
    "N00031319896021",
    "N00031401639721",
    "N00031406634921",
    "N00031406635321"
  )
  bp_saeed = filter(bp, (Donor %in% saeed_donors), Format=="bigWig", Experiment %in% assays)
  bp_saeed_summ = bp_saeed %>% group_by_(.dots=group_vars) %>% summarise (
    Data = all(assays %in% Experiment)
  )
  
  chip_with_assays_filtered = full_join(chip_with_assays_filtered, bp_novakovic_summ)
  chip_with_assays_filtered = full_join(chip_with_assays_filtered, bp_saeed_summ)
  
  out_dat = data.frame() # initialise output data frame
  
  # pull out file names etc.
  for(i in 1:dim(chip_with_assays_filtered)[1]) {
    to_add = filter(bp,
                    `Sub-group` == chip_with_assays_filtered$`Sub-group`[i],
                    Name == chip_with_assays_filtered$Name[i],
                    Sex == chip_with_assays_filtered$Sex[i],
                    Group == chip_with_assays_filtered$Group[i],
                    `Cell type` == chip_with_assays_filtered$`Cell type`[i],
                    Tissue == chip_with_assays_filtered$Tissue[i],
                    `Cell line` == chip_with_assays_filtered$`Cell line`[i],
                    Disease == chip_with_assays_filtered$Disease[i],
                    Donor == chip_with_assays_filtered$Donor[i],
                    Experiment %in% assays,
                    Format=="bigWig"
    )
    
    to_add$N = dim(to_add)[1]
    if(to_add$`Cell type`[1]=="") {
      to_add$Label = paste(to_add$Donor[1], to_add$`Cell line`[1], sep="_")
    } else {
      to_add$Label = paste(to_add$Donor[1], to_add$`Cell type`[1], sep="_")
    }
    out_dat = rbind(out_dat, to_add)
  }
  
  # add stimulus conditions to novakovic and saeed data
  out_dat$Label[grepl("(SANQUIN|N00031)", out_dat$Name)] =
    paste(
      out_dat$Donor[grepl("(SANQUIN|N00031)", out_dat$Name)],
      out_dat$`Sub-group`[grepl("(SANQUIN|N00031)", out_dat$Name)], sep="_"
    )
  
  # what chip files are missing?
  files = list.files(root)
  new_dat = str_extract(out_dat$URL, "[[:alnum:]\\.\\_\\-]+$")
  out_dat$Available = new_dat %in% files
  table(out_dat$Available)
  
  out_dat$Bigwig = paste0(root, new_dat)
  out_dat$Bigwig[!out_dat$Available] = NA
  names(out_dat)[which(names(out_dat)=="Experiment")] = "Mark"
  
  write_csv(out_dat, out_file)
  return(list(files=out_dat, samples=chip_with_assays_filtered))
  
}

