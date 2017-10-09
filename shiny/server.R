require(epiChoose)
require(shiny)
require(tidyverse)
require(ggrepel)
require(reshape2)
require(biomaRt)
require(rtracklayer)
require(Sushi)

load("gene_list_all.RData")
load("t_list.RData")
load("dat.RData")
load("msig_go_bp.RData")

shinyServer(function(input, output) {
  
  my_choices <- reactive(
    {
      
      my_choices = list()
      
      gene_choice <- renderPrint({input$gene_choice})()
      my_choices$gene_choice = unlist(regmatches(gene_choice, gregexpr("[A-Z][[:alnum:]]+", gene_choice)))
      
      go_choice <- renderPrint({input$go_choice})()
      my_choices$go_choice = str_replace(go_choice, "^.*\"([[:alnum:]\\s]+)\"", "\\1")
      # print(my_choices$go_choice)
      
      cell_target_choice <- renderPrint({input$cell_target_choice})()
      my_choices$cell_target_choice = unlist(regmatches(cell_target_choice, gregexpr("[[:alnum:]_]+", cell_target_choice)))[-1] # remove na added
      
      cell_candidate_choice <- renderPrint({input$cell_candidate_choice})()
      my_choices$cell_candidate_choice = unlist(regmatches(cell_candidate_choice, gregexpr("[[:alnum:]_]+", cell_candidate_choice)))[-1] # remove na added
      
      return(my_choices)
    }
  )
  
  output$dist_plot_1 <- renderPlot({
    
    if(
      (my_choices()$gene_choice[1]!="NULL" | my_choices()$go_choice[1]!="NULL") &
      my_choices()$cell_target_choice[1]!="NULL" &
      my_choices()$cell_candidate_choice[1]!="NULL"
    ) {
      
      if(my_choices()$go_choice[1]!="NULL") {
        genes = msig_go_bp[[which(names(msig_go_bp)==my_choices()$go_choice)]]
      } else {
        genes = my_choices()$gene_choice
      }
      # genes = go_genes[[2]]
      
      c_cells = my_choices()$cell_candidate_choice
      # c_cells = c("A549_BR1_Baseline","A549_BR2_Baseline","BEAS2B_BR1_Baseline","BEAS2B_BR2_Baseline")
      
      t_cells = my_choices()$cell_target_choice
      # t_cells = "NHBE_BR1_Baseline"
      
      col_ix = which(colnames(dat[[1]]$res) %in% genes)
      # print(col_ix)
      
      # slice matrices if necessary
      dat_plot = dat
      for(j in 1:length(dat_plot)) { # each data type
        dat_plot[[j]]$res = dat[[j]]$res[,col_ix,drop=FALSE]
      }
      
      single_labels = rownames(dat_plot[[1]]$res)
      
      c_ix = match(c_cells, rownames(dat[[1]]$res))
      # c_ix = c(1:14,18:21)
      t_ix = match(t_cells, rownames(dat[[1]]$res))
      
      # print(t_ix)
      # print(c_ix)
      
      res = dist_mat(dat_plot, comp_ix=list(c_ix, t_ix), labels=single_labels, plot_labels=c("BEAS2B","A549","NHLF"), plot_res=TRUE, use_corr=TRUE, font_size=25, label_size=input$label.size)
      
    } else {
      res = dist_mat(x=dat, comp_ix=list(), labels=single_labels, plot_res=TRUE, use_corr=TRUE, font_size=25, label_size=input$label.size, plot_blank=TRUE)
    }
  })
  
  output$dist_plot_2 <- renderPlot({
    qplot(1,1) + theme_thesis(20)
  })
  
  output$plot1 <- renderPlot({
    
    dat_out = spotfire_view(dat)
    
    if(my_choices()$gene_choice[1]!="NULL" | my_choices()$go_choice[1]!="NULL") {
      if(my_choices()$go_choice[1]!="NULL") {
        genes = msig_go_bp[[which(names(msig_go_bp)==my_choices()$go_choice)]]
      } else {
        genes = my_choices()$gene_choice
      }
      dat_out$color = ifelse(dat_out$Gene %in% genes, 2, 1)
    } else {
      dat_out$color = 1
    }
    
    dat_out$color = factor(dat_out$color)
    dat_out = tbl_df(dat_out)
    ggplot(dat_out, aes_string(names(dat_out)[2], names(dat_out)[3])) + geom_point(size=2, aes(color=color, alpha=color)) + theme_thesis() + scale_color_discrete(guide=FALSE) + scale_alpha_manual(guide=FALSE, values=c(0.3,0.8))
    
  })
  
  output$click_info <- renderPrint({
    # Because it's a ggplot2, we don't need to supply xvar or yvar; if this
    # were a base graphics plot, we'd need those.
    nearPoints(dat_out, input$plot1_click, addDist=TRUE)
  })
  
  output$brush_info <- renderPrint({
    brushedPoints(dat_out, input$plot1_brush)
  })
  
  output$download_table <- downloadHandler("selected_genes.csv", content=function(file) {
    write.csv(out$brush_info, file, quote=TRUE, row.names=FALSE)
  })
  
  output$sushi<- renderPlot({
    
    sample_ix = c(1,3,15)
    data_type = 1
    win = 20000
    
    if(
      (my_choices()$gene_choice[1]!="NULL" | my_choices()$go_choice[1]!="NULL") &
      my_choices()$cell_target_choice[1]!="NULL" &
      my_choices()$cell_candidate_choice[1]!="NULL"
    ) {
      
      if(my_choices()$go_choice[1]!="NULL") {
        print(my_choices()$go_choice)
        genes = msig_go_bp[[which(names(msig_go_bp)==my_choices()$go_choice)]]
      } else {
        genes = my_choices()$gene_choice
      }
      
      col_ix = which(colnames(dat[[1]]$res) %in% genes)
      
      roi = gene_list_all[col_ix]
      start(roi) = start(roi) - win
      end(roi) = end(roi) + win
      
      my_tracks = sapply(str_replace(dat[[data_type]]$annot$Bigwig[sample_ix], "/GWD/bioinfo/projects/", "z:/links/"), function(x) import.bw(x, which=roi))
      
      mart_1 = useMart("ensembl", dataset="hsapiens_gene_ensembl")
      t_list = getBM(attributes=c("chromosome_name","exon_chrom_start","exon_chrom_end","ensembl_transcript_id","strand","ensembl_gene_id"), filters='hgnc_symbol', values=roi$hgnc_symbol, mart=mart_1)
      t_list$type = "exon"
      t_list$chromosome_name = paste0("chr", t_list$chromosome_name)
      
      my_tracks_df = lapply(my_tracks, function(x) as.data.frame(x)[,c(1:3,6)])
      
      par(mfcol=c(length(sample_ix)+1,length(roi)), mar=c(4,4,2,2))
      
      for(g_ix in 1:length(roi)) {
        
        chrom = as.character(seqnames(roi)[g_ix])
        chromstart = start(roi)[g_ix]
        chromend = end(roi)[g_ix]
        
        for(i in 1:length(sample_ix)) {
          
          if(g_ix==1) {
            plotBedgraph(my_tracks_df[[i]], chrom, chromstart, chromend, transparency=.2, color=SushiColors(2)(length(sample_ix))[i], main=rownames(dat[[1]]$res)[sample_ix[i]])
          } else{
            plotBedgraph(my_tracks_df[[i]], chrom, chromstart, chromend, transparency=.2, color=SushiColors(2)(length(sample_ix))[i])
          }
          labelgenome(chrom, chromstart, chromend, n=10, scale="Mb")
          axis(side=2, las=2, tcl=.2)
          
        }
        
        plotGenes(t_list, chrom, chromstart, chromend, types=t_list$type, labeltext=TRUE, maxrows=50, height=0.4, plotgenetype="box")
        labelplot(title=roi$hgnc_symbol[g_ix])
        
      }
    }
  })
  
})

