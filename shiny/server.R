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

single_labels = rownames(dat[[1]]$res)
group_labels = "GSK"
pca_data = prep_for_plot(dat, annot_1=group_labels, annot_2=single_labels, marks=names(dat), plot_type="mds")

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
  
  model_choice <- eventReactive(input$do_model,
    {
      
      model_choice = list()
      
      c_choice = renderPrint({input$candidate_choice})()
      c_choice = unlist(regmatches(c_choice, gregexpr("[[:alnum:]_]+", c_choice)))[-1] # remove na added
      c_ix = which(rownames(dat[[1]]$res) %in% c_choice)
      
      a_choice = renderPrint({input$alt_choice})()
      a_choice = unlist(regmatches(a_choice, gregexpr("[[:alnum:]_]+", a_choice)))[-1] # remove na added
      a_ix = which(rownames(dat[[1]]$res) %in% a_choice)
      
      t_choice = renderPrint({input$target_choice})()
      t_choice = unlist(regmatches(t_choice, gregexpr("[[:alnum:]_]+", t_choice)))[-1] # remove na added
      t_ix = which(rownames(dat[[1]]$res) %in% t_choice)
      
      x_axis_choice <- renderPrint({input$x_axis})()
      x_axis_choice = str_replace(x_axis_choice, "^.*\"([[:alnum:]\\s]+)\"", "\\1") 
      x_axis_ix = which(names(dat)==x_axis_choice)
      
      x_axis_choice <- renderPrint({input$x_axis})()
      x_axis_choice = str_replace(x_axis_choice, "^.*\"([[:alnum:]\\s]+)\"", "\\1") 
      
      y_axis_choice <- renderPrint({input$y_axis})()
      y_axis_choice = str_replace(y_axis_choice, "^.*\"([[:alnum:]\\s]+)\"", "\\1")
      
      if(
        length(c_ix) &
        length(a_ix) &
        length(t_ix) &
        x_axis_choice!="NULL" &
        y_axis_choice!="NULL"
      )  {
        
        model_choice$my_df = spotfire_view(dat, x_axis=x_axis_choice, y_axis=y_axis_choice, comp_ix=list(c_ix, a_ix, t_ix))
        rownames(model_choice$my_df) = NULL
        model_choice$my_df = tbl_df(model_choice$my_df)
      }
      
      return(model_choice)
      
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
      
      c_cells = my_choices()$cell_candidate_choice
      t_cells = my_choices()$cell_target_choice
      
      col_ix = which(colnames(dat[[1]]$res) %in% genes)
      
      # slice matrices if necessary
      dat_plot = dat
      for(j in 1:length(dat_plot)) { # each data type
        dat_plot[[j]]$res = dat[[j]]$res[,col_ix,drop=FALSE]
      }
      
      single_labels = rownames(dat_plot[[1]]$res)
      
      c_ix = match(c_cells, rownames(dat[[1]]$res))
      t_ix = match(t_cells, rownames(dat[[1]]$res))
      
      res = dist_mat(dat_plot, comp_ix=list(c_ix, t_ix), labels=single_labels, plot_labels=c("BEAS2B","A549","NHLF"), plot_res=TRUE, use_corr=TRUE, font_size=15, label_size=input$label.size.1)
      
    } else {
      res = dist_mat(x=dat, comp_ix=list(), labels=single_labels, plot_res=TRUE, use_corr=TRUE, font_size=15, label_size=input$label.size.1, plot_blank=TRUE)
    }
  })
  
  output$dist_plot_2 <- renderPlot({
    
    print(ggplot(pca_data, aes(x=x, y=y, color=annot_1)) + geom_point(size=3, shape=17) + theme_thesis(20) + geom_text_repel(aes(label=annot_2), fontface="bold", size=input$label.size.2, force=0.5) + facet_wrap(~mark, ncol=2, scales="free"))
    
  })
  
  output$plot1 <- renderPlot({
    
    dat_out = model_choice()$my_df
    
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
  
  selected_data <- reactive({
    brushedPoints(model_choice()$my_df, input$plot1_brush)
  })
  
  output$brush_info <- renderPrint({
    selected_data()
  })
  
  output$download_table <- downloadHandler(
    filename = function() { 
      paste('selected_rows', '.csv', sep='') 
    },
    content = function(file) {
      write_csv(selected_data(), file)
    }
  )
  
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

