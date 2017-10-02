library(shiny)
mtcars2 <- mtcars[, c("mpg", "cyl", "disp", "hp", "wt", "am", "gear")]

shinyServer(function(input, output) {
  
  my_choices <- reactive(
    {
      
      my_choices = list()
      
      gene_choice <- renderPrint({input$gene_choice})()
      my_choices$gene_choice = unlist(regmatches(gene_choice, gregexpr("[A-Z][[:alnum:]]+", gene_choice)))
      
      cell_target_choice <- renderPrint({input$cell_target_choice})()
      my_choices$cell_target_choice = unlist(regmatches(cell_target_choice, gregexpr("[[:alnum:]_]+", cell_target_choice)))[-1] # remove na added
      
      cell_candidate_choice <- renderPrint({input$cell_candidate_choice})()
      my_choices$cell_candidate_choice = unlist(regmatches(cell_candidate_choice, gregexpr("[[:alnum:]_]+", cell_candidate_choice)))[-1] # remove na added
      
      return(my_choices)
    }
  )
  
  output$dist_plot <- renderPlot({
    
    if(
      my_choices()$gene_choice[1]!="NULL" &
      my_choices()$cell_target_choice[1]!="NULL" &
      my_choices()$cell_candidate_choice[1]!="NULL"
    ) {
      
      genes = my_choices()$gene_choice
      # genes = go_genes[[2]]
      c_cells = my_choices()$cell_candidate_choice
      # c_cells = c("A549_BR1_Baseline","A549_BR2_Baseline","BEAS2B_BR1_Baseline","BEAS2B_BR2_Baseline")
      t_cells = my_choices()$cell_target_choice
      # t_cells = "NHBE_BR1_Baseline"
      
      col_ix = which(gene_list_all$hgnc_symbol %in% genes)
      print(col_ix)
      
      # slice matrices if necessary
      dat_plot = dat
      for(j in 1:length(dat_plot)) { # each data type
        dat_plot[[j]]$res = dat[[j]]$res[,col_ix,drop=FALSE]
      }
      
      single_labels = rownames(dat_plot[[1]]$res)
      
      c_ix = match(c_cells, rownames(dat[[1]]$res))
      # c_ix = c(1:14,18:21)
      t_ix = match(t_cells, rownames(dat[[1]]$res))
      
      print(t_ix)
      print(c_ix)
      
      res = dist_mat(dat_plot, comp_ix=list(c_ix, t_ix), labels=single_labels, plot_labels=c("BEAS2B","A549","NHLF"), plot_res=TRUE, use_corr=TRUE, font_size=30)
      
    }
  })
  
  output$plot1 <- renderPlot({
    ggplot(dat_out, aes_string(names(dat_out)[2], names(dat_out)[3])) + geom_point(alpha=0.2) + theme_thesis()
  })
  
  output$click_info <- renderPrint({
    # Because it's a ggplot2, we don't need to supply xvar or yvar; if this
    # were a base graphics plot, we'd need those.
    nearPoints(dat_out, input$plot1_click, addDist=TRUE)
  })
  
  output$brush_info <- renderPrint({
    brushedPoints(dat_out, input$plot1_brush)
  })
  
})

