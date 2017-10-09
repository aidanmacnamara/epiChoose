require(shiny)
require(stringr)

load("dat.RData")
load("msig_go_bp.RData")

cells = rownames(dat[[1]]$res)
genes = colnames(dat[[1]]$res)
# names(msig_go_bp) = str_replace(str_replace_all(names(msig_go_bp), "_", " "), "GO\\s+", "")

ui <- fluidPage(
  
  titlePanel("epiChoose: Data-driven cell model choice"),
  
  sidebarLayout(
    
    sidebarPanel(
      
      selectInput("cell_target_choice", label="Cell Target Choice", 
                  choices = as.list(sort(cells)),
                  multiple=TRUE),
      
      helpText(
        "Choose the target cell type"
      ),
      
      h3(""),
      h3(""),
      
      selectInput("cell_candidate_choice", label="Cell Model Choice", 
                  choices = as.list(sort(cells)),
                  multiple=TRUE),
      
      helpText(
        "Choose the possible cell models"
      ),
      
      h3(""),
      h3(""),
      
      selectInput("gene_choice", label="Gene Choice", 
                  choices = as.list(sort(genes)),
                  multiple=TRUE),
      
      helpText(
        "Choose the genes over which to make the comparison"
      ),
      
      h3(""),
      h3(""),
      
      selectInput("go_choice", label="Gene Ontology Choice", 
                  choices = as.list(sort(names(msig_go_bp))),
                  multiple=TRUE),
      
      helpText(
        "Choose the GO biological process over which to make the comparison"
      ),
      
      h3(""),
      h3(""),
      
      selectInput("data_type_choice", label="Data Type Choice", 
                  choices = as.list(sort(names(dat))),
                  multiple=FALSE),
      
      helpText(
        "Choose which data type to use"
      ),
      
      width=3
    ),
    
    mainPanel(
      tabsetPanel(
        
        id = "tabs",
        
        tabPanel("Gene Choice",
                 
                 h3(""),
                 h3(""),
                 
                 column(width=6,
                        plotOutput("dist_plot_1", height=600),
                        sliderInput("label.size", NULL, min=1, max=10, value=4, step=1)
                 ),
                 column(width=6,
                        plotOutput("dist_plot_2", height=600)
                 )
        ),
        
        tabPanel("Model Choice",
                 column(width=6,
                        plotOutput("plot1", height=600,
                                   click = "plot1_click",
                                   brush = brushOpts(id="plot1_brush")
                        ),
                        
                        h3(""),
                        h3(""),
                        
                        column(width=4,     
                               selectInput("candidate_choice", label="Candidate Choice", 
                                           choices = as.list(LETTERS[1:3]),
                                           multiple=TRUE),
                               
                               helpText(
                                 "Choose the model candidate"
                               )
                        ),
                        
                        column(width=4,     
                               selectInput("alt_choice", label="Alternative Choice", 
                                           choices = as.list(LETTERS[4:6]),
                                           multiple=TRUE),
                               
                               helpText(
                                 "Choose the model alternative"
                               )
                        ),
                        
                        column(width=2,     
                               selectInput("x_axis", label="X-axis", 
                                           choices = as.list(LETTERS[4:6]),
                                           multiple=TRUE),
                               
                               helpText(
                                 "Choose the data type for the x-axis"
                               )
                        ),
                        
                        column(width=2,     
                               selectInput("y_axis", label="Y-axis", 
                                           choices = as.list(LETTERS[4:6]),
                                           multiple=TRUE),
                               
                               helpText(
                                 "Choose the data type for the y-axis"
                               )
                        )
                 ),
                 column(width=6,
                        h4("Brushed points"),
                        verbatimTextOutput("brush_info")
                 ),
                 downloadButton(
                   'downloadData',
                   'Download'
                 )
        ),
        
        tabPanel("Gene Browser",
                 plotOutput("sushi", height=800)
        )
      ),
      
      width=9
      
    )
  )
)

