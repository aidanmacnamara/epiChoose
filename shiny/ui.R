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
                  multiple=TRUE,
                  selected=list("NHBE_BR1_Baseline")),
      
      bsTooltip("cell_target_choice", "something", placement="bottom", trigger="hover", options=NULL),
      
      helpText(
        "Choose the target cell type"
      ),
      
      h3(""),
      h3(""),
      
      selectInput("cell_candidate_choice", label="Cell Model Choice", 
                  choices = as.list(sort(cells)),
                  multiple=TRUE,
                  selected=list(
                    "A549_BR1_Baseline",
                    "A549_BR2_Baseline",
                    "BEAS2B_BR1_Baseline",
                    "BEAS2B_BR2_Baseline",
                    "A549_Broad"
                  )),
      
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
      
      width=2
    ),
    
    mainPanel(
      tabsetPanel(
        
        id = "tabs",
        
        tabPanel("Gene Choice",
                 
                 h3(""),
                 h3(""),
                 
                 column(width=5,
                        plotOutput("dist_plot_1", height=600),
                        sliderInput("label.size.1", NULL, min=1, max=10, value=4, step=1)
                 ),
                 column(width=7,
                        plotOutput("dist_plot_2", height=800),
                        sliderInput("label.size.2", NULL, min=1, max=10, value=4, step=1)
                 )
        ),
        
        tabPanel("Model Choice",
                 fluidRow(
                   column(width=7,
                          plotOutput("plot1", height=600,
                                     click = "plot1_click",
                                     brush = brushOpts(id="plot1_brush")
                          )
                   ),
                   
                   column(width=5,
                          h4("Brushed points"),
                          verbatimTextOutput("brush_info"),
                          
                          downloadButton(
                            'download_table',
                            'Download'
                          )
                   )
                 ),
                 
                 h3(""),
                 h3(""),
                 
                 column(width=3,   
                        
                        selectInput("target_choice", label="Target Choice", 
                                    choices = as.list(sort(cells)),
                                    multiple=TRUE,
                                    selected=list("NHBE_BR1_Baseline")),
                        
                        
                        selectInput("x_axis", label="X-axis", 
                                    choices = as.list(sort(names(dat))),
                                    multiple=FALSE,
                                    selected=list("H3K27ac"))
                 ),
                 
                 column(width=3,
                        
                        selectInput("candidate_choice", label="Candidate Choice", 
                                    choices = as.list(sort(cells)),
                                    multiple=TRUE,
                                    selected=list("A549_BR1_Baseline")),
                        
                        selectInput("y_axis", label="Y-axis", 
                                    choices = as.list(sort(names(dat))),
                                    multiple=FALSE,
                                    selected=list("RNA"))
                 ),
                 
                 
                 column(width=3,
                        
                        selectInput("alt_choice", label="Alternative Choice", 
                                    choices = as.list(sort(cells)),
                                    multiple=TRUE,
                                    selected=list("BEAS2B_BR1_Baseline")),
                        
                        actionButton("do_model", "Run")
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

