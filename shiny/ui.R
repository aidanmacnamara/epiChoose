library(shiny)

# load("dat.RData")
# load_all("../")
cells = rownames(dat[[1]]$res)
genes = colnames(dat[[1]]$res)


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
      
      width=3
    ),
    
    mainPanel(
      tabsetPanel(
        id = "tabs",
        tabPanel("Gene Choice", plotOutput("dist_plot", height=800)),
        tabPanel("Model Choice",
                 column(width=6,
                        plotOutput("plot1", height=600,
                                   click = "plot1_click",
                                   brush = brushOpts(id="plot1_brush")
                        )
                 ),
                 column(width=6,
                        h4("Brushed points"),
                        verbatimTextOutput("brush_info")
                 )
        )
      ),
      width=9
    )
  )
)

