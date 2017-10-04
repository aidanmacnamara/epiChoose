library(shiny)

# load("dat.RData")

# load_all("../")
cells = rownames(dat[[1]]$res)
genes = colnames(dat[[1]]$res)
names(msig_go_bp) = str_replace(str_replace_all(names(msig_go_bp), "_", " "), "GO\\s+", "")


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
      
      width=3
    ),
    
    mainPanel(
      tabsetPanel(
        id = "tabs",
        tabPanel("Gene Choice",
                 plotOutput("dist_plot", height=800),
                 sliderInput("label.size", NULL, min=1, max=10, value=4, step=1)
        ),
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
                 ),
                 downloadButton(
                   'downloadData',
                   'Download'
                 )
        )
      ),
      width=9
    )
  )
)

