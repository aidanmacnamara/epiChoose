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
        tabPanel("Something", plotOutput("dist_plot", height=800))
      ),
      width=9
    )
  )
)
  
  