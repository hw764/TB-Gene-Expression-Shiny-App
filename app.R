
library(shiny)

library(tidyverse)

library(mclust)

# Load the TB Nanostring dataset
TBnanostring <- readRDS("TBnanostring.rds")

# Separate gene expression measurements from TB status
gene_data <- TBnanostring %>%
  select(-TB_Status)

# Save the actual TB status labels
tb_status <- TBnanostring$TB_Status

# Check that the data loaded correctly
print(dim(gene_data))
print(table(tb_status))

# Load packages for dimensionality reduction
library(uwot)
library(NMF)

# Convert gene expression data to a matrix
gene_matrix <- as.matrix(gene_data)

# PCA: reduce the 107 genes to 2 components
pca_result <- prcomp(gene_matrix, center = TRUE, scale. = FALSE)

pca_coords <- pca_result$x[, 1:2]

# UMAP: reduce the 107 genes to 2 dimensions
set.seed(123)
umap_coords <- uwot::umap(gene_matrix, n_components = 2)

# NMF requires non-negative data
# Shift each gene so its minimum is zero
nmf_matrix <- sweep(
  gene_matrix,
  2,
  apply(gene_matrix, 2, min),
  FUN = "-"
)

# NMF: reduce the data to 2 components
set.seed(123)
nmf_result <- NMF::nmf(nmf_matrix, rank = 2, nrun = 10)

nmf_coords <- NMF::basis(nmf_result)

# Check the dimensions
print(dim(pca_coords))
print(dim(umap_coords))
print(dim(nmf_coords))



ui <- fluidPage(
  
  titlePanel("TB Gene Expression Clustering App"),
  
  sidebarLayout(
    
    sidebarPanel(
      
      # Choose dimensionality reduction method
      selectInput(
        "method",
        "Dimension Reduction Method",
        choices = c("PCA", "UMAP", "NMF"),
        selected = "PCA"
      ),
      
      # Choose number of K-means clusters
      numericInput(
        "clusters",
        "Cluster Count (k)",
        value = 2,
        min = 2,
        max = 10
      )
      
    ),
    
    mainPanel(
      plotOutput("plot1"),
      
      h4("Clusters vs. TB Status"),
      tableOutput("cluster_table"),
      
      h4("Adjusted Rand Index"),
      textOutput("ari_value")
    )
    
  )
  
)


server <- function(input, output, session) {
  
  # Choose which dimensionality reduction results to use
  selectedData <- reactive({
    
    if (input$method == "PCA") {
      coords <- pca_coords
      
    } else if (input$method == "UMAP") {
      coords <- umap_coords
      
    } else {
      coords <- nmf_coords
    }
    
    # Create a dataset containing the coordinates and TB status
    data.frame(
      Dim1 = coords[, 1],
      Dim2 = coords[, 2],
      TB_Status = tb_status
    )
    
  })
  
  # Run K-means on the selected two dimensions
  clusters <- reactive({
    
    set.seed(123)
    
    kmeans(
      selectedData()[, c("Dim1", "Dim2")],
      centers = input$clusters,
      nstart = 25
    )
    
  })
  
  # Plot the reduced dimensions
  output$plot1 <- renderPlot({
    
    plot_data <- selectedData()
    
    # Add cluster assignments to the dataset
    plot_data$Cluster <- factor(clusters()$cluster)
    
    # Create the scatterplot
    ggplot(
      plot_data,
      aes(
        x = Dim1,
        y = Dim2,
        color = Cluster,
        shape = TB_Status
      )
    ) +
      geom_point(size = 2.5, alpha = 0.8) +
      
      # Add cluster centers as large X symbols
      geom_point(
        data = as.data.frame(clusters()$centers),
        aes(x = Dim1, y = Dim2),
        inherit.aes = FALSE,
        shape = 4,
        size = 5,
        stroke = 1.5,
        color = "black"
      ) +
      
      labs(
        title = paste(input$method, "and K-means Clustering"),
        x = paste(input$method, "Dimension 1"),
        y = paste(input$method, "Dimension 2"),
        color = "K-means Cluster",
        shape = "TB Status"
      ) +
      
      theme_minimal()
    
  })
  
  # Create a table comparing clusters with TB status
  output$cluster_table <- renderTable({
    
    table(
      Cluster = clusters()$cluster,
      TB_Status = tb_status
    )
    
  }, rownames = TRUE)
  
  # Calculate the Adjusted Rand Index
  output$ari_value <- renderText({
    
    ari <- mclust::adjustedRandIndex(
      clusters()$cluster,
      tb_status
    )
    
    paste("ARI =", round(ari, 4))
    
  })
  
}


shinyApp(ui = ui, server = server)

