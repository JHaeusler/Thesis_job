library(shiny)
library(DT)

ui <- fluidPage(
  titlePanel("Formulario con Matriz de Selección"),
  
  sidebarLayout(
    sidebarPanel(
      helpText("Selecciona una celda de la matriz para registrar tu respuesta."),
      hr(),
      verbatimTextOutput("info_seleccion")
    ),
    
    mainPanel(
      tags$h4("Matriz de Decisión"),
      DTOutput("matriz_input")
    )
  )
)

server <- function(input, output, session) {
  
  # 1. Definimos los datos de la matriz
  df_matriz <- data.frame(
    Bajo = c("A1", "B1", "C1"),
    Medio = c("A2", "B2", "C2"),
    Alto = c("A3", "B3", "C3"),
    row.names = c("Fila A", "Fila B", "Fila C")
  )
  
  # 2. Renderizamos la tabla interactiva
  output$matriz_input <- renderDT({
    datatable(
      df_matriz,
      selection = list(mode = 'single', target = 'cell'), # Selección única de celda
      options = list(dom = 't', ordering = FALSE),        # Oculta buscador y filtros
      class = 'cell-border stripe'
    )
  })
  
  # 3. Capturamos la selección del usuario
  output$info_seleccion <- renderPrint({
    sel <- input$matriz_input_cells_selected # Devuelve [fila, columna]
    
    if (length(sel) == 0) {
      cat("No se ha seleccionado ninguna combinación.")
    } else {
      # Extraemos nombres de fila y columna
      fila_idx <- sel[1, 1]
      col_idx  <- sel[1, 2] # Nota: la columna 0 suele ser el nombre de la fila en DT
      
      nombre_fila <- rownames(df_matriz)[fila_idx]
      nombre_col  <- colnames(df_matriz)[col_idx]
      
      cat(paste("Seleccionaste:\nFila:", nombre_fila, "\nColumna:", nombre_col))
    }
  })
}

shinyApp(ui, server)