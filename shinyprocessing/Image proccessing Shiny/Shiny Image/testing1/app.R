library(shiny)
library(shinyWidgets)
library(EBImage)
library(magrittr)
library(magick)
library(DT)

colors <- rainbow(64)

ui <- fluidPage(
  tags$style("
    .select2-results__option {
      width: 30px;
      height: 30px;
      display: inline-block;
      margin-right: 5px;
      margin-bottom: 5px;
    }
    .select2-results__group {
      display: grid;
      grid-template-columns: repeat(auto-fill, minmax(30px, 1fr));
    }
  "),
  sidebarLayout(
    sidebarPanel(
      fileInput('file', 'Choose PNG Image',
                accept = c('image/png')),
      colorSelectorInput('col', 'Choose a color', choices = colors),
      textInput('label', 'Enter label for color'),
      actionButton("save", "Save Label"),
      downloadButton("download_img", "Download Image"),
      downloadButton("download_data", "Download Data"),
      DTOutput("label_table")
    ),
    mainPanel(
      imageOutput("original_img"),
      plotOutput("grayscale_img", click = "grayscale_img_click")
    )
  )
)

server <- function(input, output, session) {
  savedColorAndLabel <- reactiveValues(color = NULL, label = NULL, data = data.frame())
  
  observeEvent(input$save, {
    savedColorAndLabel$color <- input$col
    savedColorAndLabel$label <- input$label
    new_row <- data.frame(color = input$col, label = input$label)
    savedColorAndLabel$data <- rbind(savedColorAndLabel$data, new_row)
    output$label_table <- renderDT(savedColorAndLabel$data, options = list(pageLength = 25))
  })
  
  img_data <- reactive({
    req(input$file)
    validate(need(tools::file_ext(input$file$datapath) == "png", "Please upload a PNG file."))
    
    # Use EBImage to read and process the image
    img <- EBImage::readImage(input$file$datapath)
    
    if (length(dim(img)) == 2 || dim(img)[3] == 1) {
      print("The image is already grayscale.")
    } else {
      print("Converting the image to grayscale.")
      img <- channel(img, 'gray')
    }
    
    bw_img <- normalize(img)
    threshold <- otsu(bw_img)
    thresh_img <- bw_img > threshold
    dist_map <- distmap(thresh_img)
    ws <- watershed(dist_map)
    ws_binary <- ifelse(ws > 1, 1, 0)
    cc <- bwlabel(ws_binary)
    cc
  })
  
  output$original_img <- renderImage({
    img <- EBImage::readImage(input$file$datapath)
    cc <- img_data()
    EBImage::display(EBImage::colorLabels(cc))
    list(src = input$file$datapath,
         contentType = 'image/png',
         alt = "Original Image",
         width = 400,
         height = 400)
  }, deleteFile = FALSE)
  
  output$grayscale_img <- renderPlot({
    plot(as.raster(img_data()), axes = FALSE)
  })
  
  output$download_img <- downloadHandler(
    filename = function() {
      paste("labeled_image", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      png(file)
      plot(as.raster(img_data()), axes = FALSE)
      dev.off()
    }
  )
  
  output$download_data <- downloadHandler(
    filename = function() {
      paste("color_label_data", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(savedColorAndLabel$data, file, row.names = FALSE)
    }
  )
}

shinyApp(ui = ui, server = server)
