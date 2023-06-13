library(shiny)
library(shinyWidgets)
library(png)
library(EBImage)
library(DT)

# Create a list of colors for the color selector
colors <- rainbow(64)

ui <- fluidPage(
  # Inject custom CSS
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
  fileInput('file', 'Choose PNG Image',
            accept = c('image/png')),
  colorSelectorInput('col', 'Choose a color', choices = colors),
  textInput('label', 'Enter label for color'),
  imageOutput("original_img"),
  plotOutput("grayscale_img", click = "grayscale_img_click"),
  actionButton("save", "Save Image"),
  DTOutput("label_table")
)

server <- function(input, output, session) {
  
  colorLabelStore <- reactiveValues(color=NULL, label=NULL)
  
  observeEvent(input$file, {
    req(input$file)
    img <- readPNG(input$file$datapath)
    # img should be processed here...
    output$original_img <- renderImage({
      list(src = input$file$datapath,
           contentType = input$file$type,
           width = 400)
    })
    # grayscale image here is just the original image for demonstration
    # it should be the processed grayscale image instead
    output$grayscale_img <- renderPlot({
      plot(as.raster(img), axes = FALSE)
    })
  })
  
  observeEvent(input$grayscale_img_click, {
    # update color and label
    colorLabelStore$color <- input$col
    colorLabelStore$label <- input$label
    # process the clicked segment with colorLabelStore$color and colorLabelStore$label here...
  })
  
  observeEvent(input$save, {
    # save the image and generate the dataframe here...
    label_table <- data.frame(
      # just an example dataframe
      # this should be generated according to your specific needs
      Segment = c(1, 2, 3),
      Label = c("label1", "label2", "label3")
    )
    output$label_table <- renderDT({
      datatable(label_table)
    })
  })
  
}

shinyApp(ui, server)
