# ── Install missing packages if needed ────────────────────────────────────────

required_packages <- c("shiny", "ggplot2", "tidyr", "dplyr", "here")

new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if (length(new_packages) > 0) install.packages(new_packages)

# DESeq2 is from Bioconductor so needs separate handling
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2")

# ── Load packages ─────────────────────────────────────────────────────────────

library(shiny)
library(ggplot2)
library(tidyr)
library(dplyr)
library(readr)
library(here)

# ── Data prep (runs once on startup) ──────────────────────────────────────────

coldata.table <- read.table(here("data-files/coldata_v2.3.2.txt"),
                            header = TRUE, sep = "\t",
                            stringsAsFactors = FALSE, row.names = "sample")

vsd_normalized_counts <- read_tsv("https://raw.githubusercontent.com/chloevdb/spotty-sex-change/main/shiny-app-gene-expn/data-files/vsd_normalized_counts.tsv")

data.long <- pivot_longer(vsd_normalized_counts,
                          cols = -GeneID,
                          names_to  = "Sample",
                          values_to = "expression")

data.long.coldata <- merge(data.long, coldata.table,
                           by.x = "Sample", by.y = "samplenames")

custom_orderv22  <- c("F", "ET", "MT", "LT", "TPM", "IPM")
custom_colorsv22 <- c("#F564E3", "#984EA3", "#00BFC4",
                      "#4DAF4A", "#619CFF", "#CD9600")

custom_orderv2  <- c("NBF", "BF", "ET", "MT", "LT", "TPM", "IPM")
custom_colorsv2 <- c("#F564E3", "#C77CFF", "#984EA3", "#00BFC4",
                     "#4DAF4A", "#619CFF", "#CD9600")

data.long.coldata$histov22 <- factor(data.long.coldata$histov22,
                                     levels = custom_orderv22)

# ── Gene list ─────────────────────────────────────────────────────────────────

gene_list <- unique(data.long.coldata$GeneID)

# ── UI ────────────────────────────────────────────────────────────────────────

ui <- fluidPage(
  titlePanel("Spotty sex-change gonad gene expression browser"),
  sidebarLayout(
    sidebarPanel(
      selectInput(
        inputId  = "genes",
        label    = "Select gene(s)",
        choices  = gene_list,
        selected = gene_list[1],
        multiple = TRUE
      ),
      helpText("Hold Ctrl / ⌘ to select multiple genes."),
      checkboxInput(
        inputId = "free_scales",
        label   = "Free axis scales (each gene scaled independently)",
        value   = TRUE
      ),
      checkboxGroupInput(
        inputId  = "selected_groups",
        label    = "Sample groups to display",
        choices  = custom_orderv22,
        selected = custom_orderv22
      ),
      checkboxInput(
        inputId = "split_female",
        label   = "Split F into NBF and BF",
        value   = FALSE
      ),
      checkboxGroupInput(
        inputId  = "plot_type",
        label    = "Plot type",
        choices  = c("Boxplot" = "boxplot", "Violin" = "violin"),
        selected = "boxplot"
      ),
      sliderInput(
        inputId = "ncol",
        label   = "Number of columns",
        min     = 1, max = 6, value = 3, step = 1
      ),
      sliderInput(
        inputId = "plot_height",
        label   = "Plot height (px)",
        min     = 300, max = 1200, value = 600, step = 50
      ),
      hr(),
      selectInput(
        inputId  = "file_format",
        label    = "Download format",
        choices  = c("PNG" = "png", "PDF" = "pdf", "SVG" = "svg"),
        selected = "png"
      ),
      downloadButton("download_plot", "Download plot")
    ),
    mainPanel(
      uiOutput("boxplot_ui")
    )
  )
)

# ── Server ────────────────────────────────────────────────────────────────────

server <- function(input, output, session) {
  
  plot_data <- reactive({
    df <- data.long.coldata %>%
      filter(GeneID %in% input$genes) %>%
      filter(histov22 %in% input$selected_groups)
    
    if (input$split_female) {
      df <- df %>%
        mutate(plot_group = ifelse(histov22 == "F", histov2, as.character(histov22)))
      group_order  <- custom_orderv2
      group_colors <- custom_colorsv2
    } else {
      df <- df %>%
        mutate(plot_group = as.character(histov22))
      group_order  <- custom_orderv22
      group_colors <- custom_colorsv22
    }
    
    df$plot_group <- factor(df$plot_group, levels = group_order)
    list(df = df, order = group_order, colors = group_colors)
  })
  
  make_plot <- reactive({
    req(input$genes, input$plot_type, length(input$selected_groups) > 0)
    validate(
      need(length(input$plot_type) > 0, "Please select at least one plot type."),
      need(length(input$selected_groups) > 0, "Please select at least one sample group.")
    )
    
    pd <- plot_data()
    
    pd$df %>%
      ggplot(aes(x = plot_group, y = expression, fill = plot_group)) +
      (if ("violin" %in% input$plot_type) geom_violin(alpha = 0.6) else NULL) +
      (if ("boxplot" %in% input$plot_type) geom_boxplot(size = 0.3, outlier.size = 0.3, width = 0.6) else NULL) +
      scale_fill_manual(values = pd$colors) +
      labs(y = "Variance stabilised counts", x = "Sample") +
      facet_wrap(~ GeneID,
                 scales      = if (input$free_scales) "free" else "fixed",
                 ncol        = input$ncol,
                 axes        = "all",
                 axis.labels = "all") +
      theme_classic() +
      theme(
        axis.line        = element_line(colour = "black"),
        axis.title.y     = element_text(vjust = 3),
        axis.title.x     = element_text(vjust = -1),
        axis.text.x      = element_text(size = 12, angle = 45, hjust = 1, colour = "black"),
        axis.text.y      = element_text(size = 12),
        plot.title       = element_text(size = 18, hjust = 0.5, vjust = 3, face = "bold.italic"),
        strip.text       = element_text(size = 12, face = "italic"),
        strip.background = element_blank(),
        legend.position  = "none"
      )
  })
  
  output$boxplot_ui <- renderUI({
    plotOutput("boxplot", height = paste0(input$plot_height, "px"))
  })
  
  output$boxplot <- renderPlot({ make_plot() })
  
  output$download_plot <- downloadHandler(
    filename = function() {
      genes_str <- paste(input$genes, collapse = "_")
      paste0("gene_expression_", genes_str, ".", input$file_format)
    },
    content = function(file) {
      ggsave(file,
             plot   = make_plot(),
             device = input$file_format,
             width  = 10,
             height = 6,
             dpi    = 300)
    }
  )
}

# ── Run ───────────────────────────────────────────────────────────────────────

shinyApp(ui, server)