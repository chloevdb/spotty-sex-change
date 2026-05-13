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
library(DESeq2)
library(here)

# ── Data prep (runs once on startup) ──────────────────────────────────────────

counts.table <- read.table(here("data-files/readspergene_v2.matrix"),
                           header = TRUE, sep = "\t",
                           stringsAsFactors = FALSE, row.names = "GeneID")

coldata.table <- read.table(here("data-files/coldata_v2.3.2.txt"),
                            header = TRUE, sep = "\t",
                            stringsAsFactors = FALSE, row.names = "sample")

dds <- DESeqDataSetFromMatrix(countData = counts.table,
                              colData   = coldata.table,
                              design    = ~ histov22)
dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind = TRUE)

vsd_normalized_counts <- as.data.frame(assay(vsd))
vsd_normalized_counts$GeneID <- rownames(vsd_normalized_counts)
rownames(vsd_normalized_counts) <- NULL
vsd_normalized_counts <- vsd_normalized_counts[, c(97, 1:96)]

data.long <- pivot_longer(vsd_normalized_counts,
                          cols = -GeneID,
                          names_to  = "Sample",
                          values_to = "expression")

data.long.coldata <- merge(data.long, coldata.table,
                           by.x = "Sample", by.y = "samplenames")

custom_order  <- c("F", "ET", "MT", "LT", "TPM", "IPM")
custom_colors <- c("#F564E3", "#984EA3", "#00BFC4",
                   "#4DAF4A", "#619CFF",  "#CD9600")

data.long.coldata$histov22 <- factor(data.long.coldata$histov22,
                                     levels = custom_order)

# ── Gene list  ─────────────────────
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
      plotOutput("boxplot", height = "600px")
    )
  )
)

# ── Server ────────────────────────────────────────────────────────────────────

server <- function(input, output, session) {
  
  make_plot <- reactive({
    req(input$genes)
    data.long.coldata %>%
      filter(GeneID %in% input$genes) %>%
      mutate(histov22 = factor(histov22, levels = custom_order)) %>%
      ggplot(aes(x = histov22, y = expression, fill = histov22)) +
      geom_boxplot(size = 0.3, outlier.size = 0.3) +
      scale_fill_manual(values = custom_colors) +
      labs(y = "Variance stabilised counts", x = "Sample") +
      facet_wrap(~ GeneID, scales = "free") +
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

