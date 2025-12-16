#Necessary libraries
library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(pheatmap)
library(ggpubr)
library(tibble)
library(commonmark)
library(matrixStats)
library(ggrepel)
library(DT)
library(ggforce)
if (FALSE) {
  library(munsell)
}


# Find experiments (directories inside "data/")
experiment_dirs <- list.dirs("data", recursive = FALSE)
names(experiment_dirs) <- basename(experiment_dirs)

# Layout of shiny app, we will divide it into tabs (for plots, heatmaps and PCA), each tab has sidebars with select options for Experiment, Genotype, Genes

# UI interface
##### 
ui <- navbarPage("Transcriptomic Explorer",
                 
                 tabPanel("Expression",
                          sidebarLayout(
                            sidebarPanel(
                              selectInput("experiment_line", "Select experiment:", choices = names(experiment_dirs)),
                              uiOutput("genotype_buttons"),
                              selectizeInput("genes", "Select gene(s):", choices = NULL, multiple = TRUE)
                            ),
                            mainPanel(
                              htmlOutput("dataset_info_line"),
                              # tabset panel to divide Expression tab into smaller tabs with these subtabs
                              tabsetPanel(
                                tabPanel("Expression Line Plot",
                                         plotOutput("line_plot"),
                                         textOutput("message")),
                                tabPanel("Expression Box Plot",
                                         plotOutput("box_plot"),
                                         textOutput("message")),
                                tabPanel("Expression Data Table",
                                         downloadButton("download_data", "Download Table"),
                                         DTOutput("line_data_table"))
                              )
                            )
                          )
                 ),
                 
                 tabPanel("Heatmap",
                          sidebarLayout(
                            sidebarPanel(
                              selectInput("experiment_heatmap", "Select experiment:", choices = names(experiment_dirs)),
                              
                              uiOutput("genotype_selector_heatmap"),
                              
                              radioButtons("heatmap_scale", "Heatmap scale:",
                                           choices = c("TPM", "Z-score"), selected = "Z-score", inline = TRUE),
                              
                              radioButtons("gene_selection_mode", "Select genes:",
                                           choices = c("Select manually", "Use named gene list"),
                                           selected = "Select manually"),
                              
                              conditionalPanel(
                                condition = "input.gene_selection_mode == 'Select manually'",
                                selectizeInput("genes_heat", "Select gene(s):", choices = NULL, multiple = TRUE)
                              ),
                              
                              conditionalPanel(
                                condition = "input.gene_selection_mode == 'Use named gene list'",
                                selectInput("gene_set_choice", "Choose gene list:", choices = NULL)
                              )
                            ),
                            mainPanel(
                              htmlOutput("dataset_info_heatmap"),
                              plotOutput("heatmap_plot")
                            )
                          )
                 ),
                 
                 tabPanel("PCA",
                          sidebarLayout(
                            sidebarPanel(
                              selectInput("experiment_pca", "Select experiment:", choices = names(experiment_dirs)),
                              uiOutput("genotype_pca"),
                              uiOutput("treatments_pca"),
                              checkboxInput("show_labels", "Show Labels", value = TRUE)
                            ),
                            mainPanel(
                              plotOutput("pca_plot", height = "600px")
                            )
                          )
                 )
                 
)
 

# SERVER
#####
# Server functions
server <- function(input, output, session) {
  
  
  # Description of experimental conditions for each dataset
  description_table <- read.delim("descriptions.txt", sep = ",")
  
  get_description <- function(dataset_name) {
    row <- description_table[description_table$dataset == dataset_name, ]
    if (nrow(row) == 1) return(row$description)
    return("No description available.")
  }
  
  
  # Table of gene names to allow searching by AGI or gene symbol
  gene_symbol_map <- read_csv("gene_symbols.csv") |>
    mutate(display = if_else(GeneSymbol != "", paste0(AGI, " (", GeneSymbol, ")"), AGI)) |>
    distinct()
  
  
  # Load named gene sets (available only for Heatmap tab)
  gene_set_files <- list.files("gene_sets", pattern = "\\.txt$", full.names = TRUE)
  gene_sets <- lapply(gene_set_files, function(f) trimws(read_lines(f)))
  names(gene_sets) <- tools::file_path_sans_ext(basename(gene_set_files))
  
  
  # Load replicate-level data for Expression tab
  replicate_data <- reactive({
    req(input$experiment_line)
    path <- experiment_dirs[[input$experiment_line]]
    
    full_expr <- read_table(file.path(path, "tpm.tpm"))
    meta <- read_table(file.path(path, "meta.txt")) |>
      mutate(Replicate = as.character(Replicate))
    
    validate(
      need(all(c("Group", "Genotype", "Treatment", "Organ") %in% colnames(meta)),
           "Metadata must contain Group, Genotype, Treatment, Organ columns.")
    )
    
    list(full_expr = full_expr, metadata = meta)
  })
  
  
  # Load averaged data for Heatmap tab
  avg_data <- reactive({
    req(input$experiment_heatmap)
    path <- experiment_dirs[[input$experiment_heatmap]]
    avg_data <- read_table(file.path(path, "tpm_avg.txt"))
  })
  
  # Load replicate data for PCA tab (reuse replicate_data reactive but separate for clarity)
  replicate_data_pca <- reactive({
    req(input$experiment_pca)
    path <- experiment_dirs[[input$experiment_pca]]
    
    full_expr <- read_table(file.path(path, "tpm.tpm"))
    meta <- read_table(file.path(path, "meta.txt")) |>
      mutate(Replicate = as.character(Replicate))
    
    list(full_expr = full_expr, metadata = meta)
  })
  
  # Merge replicate expression with metadata for Expression tab
  full_data <- reactive({
    dat <- replicate_data()
    df <- dat$full_expr |> 
      pivot_longer(-Gene, names_to = "Replicate", values_to = "TPM") |> 
      mutate(Replicate = as.character(Replicate)) |> 
      left_join(dat$metadata, by = "Replicate")
    
    validate(need(all(c("Group", "Genotype", "Treatment", "Organ") %in% colnames(df)),
                  "Metadata must contain Group, Genotype, Treatment, Organ columns."))
    
    # Define custom Group order per dataset
    group_order <- switch(input$experiment_line,
                          "LB + LRFR - Autophagy" = c("Cot.Col0.WL",
                                                      "Cot.Col0.LB",
                                                      "Cot.Col0.LRFR",
                                                      "Cot.pif457.WL",
                                                      "Cot.pif457.LB",
                                                      "Cot.pif457.LRFR",
                                                      "Cot.yuc2589.WL",
                                                      "Cot.yuc2589.LB",
                                                      "Cot.yuc2589.LRFR",
                                                      "Cot.smt2.WL",
                                                      "Cot.smt2.LB",
                                                      "Cot.smt2.LRFR",
                                                      "Hyp.Col0.WL",
                                                      "Hyp.Col0.LB",
                                                      "Hyp.Col0.LRFR",
                                                      "Hyp.pif457.WL",
                                                      "Hyp.pif457.LB",
                                                      "Hyp.pif457.LRFR",
                                                      "Hyp.yuc2589.WL",
                                                      "Hyp.yuc2589.LB",
                                                      "Hyp.yuc2589.LRFR",
                                                      "Hyp.smt2.WL",
                                                      "Hyp.smt2.LB",
                                                      "Hyp.smt2.LRFR",
                                                      "WL", "LB", "LRFR"),
                          "LRFR - PIF7 KI lines" = c("HRFR", "LRFR", "28C", "LRFR_28C"),
                          "LRFR - time course (15-90 min)" = c("c0_w",
                                                               "c15_w",
                                                               "c45_w",
                                                               "c90_w",
                                                               "c180_w",
                                                               "c15_fr",
                                                               "c45_fr",
                                                               "c90_fr",
                                                               "c180_fr",
                                                               "h0_w",
                                                               "h15_w",
                                                               "h45_w",
                                                               "h90_w",
                                                               "h180_w",
                                                               "h15_fr",
                                                               "h45_fr",
                                                               "h90_fr",
                                                               "h180_fr",
                                                               "0_min", "15_min", "45_min", "90_min", "180_min"),
                          "LRFR - wounding" = c("Col-0_HRFR_1.5h_nw",
                                                "Col-0_HRFR_1.5h_w",
                                                "Col-0_LRFR_1.5h_nw",
                                                "Col-0_LRFR_1.5h_w",
                                                "Col-0_HRFR_4h_nw",
                                                "Col-0_HRFR_4h_w",
                                                "Col-0_LRFR_4h_nw",
                                                "Col-0_LRFR_4h_w",
                                                "aos_HRFR_1.5h_nw",
                                                "aos_HRFR_1.5h_w",
                                                "aos_LRFR_1.5h_nw",
                                                "aos_LRFR_1.5h_w",
                                                "aos_HRFR_4h_nw",
                                                "aos_HRFR_4h_w",
                                                "aos_LRFR_4h_nw",
                                                "aos_LRFR_4h_w",
                                                "HRFR_1.5h", "HRFR_4h", "LRFR_1.5h", "LRFR_4h"),
                          unique(df$Group))
    
    df$Group <- factor(df$Group, levels = group_order)
    
    df
    
    #Define custom order of Treatment
    treatment_order <- switch(input$experiment_line,
                          "LB + LRFR - Autophagy" = c("WL", "LB", "LRFR"),
                          "LRFR - PIF7 KI lines" = c("HRFR", "LRFR", "28C", "LRFR_28C"),
                          "LRFR - time course (15-90 min)" = c("0_min", "15_min", "45_min", "90_min", "180_min"),
                          "LRFR - wounding" = c("HRFR_1.5h", "HRFR_4h", "LRFR_1.5h", "LRFR_4h"),
                          "WL + R + RB - Leaf" = c("W_1", "W_7", "W_25", "W_73",
                                                   "R_1", "R_7", "R_25", "R_73",
                                                   "RB_1", "RB_7", "RB_25", "RB_73"),
                          unique(df$Treatment))
    
    df$Treatment <- factor(df$Treatment, levels = treatment_order)
    df
  })
  
  
  # Averaged TPM for Heatmap tab
  expr_matrix <- reactive({
    mat <- avg_data()
    meta <- read_table(file.path(experiment_dirs[[input$experiment_heatmap]], "meta.txt")) |>
      dplyr::select(-Replicate) |>
      distinct(Group, Genotype, Treatment, Organ)
    
    #join with metadata
    mat_long <- mat |>
      pivot_longer(-Gene, names_to = "Group", values_to = "TPM")
    
    mat_joined <- left_join(mat_long, meta, by = "Group")
    
    if (!is.null(input$genotype_heatmap)) {
      mat_joined <- mat_joined |> filter(Genotype %in% input$genotype_heatmap)
    }
    
    mat_filtered <- mat_joined |>
      dplyr::select(Gene, Group, TPM) |>
      pivot_wider(names_from = Group, values_from = TPM) |>
      column_to_rownames("Gene") |>
      as.matrix()
    
    mat_filtered[is.na(mat_filtered)] <- 0
    mat_filtered
  })
  
  
  # Z-score for Heatmap tab
  zscore_matrix <- reactive({
    mat <- expr_matrix()
    z_mat <- t(scale(t(mat)))
    z_mat[is.na(z_mat)] <- 0
    z_mat
  })
  
  
  # Button selection for Genotype in Heatmap tab
  output$genotype_buttons <- renderUI({
    req(full_data())
    choices <- sort(unique(full_data()$Genotype))
    checkboxGroupInput("genotype", "Select genotype(s):", choices = choices, selected = choices)
  })
  
  
  # Heatmap tab - output
  output$genotype_selector_heatmap <- renderUI({
    req(input$experiment_heatmap)
    path <- experiment_dirs[[input$experiment_heatmap]]
    meta <- read_table(file.path(path, "meta.txt")) |>
      mutate(Replicate = as.character(Replicate))
    
    validate(need("Genotype" %in% colnames(meta), "Metadata must include 'Genotype' column."))
    
    genotypes <- sort(unique(meta$Genotype))
    checkboxGroupInput("genotype_heatmap", "Select genotype(s):", choices = genotypes, selected = genotypes)
  })
  
  
  # PCA tab - output
  # Reactive PCA tab data filtered by selected Treatment(s)
  pca_data <- reactive({
    req(input$experiment_pca, input$treatments_pca)
    
    # Load data
    dat <- replicate_data_pca()
    expr_df <- dat$full_expr
    meta_df <- dat$metadata |> 
      filter(Treatment %in% input$treatments_pca)
    
    # Match expression data with selected replicates
    valid_reps <- meta_df$Replicate
    expr_mat <- expr_df |> 
      dplyr::select(Gene, all_of(valid_reps)) |> 
      column_to_rownames("Gene") |> 
      as.matrix()
    
    # Log transform (asinh for negatives, log2 otherwise)
    has_negative <- any(expr_mat < 0, na.rm = TRUE)
    transformed <- if (has_negative) {
      t(asinh(expr_mat))
    } else {
      t(log2(expr_mat + 1))
    }
    
    # Remove NAs
    transformed[is.na(transformed)] <- 0
    
    # Remove zero-variance genes
    gene_vars <- matrixStats::rowVars(t(transformed))
    nonzero_var_idx <- which(gene_vars > 0)
    
    if (length(nonzero_var_idx) == 0) {
      # Abort PCA if no variable genes
      return(data.frame())
    }
    
    # Select top 500 most variable genes (from non-zero variance only)
    top_genes <- order(gene_vars[nonzero_var_idx], decreasing = TRUE)[1:min(500, length(nonzero_var_idx))]
    selected_genes <- nonzero_var_idx[top_genes]
    filtered <- transformed[, selected_genes, drop = FALSE]
    
    # PCA
    pca_res <- prcomp(filtered, scale. = TRUE)
    
    # Format PCA output
    df_out <- data.frame(
      PC1 = pca_res$x[, 1],
      PC2 = pca_res$x[, 2],
      Replicate = rownames(pca_res$x)
    ) |> 
      left_join(meta_df, by = "Replicate") |>
      filter(Genotype %in% input$genotype_pca,
             Treatment %in% input$treatments_pca)
    
    df_out
  })
    
  
  # Update selectors after PCA loads
  output$genotype_pca <- renderUI({
    req(input$experiment_pca)
    path <- experiment_dirs[[input$experiment_pca]]
    meta <- read_table(file.path(path, "meta.txt")) |>
      mutate(Replicate = as.character(Replicate))
    
    validate(need("Genotype" %in% colnames(meta), "Metadata must include 'Genotype' column."))
    
    genotypes <- sort(unique(meta$Genotype))
    checkboxGroupInput("genotype_pca", "Select genotype(s):", choices = genotypes, selected = genotypes)
  })
  
  output$treatments_pca <- renderUI({
    req(input$experiment_pca)
    path <- experiment_dirs[[input$experiment_pca]]
    meta <- read_table(file.path(path, "meta.txt")) |>
      mutate(Replicate = as.character(Replicate))
    
    validate(need("Treatment" %in% colnames(meta), "Metadata must include 'Treatment' column."))
    
    treatment <- sort(unique(meta$Treatment))
    checkboxGroupInput("treatments_pca", "Select treatment(s):", choices = treatment, selected = treatment)
  })
  
  
  # Option to select Genes in Heatmap tab
  observeEvent(full_data(), {
    valid_genes <- unique(full_data()$Gene)
    gene_choices <- gene_symbol_map |>
      filter(AGI %in% valid_genes) |>
      mutate(label = ifelse(GeneSymbol != "", paste0(AGI, " (", GeneSymbol, ")"), AGI)) |>
      dplyr::select(label, AGI) |>
      dplyr::distinct()
    
    choices_vector <- setNames(gene_choices$AGI, gene_choices$label)
    
    updateSelectizeInput(session, "genes", choices = choices_vector, server = TRUE)
  })
  
  
  # Function to join AGI gene denominators and gene names, to select Gene by AGI or SYMBOL
  observeEvent(avg_data(), {
    valid_genes <- trimws(avg_data()$Gene)
    gene_choices <- gene_symbol_map |>
      filter(AGI %in% valid_genes) |>
      mutate(label = ifelse(GeneSymbol != "", paste0(AGI, " (", GeneSymbol, ")"), AGI)) |>
      dplyr::select(label, AGI) |>
      distinct()
    
    choices_vector <- setNames(gene_choices$AGI, gene_choices$label)
    
    updateSelectizeInput(session, "genes_heat", choices = choices_vector, server = TRUE)
  })
  
  
  # Option to select Gene list (Heatmap tab)
  observeEvent(input$experiment_heatmap, {
    updateSelectInput(session, "gene_set_choice",
                      choices = names(gene_sets),
                      selected = names(gene_sets)[1])
  })
  
  
  # Message to select at least two genes for Heatmap tab
  output$message <- renderText({
    if (is.null(input$genes) || length(input$genes) == 0) {
      "Select at least one gene to view plots."
    } else NULL
  })
  
  
  # Information about datasets, displayed above Expression tab - output
  output$dataset_info_line <- renderUI({
    req(input$experiment_line)
    desc <- get_description(input$experiment_line)
    HTML(paste0("<div style='margin-bottom: 20px;'>", commonmark::markdown_html(desc), "</div>"))
  })
  
  
  # Information about datasets, displayed above Heatmap tab - output
  output$dataset_info_heatmap <- renderUI({
    req(input$experiment_heatmap)
    desc <- get_description(input$experiment_heatmap)
    HTML(paste0("<div style='margin-bottom: 20px;'>", commonmark::markdown_html(desc), "</div>"))
  })
  
  
  # Render Line Plot in Expression tab
  output$line_plot <- renderPlot({
    req(input$genes, input$genotype)
    df <- full_data() |> 
      filter(Gene %in% input$genes, Genotype %in% input$genotype) |>
      left_join(gene_symbol_map, by = c("Gene" = "AGI"))
    
    # Base ggplot
    p <- ggplot(df, aes(x = Treatment, y = TPM, color = Genotype)) +
      geom_line(aes(group = Genotype),
                stat = "summary", fun = "mean",
                linetype = 2, size = 0.75) +
      geom_point(aes(fill = Genotype),
                 size = 4, shape = 21, color = "grey10") +
      theme_pubclean() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            strip.background = element_blank(),
            strip.text = element_text(size = 12, face = "bold")) +
      labs(title = "Gene Expression (TPM/CPM)", x = "Group", y = "TPM")
    
    # Conditional facetting logic
    if (input$experiment_line == "WL + R + RB - Leaf") {
      # This dataset uses a special faceting rule
      p <- p + facet_grid(display ~ factor(Treatment_group, c("W", "R", "RB")),
                          scales = "free"
                          #, ncol = 3
                          )
    } else {
      # Default faceting
      p <- p + facet_wrap(display ~ Organ, scales = "free", ncol = 2)
    }
    
    p
  },
  height = function() {
    base_height <- 400
    per_row <- 400
    
    n_genes <- length(input$genes)
    n_organs <- length(unique(full_data()$Organ))
    n_facets <- n_genes * n_organs
    ncol <- 2
    n_rows <- ceiling(n_facets / ncol)
    
    base_height + per_row * max(0, n_rows - 1)
  })
  
  
  # Render Box Plot in Expression tab
  output$box_plot <- renderPlot({
    req(input$genes, input$genotype)
    df <- full_data() |> 
      filter(Gene %in% input$genes, Genotype %in% input$genotype) |>
      left_join(gene_symbol_map, by = c("Gene" = "AGI"))
    
    #ggplot line plot
    ggplot(df, aes(x = Treatment, y = TPM, fill = Genotype)) +
      geom_boxplot(position = position_dodge(0.7), 
               width = 0.5, alpha = 0.5) +
      geom_point(aes(color = Genotype),
                 position = position_dodge(0.7), size = 3, shape = 21, color = "grey10") +
      facet_wrap(display ~ Organ, scales = "free", ncol = 2) +
      theme_pubclean() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            strip.background = element_blank(),
            strip.text = element_text(size = 12, face = "bold")) +
      labs(title = "Gene Expression (TPM/CPM)", x = "Group", y = "TPM")
  }, 
  
  height = function() {
    base_height <- 400
    per_row <- 400
    
    # Approximate the number of facets (Gene x Organ combos)
    n_genes <- length(input$genes)
    n_organs <- length(unique(full_data()$Organ))
    n_facets <- n_genes * n_organs
    
    ncol <- 2
    n_rows <- ceiling(n_facets / ncol)
    
    base_height + per_row * max(0, n_rows - 1)
  })
  
  
  # Render Data Table of gene expression in Expression tab
  output$line_data_table <- renderDT({
    req(input$genes, input$genotype)
    
    path <- experiment_dirs[[input$experiment_line]]
    meta <- read_table(file.path(path, "meta.txt")) |>
      dplyr::select(Group, Genotype, Treatment, Organ) |>
      dplyr::distinct(Group, Genotype, Treatment, Organ) |>
      mutate(Group = trimws(Group))
    avg_data <- read_table(file.path(path, "tpm_avg.txt")) |>
      pivot_longer(-Gene, names_to = "Group", values_to = "Expression (TPM/CPM)") |>
      mutate(Group = trimws(Group)) |>
      left_join(meta, by = "Group") |>
      left_join(gene_symbol_map, by = c("Gene" = "AGI")) |>
      dplyr::select(Gene, GeneSymbol, Genotype, Organ, Treatment, "Expression (TPM/CPM)") |>
      filter(Gene %in% input$genes, Genotype %in% input$genotype) 
    
    datatable(avg_data, options = list(pageLength = 10))
    
  })
  
  
  # Download Expression table of selected genes
  output$download_data <- downloadHandler(
    filename = function() {
      paste0("expression_data_", Sys.Date(), ".txt")
    },
    content = function(file) {
      req(input$genes, input$genotype)
      path <- experiment_dirs[[input$experiment_line]]
      meta <- read_table(file.path(path, "meta.txt")) |>
        dplyr::select(Group, Genotype, Treatment, Organ)
      df <- avg_data() |>
        pivot_longer(-Gene, names_to = "Group", values_to = "Expression (TPM/CPM)") |>
        left_join(meta, by = "Group") |>
        dplyr::select(Gene, Genotype, Organ, Treatment, "Expression (TPM/CPM)") |>
        filter(Gene %in% input$genes, Genotype %in% input$genotype)
      write.table(df, file, row.names = FALSE, quote = FALSE, sep = "\t")
    }
  )
  


  # Render Heatmap in Heatmap tab
  output$heatmap_plot <- renderPlot({
    mat <- if (input$heatmap_scale == "Z-score") {
      zscore_matrix()
    } else {
      expr_matrix()
    }
    
    rownames(mat) <- trimws(rownames(mat))
    
    # Select genes
    selected_genes <- if (input$gene_selection_mode == "Select manually") {
      trimws(input$genes_heat)
    } else {
      req(input$gene_set_choice)
      gene_sets[[input$gene_set_choice]]
    }
    
    selected_genes <- intersect(selected_genes, rownames(mat))
    mat_sub <- mat[selected_genes, , drop = FALSE]
    
    
    # Check if there are enough genes
    validate(
      need(nrow(mat_sub) >= 2, "Select at least two genes")
    )
    
    if (ncol(mat_sub) < 1) {
      plot.new()
      text(0.5, 0.5, "No columns to display")
      return()
    }
    
    # Load metadata
    meta <- read_table(file.path(experiment_dirs[[input$experiment_heatmap]], "meta.txt")) |>
      mutate(Group = trimws(Group))
    
    # Match columns in mat_sub with metadata
    meta <- meta[match(colnames(mat_sub), meta$Group), ]
    
    # Create annotation for Organ and Genotype
    annotation_col <- data.frame(
      Organ = meta$Organ,
      Genotype = meta$Genotype
    )
    rownames(annotation_col) <- meta$Group
    
    # Order columns by Organ, then Genotype
    col_order <- order(annotation_col$Organ, annotation_col$Genotype)
    mat_sub <- mat_sub[, col_order, drop = FALSE]
    annotation_col <- annotation_col[col_order, , drop = FALSE]
    
    # Compute gaps between Organ groups first, then Genotype within Organ
    split_factor <- paste(annotation_col$Organ, annotation_col$Genotype, sep = "_")
    gap_positions <- which(diff(as.numeric(factor(split_factor))) != 0)
    
    # Draw heatmap with annotation bar
    pheatmap(mat_sub,
             cluster_rows = TRUE,
             cluster_cols = FALSE,
             gaps_col = gap_positions,
             angle_col = 45,
             annotation_col = annotation_col,
             show_colnames = TRUE,
             show_rownames = TRUE)
  })
 

  # Render PCA in PCA tab
  output$pca_plot <- renderPlot({
    df <- pca_data()
    req(nrow(df) > 0)
    
    df <- df |>
      filter(Genotype %in% input$genotype_pca,
             Treatment %in% input$treatments_pca)
    
    # Count number of unique groups
    group_count <- df |> 
      dplyr::count(Group) |> 
      filter(n >= 3) |> 
      pull(Group)
    
    df_labels <- df |> group_by(Group) |> slice(1) |> ungroup()
    df_labels_organ <- df |> group_by(Organ) |> slice(1) |> ungroup()
    
    # Start base plot
    p <- ggplot(df, aes(x = PC1, y = PC2, shape = Genotype, color = Treatment)) +
      geom_point(size = 5, alpha = 0.75) +
      theme_pubclean() +
      labs(title = paste("PCA plot of", input$experiment_pca),
           x = "PC1", y = "PC2") +
      theme(
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11)
      )
    
    
    # Add encircling if more than 1 group
    p <- p + geom_mark_ellipse(
      data = filter(df, Group %in% group_count),
      aes(group = Group, color = Treatment),
      linetype = "dashed", linewidth = 0.3,
      expand = unit(2, "mm"),
      inherit.aes = TRUE)
    

        # Add labels if checkbox is selected
        if (isTRUE(input$show_labels)) {
          p <- p +
            geom_label_repel(data = df_labels, 
                             mapping = aes(x = PC1, y = PC2, label = Group), 
                             fill = "white", alpha = 0.85, 
                             size = 5, force = 2, box.padding = 1) +
            geom_text_repel(data = df_labels_organ, 
                            mapping = aes(x = PC1, y = PC2, label = Organ), 
                            color = "grey10", alpha = 0.25, 
                            size = 10, box.padding = 3)
        }
    
    p
    
  })

}

shinyApp(ui, server)

