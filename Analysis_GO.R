# -------------------------------------------------------------------
# Gene Ontology Enrichment Analysis Pipeline
# -------------------------------------------------------------------

# 1. Safely verify and load required dependencies
required_pkgs <- c("clusterProfiler", "org.Hs.eg.db")
missing_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]

if (length(missing_pkgs) > 0) {
  stop(sprintf(
    "Missing required package(s): %s. Please install via Conda or BiocManager.", 
    paste(missing_pkgs, collapse = ", ")
  ))
}

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

# 2. Define the main analysis function
run_go_pipeline <- function(input_list_path = "hg38_GCGC_target_genes.tsv",
                            out_csv_path    = "hg38_GCGC_GO_results.csv",
                            out_plot_path   = "hg38_GCGC_GO_dotplot.pdf") {
  
  # Check for file existence
  if (!file.exists(input_list_path)) {
    stop(sprintf("Fatal Error: The target file '%s' could not be located.", input_list_path))
  }
  
  # Read target dataset
  raw_input <- read.table(input_list_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  target_symbols <- raw_input$V1
  message(sprintf("SUCCESS: Loaded %d gene symbols for downstream analysis.", length(target_symbols)))
  
  # Execute the GO mapping (Biological Process)
  message("STATUS: Initiating GO Enrichment Analysis (BP)...")
  go_enrichment_res <- enrichGO(
    gene          = target_symbols,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.01
  )
  
  # Parse output and generate deliverables
  if (is.null(go_enrichment_res) || nrow(as.data.frame(go_enrichment_res)) == 0) {
    warning("ALERT: No significant GO terms were isolated based on the defined cutoffs.")
    write.csv(data.frame(Analysis_Status = "No significant enrichment detected"), 
              out_csv_path, row.names = FALSE)
              
  } else {
    # Export tabular results
    write.csv(as.data.frame(go_enrichment_res), out_csv_path, row.names = FALSE)
    message(sprintf("SUCCESS: Enrichment data exported successfully to '%s'.", out_csv_path))
    
    # Render and save dotplot visualization
    pdf(out_plot_path, height = 8, width = 8)
    print(dotplot(go_enrichment_res, showCategory = 20, font.size = 6))
    dev.off()
    message(sprintf("SUCCESS: Visualization rendered and saved to '%s'.", out_plot_path))
  }
}

# -------------------------------------------------------------------
# 3. Execute the pipeline
# -------------------------------------------------------------------
run_go_pipeline()
