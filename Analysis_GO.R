# 1. Install required packages if missing
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
    BiocManager::install("org.Hs.eg.db")
if (!requireNamespace("clusterProfiler", quietly = TRUE))
    BiocManager::install("clusterProfiler")

# 2. Load the libraries
library(clusterProfiler)
library(org.Hs.eg.db)

# 3. Read the target gene list
genes <- read.table("hg38_GCGC_target_genes.tsv", 
                    sep = "\t", header = FALSE, stringsAsFactors = FALSE)$V1

# 4. Perform GO Enrichment Analysis (Biological Process)
ego <- enrichGO(gene          = genes,
                OrgDb         = org.Hs.eg.db,
                keyType       = "SYMBOL", 
                ont           = "BP",     
                pAdjustMethod = "BH",     
                qvalueCutoff  = 0.01)

# 5. Save results and generate plot
write.csv(as.data.frame(ego), "hg38_GCGC_GO_results.csv")

pdf("hg38_GCGC_GO_dotplot.pdf", height = 8, width = 8)
dotplot(ego, showCategory = 20, font.size = 6)
dev.off()
# ==========================================
# 1. Load the libraries (No auto-install blocks)
# ==========================================
if (!requireNamespace("clusterProfiler", quietly = TRUE) || 
    !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("Missing required R packages. Please install them via Conda first using:\n",
         "conda install -c bioconda bioconductor-clusterprofiler bioconductor-org.hs.eg.db")
}

library(clusterProfiler)
library(org.Hs.eg.db)

# ==========================================
# 2. Read the target gene list
# ==========================================
gene_file <- "hg38_GCGC_target_genes.tsv"

if (!file.exists(gene_file)) {
    stop(paste("Error: Input file", gene_file, "not found in the current directory!"))
}

genes <- read.table(gene_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)$V1
cat("Loaded", length(genes), "genes for analysis.\n")

# ==========================================
# 3. Perform GO Enrichment Analysis (BP)
# ==========================================
cat("Running GO Enrichment Analysis...\n")
ego <- enrichGO(gene          = genes,
                OrgDb         = org.Hs.eg.db,
                keyType       = "SYMBOL", 
                ont           = "BP",     
                pAdjustMethod = "BH",     
                qvalueCutoff  = 0.01)

# ==========================================
# 4. Save results and generate plot
# ==========================================
if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    warning("No significantly enriched GO terms found with the current thresholds.")
    write.csv(data.frame(Message="No significant results"), "hg38_GCGC_GO_results.csv")
} else {
    # Save the CSV
    write.csv(as.data.frame(ego), "hg38_GCGC_GO_results.csv")
    cat("Results saved to 'hg38_GCGC_GO_results.csv'.\n")
    
    # Generate the PDF Plot
    pdf("hg38_GCGC_GO_dotplot.pdf", height = 8, width = 8)
    print(dotplot(ego, showCategory = 20, font.size = 6))
    dev.off()
    cat("Plot saved to 'hg38_GCGC_GO_results.pdf'.\n")
}
# ==========================================
# 1. Load the libraries safely
# ==========================================
if (!requireNamespace("clusterProfiler", quietly = TRUE) || 
    !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("Missing required R packages. Please ensure your Conda environment is set up correctly.")
}

library(clusterProfiler)
library(org.Hs.eg.db)

# ==========================================
# 2. Read the target gene list
# ==========================================
gene_file <- "hg38_GCGC_target_genes.tsv"

if (!file.exists(gene_file)) {
    stop(paste("Error: Input file", gene_file, "not found in the current directory!"))
}

genes <- read.table(gene_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)$V1
cat("Loaded", length(genes), "genes for analysis.\n")

# ==========================================
# 3. Perform GO Enrichment Analysis (BP)
# ==========================================
cat("Running GO Enrichment Analysis...\n")
ego <- enrichGO(gene          = genes,
                OrgDb         = org.Hs.eg.db,
                keyType       = "SYMBOL", 
                ont           = "BP",     
                pAdjustMethod = "BH",     
                qvalueCutoff  = 0.01)

# ==========================================
# 4. Save results and generate plot
# ==========================================
if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    warning("No significantly enriched GO terms found with the current thresholds.")
    write.csv(data.frame(Message="No significant results"), "hg38_GCGC_GO_results.csv")
} else {
    # Save the CSV
    write.csv(as.data.frame(ego), "hg38_GCGC_GO_results.csv")
    cat("Results saved to 'hg38_GCGC_GO_results.csv'.\n")
    
    # Generate the PDF Plot
    pdf("hg38_GCGC_GO_dotplot.pdf", height = 8, width = 8)
    print(dotplot(ego, showCategory = 20, font.size = 6))
    dev.off()
    cat("Plot saved to 'hg38_GCGC_GO_dotplot.pdf'.\n")
}
