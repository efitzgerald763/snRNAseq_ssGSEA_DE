library(GSVA)
library(msigdbr)
library(Seurat)

# load and pseudobulk snRNA-seq data
load("/Single_cell_studies/Female_MDD/Female_depression_harmony.RData")

# Use the aggregate expression function from seurat that presumes data is normalized
Pseudobulk <- AggregateExpression(DLPFC_harmony, 
                           group.by = c("celltype_study", "Sample"),
                           assays = 'RNA',
                           slot = "counts",
                           return.seurat = TRUE)

# Access the RNA assay counts slot
counts_data <- GetAssayData(Pseudobulk, slot = "counts", assay = "RNA")

# Split the counts data by celltype_study
celltype_list <- split(as.data.frame(counts_data), Pseudobulk@meta.data$celltype_study)

# Optional: Name each dataframe in the list by the corresponding celltype
names(celltype_list) <- unique(Pseudobulk@meta.data$celltype_study)

# Extract metadata
datMeta <- Pseudobulk@meta.data




------------------------------------------------------------------------------------
# Generate ssGSEA scores for each ontology term in each cluster

library(GSVA)
library(BiocParallel)

# Select human gene ontology sets
go_gene_sets <- msigdbr(species = "Homo sapiens", subcategory = "CP:REACTOME")
go_gene_sets <- split(go_gene_sets$gene_symbol, go_gene_sets$gs_name)
x <- msigdbr_collections()

# use multiple cores
bpparam <- MulticoreParam(workers = 10)

# Function to compute ssGSEA scores
compute_ssgsea <- function(df) {
  gsva(df, go_gene_sets,
       method = "ssgsea",
       min.sz = 10,
       max.sz = 500,
       BPPARAM = bpparam)
}

# Apply the function to the count list
ssgsea_scores_list <- lapply(celltype_list, compute_ssgsea)

# Stop the parallel cluster
stopCluster(bpparam)

------------------------------------------------------------------------------------
# Check differential expression at the pathway level using a linear model

  datMeta <- column_to_rownames(datMeta, var = "Sample")

# Define function to use limma to run a linear regression for ssGSEA scores in cntrl and MDD
run_analysis <- function(data) {
  # Design matrix 
  design <- model.matrix(~ Condition, data = datMeta)
  
  # Fit model
  fit <- lmFit(data, design)
  fit <- eBayes(fit)
  
  # Adjust p-values
  topTable <- topTable(fit, adjust="fdr")
  
  return(topTable)
}

# Apply the function to each dataframe in the list
results_list <- lapply(celltype_list, run_analysis)


