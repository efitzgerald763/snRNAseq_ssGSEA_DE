library(Seurat)
library(tidyverse)
library(harmony)
library(scDblFinder)
library(Azimuth)
library(SeuratData)


# Function to identify outliers based on median absolute deviation (MAD)
mad_outlier <- function(sobj, metric, nmads) {
  M <- sobj@meta.data[[metric]]
  median_M <- median(M, na.rm = TRUE)
  mad_M <- mad(M, na.rm = TRUE)
  outlier <- (M < (median_M - nmads * mad_M)) | (M > (median_M + nmads * mad_M))
  return(outlier)
}

# Load DLPFC dataset
DLPFC <- Read10X(data.dir = "/Single_cell_studies/Female_MDD/Barcodes/", gene.column = 1)
DLPFC <- CreateSeuratObject(counts = DLPFC, project = "DLPFC")

# Process metadata
metadata <- DLPFC@meta.data
metadata$cells <- rownames(metadata)
metadata$Sample <- word(metadata$cells, 1, sep = "\\.")

# Merge with external metadata
MET <- openxlsx::read.xlsx("Single_cell_studies/Female_MDD/Metadata.xlsx") %>%
  dplyr::select(Sample, Condition, Bank, Batch, Chemistry, Sequencing, Age, Race, Sex)

met_combined <- merge(metadata, MET, by = 'Sample')
rownames(met_combined) <- met_combined$cells

# Change names of existing cell types
met_combined$celltype_study <- case_when(
  grepl("Mic", rownames(met_combined), ignore.case = TRUE) ~ "Microlia",
  grepl("Oli", rownames(met_combined), ignore.case = TRUE) ~ "Oligodendrocytes",
  grepl("InN", rownames(met_combined), ignore.case = TRUE) ~ "Inhibitory_neurons",
  grepl("Ast", rownames(met_combined), ignore.case = TRUE) ~ "Astrocytes",
  grepl("OPC", rownames(met_combined), ignore.case = TRUE) ~ "OPCs",
  grepl("ExN", rownames(met_combined), ignore.case = TRUE) ~ "Excitatory_neurons",
  grepl("End", rownames(met_combined), ignore.case = TRUE) ~ "Endothelial",
  TRUE ~ NA_character_
)

DLPFC@meta.data <- met_combined

# Subset the data based on specific criteria
DLPFC <- DLPFC %>%
  subset(subset = Sex == 'Female') %>%
  subset(subset = Chemistry == 'v3')

# Calculate ribosomal and mitochondrial ratios
DLPFC$riboRatio <- PercentageFeatureSet(object = DLPFC, pattern = "^RP[LS]")
DLPFC$mitoRatio <- PercentageFeatureSet(object = DLPFC, pattern = "^MT")

# Filter outliers
Keep_me <- !mad_outlier(DLPFC, 'nFeature_RNA', 5) &
  !mad_outlier(DLPFC, 'nCount_RNA', 5) &
  !mad_outlier(DLPFC, 'riboRatio', 3) &
  !mad_outlier(DLPFC, 'mitoRatio', 3)

DLPFC <- subset(DLPFC, cells = which(Keep_me))

# Data normalization and dim reduction
DLPFC_harmony <- DLPFC %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(vars.to.regress = c('riboRatio', 'mitoRatio')) %>%
  RunPCA(npcs = 30)

# Run Harmony integration and clustering
DLPFC_harmony <- DLPFC_harmony %>%
  RunHarmony("Batch", plot_convergence = TRUE) %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = c(0.01, 0.1))

# Doublet detection using scDblFinder
Doublet_obj <- as.SingleCellExperiment(DLPFC_harmony)
Doublet_obj <- scDblFinder(Doublet_obj, samples = "Sample", clusters = "seurat_clusters")
Doublet_obj <- as.Seurat(Doublet_obj)

# Visualize doublet scores
FeaturePlot(object = Doublet_obj, features = 'scDblFinder.score', raster = F)
DLPFC_harmony <- subset(Doublet_obj, subset = scDblFinder.class == 'singlet')

# Re-run UMAP and clustering after doublet removal
DLPFC_harmony <- DLPFC_harmony %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = c(0.01, 0.1))

# Plot UMAP with different groupings
DimPlot(DLPFC_harmony, reduction = "umap", group.by = "Condition", raster = F)
DimPlot(DLPFC_harmony, reduction = "umap", group.by = "Sample", raster = F)
DimPlot(DLPFC_harmony, reduction = "umap", label = TRUE, repel = TRUE, raster = F)

# Stacked bar plot of cell proportions by cluster
meta <- DLPFC_harmony@meta.data
ggplot(meta, aes(fill = seurat_clusters, x = Sample)) +
  geom_bar(position = "fill", stat = "count") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 17),
        legend.title = element_blank()) +
  xlab("") +
  ylab("Proportion of cells")


# Set default assay and identify clusters
DefaultAssay(DLPFC_harmony) <- 'RNA'
Idents(DLPFC_harmony) <- "RNA_snn_res.0.01"

# Plot violin plots for selected markers
VlnPlot(DLPFC_harmony, 
        features = c('SLC1A2', 'SLC1A3', 'CUX2', 'RORB', 'SLC17A7', 
                     'GAD1', 'GAD2', 'DOCK8', 'VCAN', 'PDGFRA', 
                     'FLT1', 'ABCB1', 'MOG', 'MBP'), 
        stack = TRUE, fill.by = 'ident') +
  xlab("") + ylab('') +
  theme_classic() +
  geom_boxplot(width = 0.05, show.legend = FALSE, position = position_dodge(0.9)) +
  theme(legend.position = 'none',
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 16),
        strip.background = element_blank())

# Identify all cluster markers
markers <- FindAllMarkers(DLPFC_harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10_markers <- markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)


# Azimuth annotation
DLPFC_harmony <- RunAzimuth(DLPFC_harmony, reference = "humancortexref")

# Plot UMAP with predicted subtypes
Idents(DLPFC_harmony) <- "predicted.subclass"
DimPlot(DLPFC_harmony, reduction = "umap", label = TRUE, repel = TRUE, raster = T) +
  theme(legend.position = 'none')

# Plot UMAP with study cell types
Idents(DLPFC_harmony) <- "celltype_study"
DimPlot(DLPFC_harmony, reduction = "umap", label = TRUE, repel = TRUE, raster = T) +
  theme(legend.position = 'bottom')

# Rename clusters with more informative labels
DLPFC_harmony <- RenameIdents(object = DLPFC_harmony, 
                              "Astro" = "Astrocyte",
                              "Endo" = "Stromal",
                              "L2/3 IT" = "Ex_neuron",
                              "L5 ET" = "Ex_neuron",
                              "L5 IT" = "Ex_neuron",
                              "L5/6 NP" = "Ex_neuron",
                              "L6 CT" = "Ex_neuron",
                              "L6 IT" = "Ex_neuron",
                              "L6 IT Car3" = "Ex_neuron",
                              "L6b" = "Ex_neuron",
                              "Lamp5" = "In_neuron",
                              "Micro-PVM" = "Microglia",
                              "Oligo" = "Oligodendrocyte",
                              "OPC" = "OPC",
                              "Pvalb" = "In_neuron",
                              "Sncg" = "In_neuron",
                              "Sst" = "In_neuron",
                              "Sst Chodl" = "In_neuron",
                              "Vip" = "In_neuron",
                              "VLMC" = "Stromal")

# Save final objects
save(DLPFC_harmony, DLPFC, Doublet_obj, file = "/Single_cell_studies/Female_MDD/Female_depression_harmony.RData")

