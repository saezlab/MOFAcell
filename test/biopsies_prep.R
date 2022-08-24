library(SeuratDisk)
library(HDF5Array)
library(scuttle)
library(tidyverse)
library(liana)
library(grDevices)
source("code/utils/ccc_utils.R")

h5 <- SeuratDisk::Connect("data/biopsies/seuratobject.h5Seurat")

h5[["meta.data"]]
print(h5$index())

# Read Sobj
sobj <- SeuratDisk::LoadH5Seurat("data/biopsies/seuratobject.h5Seurat",
                                 assays = c(RNA = "counts"),
                                 graphs = FALSE)
sce <- Seurat::as.SingleCellExperiment(sobj)
sce@assays@data@listData$logcounts <- NULL
rm(sobj)

# Basic Feature Filtering
sce <- sce[rowSums(counts(sce) >= 0) >= 5, ]

# Save as HDF5 Experiment:
HDF5Array::saveHDF5SummarizedExperiment(sce, "data/biopsies/hdf5arr/", replace=TRUE)

# Read HDF5
sce <- HDF5Array::loadHDF5SummarizedExperiment("data/biopsies/hdf5arr/")

# Normalize RNA assay
sce <- scuttle::logNormCounts(sce)
gc()

### Internal Functions for LIANA to filter cell types across samples (For SC)
# 1. Filter Samples - e.g. 3 z-scores < of SUM total counts
sce <- filter_samples(sce, sample_col = "Sample")


# 2. Filter Cell types by min.cell num + min.cell by sample
# filter
sce <- filter_nonabundant_celltypes(sce,
                                    sample_col = "Sample",
                                    idents_col = "predicted.annotation.l2",
                                    min_prop = 0.25,
                                    min_cells = 20,
                                    min_samples = 3
                                    )


# after filt
cairo_pdf(filename = "plots/biopsies_ctqc.pdf",
          height = 42,
          width = 18)
print(get_abundance_summary(sce,
                      sample_col = "Sample",
                      idents_col = "predicted.annotation.l2") %>%
    plot_abundance_summary(ncol=3))

dev.off()


# # --- Test Keep only IgA and TRL
# sce$Group %>% unique()
#
# groups_of_interest <- colData(sce) %>%
#     as_tibble(rownames="barcodes") %>%
#     # filter(Group %in% c("IgA", "CTRL", "MN")) %>%
#     pull("barcodes")
# sce <- sce[,groups_of_interest]


# Run LIANA by Sample ----
context_df_dict <- liana_bysample(sce = sce,
                                  sample_col = "Sample",
                                  condition_col = "Group",
                                  idents_col = "predicted.annotation.l2",
                                  permutation.params=list(nperms=2))
saveRDS(context_df_dict, "output/biopsies/context_df_dict.RDS")

