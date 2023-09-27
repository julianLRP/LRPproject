library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)

# Load your Seurat object
processed_LRP <- readRDS('LRP.processed2.rds')

# Add a new metadata column (if not already added)
processed_LRP@meta.data$earlymappingbycelltype <- ""

processed_LRP@meta.data$earlymappingbycelltype[is.na(processed_LRP@meta.data$earlymappingbycelltype)] <- "Unassigned"

#add a new metadata column for stages
processed_LRP@meta.data$LRP_stage<- ""

# Set the identity of the Seurat object based on your new column
processed_LRP <- SetIdent(processed_LRP, value="earlymappingbycelltype")

# Set the default assay to RNA
DefaultAssay(processed_LRP) <- 'RNA'

# Define a vector with all clusters related to Pericycle
pericycle_related_clusters <- c("0", "1", "5", "8", "9", "3", "4", "15", "24")

# Subset the Seurat object for Pericycle related clusters
Pericycle_subset <- subset(x = processed_LRP, subset = seurat_clusters %in% pericycle_related_clusters)

Pericycle_cell_ids<-WhichCells(Pericycle_subset)

processed_LRP@meta.data[Pericycle_cell_ids, "earlymappingbycelltype"]<- "Pericycle"

# Define a vector with all clusters related to Endodermis
endodermis_related_clusters <- c("16")

# Subset the Seurat object for Endodermis related clusters
Endodermis_subset <- subset(x = processed_LRP, subset = seurat_clusters %in% endodermis_related_clusters)

# Get cell IDs for Endodermis cells
Endodermis_cell_ids <- WhichCells(Endodermis_subset)

# Label these cells as Endodermis in the metadata
processed_LRP@meta.data[Endodermis_cell_ids, "earlymappingbycelltype"] <- "Endodermis"


protoxylem_related_clusters <- c()

# Define a vector with all clusters related to Cortex
cortex_related_clusters <- c("28", "20")

# Subset the Seurat object for Cortex related clusters
Cortex_subset <- subset(x = processed_LRP, subset = seurat_clusters %in% cortex_related_clusters)

# Get cell IDs for Cortex cells
Cortex_cell_ids <- WhichCells(Cortex_subset)

# Label these cells as Cortex in the metadata
processed_LRP@meta.data[Cortex_cell_ids, "earlymappingbycelltype"] <- "Cortex"

# Define a vector with clusters related to Outer Layer of LRP
outer_LRP_clusters <- c("23", "29","12","7","30")

# Subset the Seurat object for these clusters
Outer_LRP_subset <- subset(x = processed_LRP, subset = seurat_clusters %in% outer_LRP_clusters)

# Get cell IDs for these clusters
Outer_LRP_cell_ids <- WhichCells(Outer_LRP_subset)

# Label these cells as Epidermis in the metadata
processed_LRP@meta.data[Outer_LRP_cell_ids, "earlymappingbycelltype"] <- "Epidermis"

#Subset the Seurat object based on specific clusters that were designated to likely be LRC and Columella
LRC <- subset(x = processed_LRP, subset = seurat_clusters %in% c("31", "10", "22", "11"))

lateral_root_cap_ids <- WhichCells (LRC)
##temporarily disabled.  
#Identify cells within this subset that meet the gene expression criteria
#columella_cell_ids <- WhichCells(LRC, expression = AT2G01420 > 0.1)

#Identify cells within this subset that do NOT meet the gene expression criteria
#lateral_root_cap_ids <- setdiff(rownames(LRC@meta.data), columella_cell_ids)

#Modify the 'earlymappingbycelltype' column in the metadata for these cells
#processed_LRP@meta.data[columella_cell_ids, "earlymappingbycelltype"] <- "Columella"
processed_LRP@meta.data[lateral_root_cap_ids, "earlymappingbycelltype"] <- "Root cap"

#subset clusters that contain Gata23 expression, they belong to XPP of early stages
Xylem_Pole_Pericycle_stage1_2_3_subset <- subset(x = processed_LRP, subset = seurat_clusters %in% c("5", "13", "4" ,"18"))

XPP_Cluster19 <-subset(x= processed_LRP, subset=seurat_clusters %in% c("19"))

XPP_Cluster19_ids <- WhichCells(XPP_Cluster19, expression = AT5G26930 > 0.1)

#Whichcells on the subset of clusters that belong to Xylem_Pole_Pericycle_stage1_2_3
Xylem_Pole_Pericycle_stage1_2_3_cell_ids<- WhichCells(Xylem_Pole_Pericycle_stage1_2_3_subset)

#Modify the 'earlymappingbycelltype' column in the metadata for these cells for QC
processed_LRP@meta.data[c(Xylem_Pole_Pericycle_stage1_2_3_cell_ids, XPP_Cluster19_ids), "earlymappingbycelltype"] <- "earlystage_primordia"

#subset clusters that belong to PIN1 AT1G73590 this will be the derivatives of the inner cell layers of stage 2 and beyond.
Inner_cell_derivatives_allLRP_stages<- subset(x = processed_LRP, subset = seurat_clusters %in% c("24", "17" ,"2" , "27" ,"18", "25" , "1"))

#Whichcells on the subset of clusters that belong to inner cell derivatives of all stages
Inner_cell_derivatives_allLRP_stages_cell_ids<- WhichCells(Inner_cell_derivatives_allLRP_stages)

#Modify the 'earlymappingbycelltype' column in the metadata for these cells for QC
processed_LRP@meta.data[Inner_cell_derivatives_allLRP_stages_cell_ids, "earlymappingbycelltype"] <- "Stele"

#Whichcells on WOX5
QC_cell_ids<- WhichCells(processed_LRP, expression = AT3G11260 > 0)

#Modify the 'earlymappingbycelltype' column in the metadata for these cells for QC
processed_LRP@meta.data[QC_cell_ids, "earlymappingbycelltype"] <- "QC"


#save the dataset
saveRDS(processed_LRP, 'LRP.processed2.rds')





