# Setup for reading files
here::i_am("Code/1 Read Data.R")
library(here)


################################################################################
# Cell Data
################################################################################
### Load Data ###
A3SS <- read.table(here("RNA Splicing Data", "Cell", "CT_A3SS_allEvents.txt"), header = T, row.names = 1)
A5SS <- read.table(here("RNA Splicing Data", "Cell", "CT_A5SS_allEvents.txt"),  header = T, row.names = 1)
MXE <- read.table(here("RNA Splicing Data", "Cell", "CT_MXE_allEvents.txt"),  header = T, row.names = 1)
RI <- read.table(here("RNA Splicing Data", "Cell", "CT_RI_allEvents.txt"),  header = T, row.names = 1)
SE <- read.table(here("RNA Splicing Data", "Cell", "CT_SE_allEvents.txt"),  header = T, row.names = 1)

### Add event column metadata
A3SS$Event <- "A3SS"
A5SS$Event <- "A5SS"
MXE$Event <- "MXE"
RI$Event <- "RI"
SE$Event <- "SE"

### Combine data
data_full <- rbind(A3SS,A5SS,RI,SE,MXE)
dim(data_full)

### Clean data of NA values
pNA <- colSums(is.na(data_full))/dim(data_full)[1]
head(pNA[order(pNA, decreasing = T)])

# Row 45 (X6105) is 85% NA - remove
# (~8000 obs if include row 45, ~100000 if not included)
data <- na.omit(data_full[-45])

data_mat <- data.matrix(data[,1:138])
cell_df <- t(data_mat)
dim(cell_df) # 138 108846

### Export & save dataframe
saveRDS(cell_df, here("RNA Splicing Data", "Cell Data.RDS"))

# Remove intermediate steps that are repeated for tissues, leave cell_df
rm(A3SS, A5SS, MXE, RI, SE, data, data_full, data_mat, pNA)
################################################################################
# Tissue Data
################################################################################
### Load Data ###
A3SS <- read.table(here("RNA Splicing Data", "Tissue", "A3SS_NoRDepth_cutOfff_allEvents.txt"), header = T, row.names = 1)
A5SS <- read.table(here("RNA Splicing Data", "Tissue", "A5SS_NoRDepth_cutOfff_allEvents.txt"),  header = T, row.names = 1)
MXE <- read.table(here("RNA Splicing Data", "Tissue", "MXE_NoRDepth_cutOfff_allEvents.txt"),  header = T, row.names = 1)
RI <- read.table(here("RNA Splicing Data", "Tissue", "RI_NoRDepth_cutOfff_allEvents.txt"),  header = T, row.names = 1)
SE <- read.table(here("RNA Splicing Data", "Tissue", "SE_NoRDepth_cutOfff_allEvents.txt"),  header = T, row.names = 1)

### Add event column metadata
A3SS$Event <- "A3SS"
A5SS$Event <- "A5SS"
MXE$Event <- "MXE"
RI$Event <- "RI"
SE$Event <- "SE"

data_full <- rbind(A3SS,A5SS,RI,SE,MXE)
dim(data_full)

# Check for missing columns
pNA <- colSums(is.na(data_full))/dim(data_full)[1]
head(pNA[order(pNA, decreasing = T)])
# No columns which are significantly NA

# Remove NA values
data <- na.omit(data_full)
dim(data) #[1] 126123     46

# Genes as columns & Tissues as rows
data_mat <- data.matrix(data[,1:45])
tissue_df <- t(data_mat)

# Export & save dataframe
saveRDS(tissue_df, here("RNA Splicing Data", "Tissue Data.RDS"))

# Remove intermediate steps leave tissue_df
rm(A3SS, A5SS, MXE, RI, SE, data, data_full, data_mat, pNA)

################################################################################
# Match with Biological source
################################################################################
# Source for cell data
bio_source <- read.csv(here("RNA Splicing Data", "RNASeq_Atlas_samples.csv"))
match <- paste0("X", bio_source$RNA_number)
Cell_source <- bio_source$Biological_source[match %in% rownames(cell_df)]
table(Cell_source) # 2 - 17 obs per system

# Source for tissue data
Tissue_source <-  bio_source$Biological_source[match %in% rownames(tissue_df)]
table(Tissue_source)

Cell_System_Types <- unique(Cell_source)
Tissue_System_Types <- unique(Tissue_source) # all unique

# Conversion from Tissue type to Cell system (to compare clustering)
conversion <- read.csv(here("RNA Splicing Data", "Tissue to Cell Sys.csv"))
map <- data.frame(RNA = rownames(tissue_df), Tissue = conversion$Tissue.Name, Map = conversion$Cell.System)
table(map$Map)

# Overlap between tissue and cells
length(Cell_source[Cell_source %in% paste0(conversion$Cell.System, " Cell System")])
