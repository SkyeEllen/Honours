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
# Source for all data
bio_source <- read.csv(here("RNA Splicing Data", "RNASeq_Atlas_samples.csv"))
bio_source$RNA_number_id <- paste0("X", bio_source$RNA_number)

bio_source_keep <- bio_source %>% select(RNA_number_id, Biological_source)
saveRDS(bio_source_keep, here("RNA Splicing Data", "Bio_source_all.RDS"))

# Source for cell data
cell_source_df <- data.frame(RNA_number_id = rownames(cell_df))
cell_source_df <- left_join(cell_source_df,
                            bio_source %>% select(RNA_number_id, Biological_source))
saveRDS(cell_source_df, here("RNA Splicing Data", "Cell source.RDS"))

# Source for tissue data
tissue_source_df <- data.frame(RNA_number_id = rownames(tissue_df))
tissue_source_df <- left_join(tissue_source_df,
                              bio_source %>% select(RNA_number_id, Biological_source))
conversion <- read.csv(here("RNA Splicing Data", "Tissue to Cell Sys.csv"))
tissue_source_df <- left_join(tissue_source_df,
                              conversion %>% select(Biological_source, Cell_system))
saveRDS(tissue_source_df, here("RNA Splicing Data", "Tissue source.RDS"))
