rm(list = ls())
set.seed(7)

library(readxl)
library(dplyr)

# Load the data
raw_data <- read_excel("./data/raw_data/results_for_Dailin.xlsx") %>%
  as.data.frame()
# put the protein descriptions into the row names
rownames(raw_data) <- raw_data$Description
raw_data <- raw_data[, -1]

# create the metadata for samples
meta <- data.frame(
  row.names = colnames(raw_data),
  Sample = colnames(raw_data),
  Condition = c(rep("WT", 6), rep("KO", 6)),
  Time = c(2, 2, 6, 6, 18, 18,
           2, 2, 6, 6, 18, 18) # in hours
)

# check the max and min
max(na.omit(raw_data)) # 100
min(na.omit(raw_data)) # 0.01

# Replace entries that are 100 or 0.01 with NAs
raw_data[raw_data == 100] <- NA
raw_data[raw_data == 0.01] <- NA

# remove rows that are all NAs
data <- raw_data[rowSums(is.na(raw_data)) != ncol(raw_data), ]

# Split the data into WT and KO groups
WT <- data[, meta$Condition == "WT"]
KO <- data[, meta$Condition == "KO"]

# Keep rows that have at least 3 non-NA values
WT <- WT[rowSums(!is.na(WT)) >= 3, ]
KO <- KO[rowSums(!is.na(KO)) >= 3, ]

meta_WT <- meta[meta$Condition == "WT", ]
meta_KO <- meta[meta$Condition == "KO", ]

all(colnames(WT) == rownames(meta_WT)) # TRUE
all(colnames(KO) == rownames(meta_KO)) # TRUE

# log_e + 1 transformation
WT <- log(WT + 1)
KO <- log(KO + 1)

# Save the data
saveRDS(list(WT = WT, KO = KO, meta_WT = meta_WT, meta_KO = meta_KO), 
        "./data/processed_data_ls.RDS")
