rm(list = ls())
set.seed(7)

library(dplyr)

source("./fcn/fcn_metric.R")
# Load the data_ls
data_ls <- readRDS("./data/processed_data_ls.RDS")

# Get the results for WT
tcc <- 26.0 # pre-specified cell cycle value for WT
WT <- data_ls$WT
WT_meta <- data_ls$meta_WT
time <- WT_meta$Time # shared by all proteins

slope_beta <- c()
T_HalfLife <- c()
R2 <- c()
CVSQ <- c()
n_points_used <- c()
for (i in 1:nrow(WT)) {
  protein_i <- WT[i, ] %>%
    unlist()
  metrics_i <- metrics(protein_i, time, tcc)
  slope_beta <- c(slope_beta, metrics_i$slope_beta)
  T_HalfLife <- c(T_HalfLife, metrics_i$T_HalfLife)
  R2 <- c(R2, metrics_i$R2)
  CVSQ <- c(CVSQ, metrics_i$CVSQ)
  n_points_used <- c(n_points_used, metrics_i$n_pt)
}
WT_res_df <- data.frame(
  Protein = rownames(WT),
  Slope = slope_beta,
  T_HalfLife = T_HalfLife,
  R2 = R2,
  CVSQ = CVSQ,
  n_points_used = n_points_used
)
write.csv(WT_res_df, "./res/2024_0912/WT_results.csv", 
          row.names = FALSE)

# Get the results for KO
tcc <- 21.7 # pre-specified cell cycle value for KO
KO <- data_ls$KO
KO_meta <- data_ls$meta_KO
time <- KO_meta$Time # shared by all proteins

slope_beta <- c()
T_HalfLife <- c()
R2 <- c()
CVSQ <- c()
n_points_used <- c()

for (i in 1:nrow(KO)) {
  protein_i <- KO[i, ] %>%
    unlist()
  metrics_i <- metrics(protein_i, time, tcc)
  slope_beta <- c(slope_beta, metrics_i$slope_beta)
  T_HalfLife <- c(T_HalfLife, metrics_i$T_HalfLife)
  R2 <- c(R2, metrics_i$R2)
  CVSQ <- c(CVSQ, metrics_i$CVSQ)
  n_points_used <- c(n_points_used, metrics_i$n_pt)
}
KO_res_df <- data.frame(
  Protein = rownames(KO),
  Slope = slope_beta,
  T_HalfLife = T_HalfLife,
  R2 = R2,
  CVSQ = CVSQ,
  n_points_used = n_points_used
)
write.csv(KO_res_df, "./res/2024_0912/KO_results.csv", 
          row.names = FALSE)
