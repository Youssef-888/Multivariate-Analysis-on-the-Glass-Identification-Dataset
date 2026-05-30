# MACT 4233 - Project #3
# PCA and MDS of the Glass Identification Data
# AI help was used only to organize and simplify the report/code.
# The analysis decisions, interpretation, and final review remain my own work.

# Data

library(readxl)

download.file(
  "https://raw.githubusercontent.com/Youssef-888/Multivariate-Analysis-on-the-Glass-Identification-Dataset/main/Youssefdat.xlsx",
  destfile = "Youssefdat.xlsx",
  mode = "wb"
)

col_names <- c("Id", "RI", "Na", "Mg", "Al", "Si", "K", "Ca", "Ba", "Fe", "Type")

glass <- as.data.frame(
  read_excel("Youssefdat.xlsx", sheet = "GlassData")
)

names(glass) <- col_names
pca_vars <- c("RI", "Na", "Mg", "Al", "Si", "K", "Ca", "Ba", "Fe")
pair_vars <- c("RI", "Na", "Mg", "Al", "Ca", "Ba")

# Pairwise plot before PCA
pairs(glass[, pair_vars], col = glass$Type, pch = 19,
      main = "Selected Pairwise Plot of Glass Measurements")

# Classical PCA on the full data
x <- scale(glass[, pca_vars])
pc <- princomp(x, cor = TRUE)

summary(pc, loadings = TRUE)
plot(pc, main = "Classical PCA Scree Plot")
pairs(pc$scores[, 1:4], main = "Classical PCA Score Pairs")

pc_sd <- pc$sd
pc_var <- pc$sd^2
pc_prop <- pc_var / sum(pc_var)
pc_cum <- cumsum(pc_prop)

classical_pca_table <- round(data.frame(pc_sd, pc_var, pc_prop, pc_cum), 4)
names(classical_pca_table) <- c("Standard_deviation", "Variance", "Proportion", "Cumulative")
classical_pca_table

pc$loadings
round(cor(pc$scores[, 1:6]), 4)
round(cov(pc$scores[, 1:6]), 4)

# BACON result used in the report
# These IDs are the within-class BACON outliers. Fe was not used in the BACON step.
bacon_outlier_ids <- c(
  1, 22, 36, 48, 53, 54, 55, 56, 57, 62,
  71, 85, 93, 103, 104, 105, 106, 107, 108, 109,
  110, 111, 112, 113, 125, 128, 129, 130, 131, 132,
  164, 165, 174, 186, 187, 188, 189, 190, 202, 208, 212
)

bacon_summary <- data.frame(
  Type = c(1, 2, 3, 5, 6, 7),
  Group_size = c(70, 76, 17, 13, 9, 29),
  Outliers = c(10, 20, 0, 3, 0, 8),
  Cutoff = c(6.1784, 6.1080, 9.6647, 13.5519, 46.5561, 9.5915)
)

bacon_outlier <- glass$Id %in% bacon_outlier_ids
bacon_summary
glass[bacon_outlier, c("Id", "Type")]

pairs(glass[, pair_vars], col = ifelse(bacon_outlier, 2, glass$Type),
      pch = ifelse(bacon_outlier, 4, 19),
      main = "Pairwise Plot with BACON Outliers Highlighted")

# PCA after removing BACON outliers
glass_clean <- glass[!bacon_outlier, ]
x_clean <- scale(glass_clean[, pca_vars])
pc_clean <- princomp(x_clean, cor = TRUE)

summary(pc_clean, loadings = TRUE)
plot(pc_clean, main = "BACON-Clean PCA Scree Plot")

clean_sd <- pc_clean$sd
clean_var <- pc_clean$sd^2
clean_prop <- clean_var / sum(clean_var)
clean_cum <- cumsum(clean_prop)

clean_pca_table <- round(data.frame(clean_sd, clean_var, clean_prop, clean_cum), 4)
names(clean_pca_table) <- c("Standard_deviation", "Variance", "Proportion", "Cumulative")
clean_pca_table

pc_clean$loadings

pca_comparison <- data.frame(
  Analysis = c("Full data PCA", "BACON-clean PCA"),
  n = c(nrow(glass), nrow(glass_clean)),
  PC1_Eigenvalue = round(c(pc_var[1], clean_var[1]), 4),
  PC1_Proportion = round(c(pc_prop[1], clean_prop[1]), 4),
  Cum_PC2 = round(c(pc_cum[2], clean_cum[2]), 4),
  Cum_PC4 = round(c(pc_cum[4], clean_cum[4]), 4),
  PCs_for_95 = c(which(pc_cum >= 0.95)[1], which(clean_cum >= 0.95)[1])
)
pca_comparison

par(mfrow = c(1, 2))
plot(pc$scores[, 1], pc$scores[, 2], col = glass$Type, pch = 19,
     xlab = "PC1", ylab = "PC2", main = "Full Data PCA")
plot(pc_clean$scores[, 1], pc_clean$scores[, 2], col = glass_clean$Type, pch = 19,
     xlab = "PC1", ylab = "PC2", main = "BACON-Clean PCA")
par(mfrow = c(1, 1))

# MDS function
mds_analysis <- function(data, groups, title) {
  x <- scale(data)
  d <- dist(x)
  out <- cmdscale(d, eig = TRUE)
  
  eigenvalues <- out$eig[out$eig > 0]
  prop <- eigenvalues / sum(eigenvalues)
  Pm <- cumsum(prop)
  m <- which(Pm >= 0.80)[1]
  
  results <- round(data.frame(
    Dimension = 1:length(eigenvalues),
    Eigenvalue = eigenvalues,
    Proportion = prop,
    P_m = Pm
  ), 4)
  
  par(mfrow = c(1, 2))
  plot(eigenvalues, type = "b", pch = 19, main = paste(title, "Eigenvalues"),
       xlab = "Dimension", ylab = "Positive eigenvalue")
  plot(Pm, type = "b", pch = 19, main = paste(title, "Goodness"),
       xlab = "Number of dimensions", ylab = "P_m", ylim = c(0, 1))
  abline(h = 0.80, lty = 2)
  par(mfrow = c(1, 1))
  
  plot(out$points[, 1], out$points[, 2], col = groups, pch = 19,
       xlab = "MDS1", ylab = "MDS2", main = paste(title, "2D Map"))
  
  list(results = results, Pm = Pm, m = m, points = out$points)
}

# MDS on full and clean data
mds_full <- mds_analysis(glass[, pca_vars], glass$Type, "Full Data MDS")
mds_clean <- mds_analysis(glass_clean[, pca_vars], glass_clean$Type, "BACON-Clean MDS")

mds_full$results
mds_clean$results

mds_comparison <- data.frame(
  Analysis = c("Full data MDS", "BACON-clean MDS"),
  n = c(nrow(glass), nrow(glass_clean)),
  P_2 = round(c(mds_full$Pm[2], mds_clean$Pm[2]), 4),
  P_4 = round(c(mds_full$Pm[4], mds_clean$Pm[4]), 4),
  Dimensions_for_80 = c(mds_full$m, mds_clean$m)
)
mds_comparison
