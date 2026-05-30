library(readxl)

data <- as.data.frame(read_excel("Youssefdat.xlsx", sheet = "GlassData"))
features <- setdiff(names(data), c("Id", "Type"))
data <- data[complete.cases(data[, c("Type", features)]), ]

X <- scale(data[, features])
k <- length(unique(data$Type))

cut_line <- function(hc, k) {
  n <- length(hc$order)
  (hc$height[n - k] + hc$height[n - k + 1]) / 2
}

plot_tree <- function(hc, title) {
  plot(hc, labels = FALSE, hang = -1,
       main = title, xlab = "Observation Index", ylab = "Height")
  abline(h = cut_line(hc, k), col = "red", lty = 2)
  legend("topright", paste("Cut at k =", k), col = "red", lty = 2, bty = "n")
}

# Figure 1
hc1 <- hclust(dist(X, method = "euclidean"), method = "ward.D2")
plot_tree(hc1, "Dendrogram: Euclidean + Ward.D2")

# Figure 2
hc2 <- hclust(dist(X, method = "euclidean"), method = "average")
plot_tree(hc2, "Dendrogram: Euclidean + Average")

# Figure 3
hc3 <- hclust(dist(X, method = "manhattan"), method = "ward.D2")
plot_tree(hc3, "Dendrogram: Manhattan + Ward.D2")

# Figure 4
hc4 <- hclust(dist(X, method = "manhattan"), method = "complete")
plot_tree(hc4, "Dendrogram: Manhattan + Complete")

# Figure 5
set.seed(123)
k_values <- 1:10
wss <- numeric(length(k_values))

for (i in seq_along(k_values)) {
  km <- kmeans(X, centers = k_values[i], nstart = 50)
  wss[i] <- km$tot.withinss
}

plot(k_values, wss, type = "b", pch = 19,
     main = "K-means Elbow Curve",
     xlab = "Number of Clusters (k)",
     ylab = "Total Within-Cluster SS (WSS)",
     xaxt = "n")
axis(1, at = k_values)
abline(v = k, col = "red", lty = 2)
legend("topright", paste("Chosen k =", k), col = "red", lty = 2, bty = "n")

# Figure 6
set.seed(129)
km_final <- kmeans(X, centers = k, nstart = 50)
pca <- prcomp(X)
scores <- pca$x[, 1:2]
var_exp <- pca$sdev^2 / sum(pca$sdev^2)

plot(scores[, 1], scores[, 2],
     col = km_final$cluster, pch = 19,
     main = paste0("K-means Cluster Plot (k = ", k,
                   ", WSS = ", round(km_final$tot.withinss, 1), ")"),
     xlab = paste0("Dim1 (", round(100 * var_exp[1], 1), "%)"),
     ylab = paste0("Dim2 (", round(100 * var_exp[2], 1), "%)"))

for (g in sort(unique(km_final$cluster))) {
  pts <- scores[km_final$cluster == g, , drop = FALSE]
  if (nrow(pts) >= 3) {
    h <- chull(pts)
    polygon(pts[c(h, h[1]), ], border = g, col = adjustcolor(g, alpha.f = 0.15))
  }
}

points(scores[, 1], scores[, 2], col = km_final$cluster, pch = 19)
legend("topright", legend = sort(unique(km_final$cluster)),
       col = sort(unique(km_final$cluster)), pch = 19, title = "cluster")
