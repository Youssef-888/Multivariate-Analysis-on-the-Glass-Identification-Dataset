library(MASS)
library(nnet)
library(robustX)
library(robustbase)
library(readxl)

cols <- c("Id","RI","Na","Mg","Al","Si","K","Ca","Ba","Fe","Type")
glass <- as.data.frame(read_excel("Youssefdat.xlsx"))
colnames(glass) <- cols

type_names <- c(
  "1" = "building_windows_float",
  "2" = "building_windows_nonfloat",
  "3" = "vehicle_windows_float",
  "5" = "containers",
  "6" = "tableware",
  "7" = "headlamps"
)

vars <- c("RI","Na","Mg","Al","Si","K","Ca","Ba")
pair_vars <- c("Na","Mg","Al","Ba")
class_levels <- sort(unique(glass$Type))

err <- function(y, p) 100 * mean(y != p)
cm <- function(y, p) table(factor(y, levels = class_levels), factor(p, levels = class_levels))

scale_xy <- function(x) scale(as.matrix(x))

run_bacon <- function(x, ids = NULL, main = "BACON Distances") {
  out <- try(mvBACON(as.matrix(x)), silent = TRUE)
  if (inherits(out, "try-error")) {
    cat(main, "- BACON failed\n")
    return(NULL)
  }
  y <- cbind(1:nrow(x), out$dis)
  colnames(y) <- c("Index", "Distance")
  plot(y, pch = 19, main = main, xlab = "Index", ylab = "BACON Distance")
  abline(h = out$limit, col = "red", lty = 2)
  points(y[!out$subset, , drop = FALSE], pch = 4, col = 2, cex = 1.5)
  if (!is.null(ids)) {
    cat("Outlier IDs:", if (any(!out$subset)) paste(ids[!out$subset], collapse = ", ") else "None", "\n")
  }
  out
}

flda <- function(x, class) {
  cat("Fisher Linear Discriminant:\n")
  a <- lda(x, class)
  d <- predict(a)
  t <- table(class, d$class)
  print(t)
  er <- 100 * (sum(t) - sum(diag(t))) / nrow(x)
  cat("Error Rate =", er, "%\n")
  invisible(d)
}

flda2 <- function(x, class) {
  if (ncol(x) != 2) {
    cat("Data should be 2-dimensional\n")
    return(NULL)
  }
  t <- factor(class)
  level <- levels(t)
  if (length(level) != 2) {
    cat("Data should have only two groups\n")
    return(NULL)
  }
  y <- x[class == level[1], ]
  x <- x[class == level[2], ]
  n1 <- nrow(x)
  n2 <- nrow(y)
  n <- n1 + n2
  xcenter <- colMeans(x)
  ycenter <- colMeans(y)
  xcov <- cov(x)
  ycov <- cov(y)
  sp <- ((n1 - 1) * xcov + (n2 - 1) * ycov) / (n - 2)
  d <- xcenter - ycenter
  m <- (xcenter + ycenter) / 2
  a <- solve(sp) %*% d
  class2 <- c(rep(1, n1), rep(2, n2))
  p <- 1
  z <- rbind(x, y)
  pred <- z - matrix(m, ncol = 2, nrow = n, byrow = TRUE)
  pred <- as.matrix(pred)
  pred <- (pred %*% a) < log(p)
  C <- (class2 != pred + 1)
  ce <- sum(C)
  cat("--------------------------------------------------\n")
  cat(" Correct Incorrect\n Class Classification Classification Total\n")
  cat("--------------------------------------------------\n")
  cd1 <- n1 - sum(C[1:n1])
  cat(" 1 ", cd1, " ", n1 - cd1, " ", n1, "\n")
  cd2 <- n2 - sum(C[(n1 + 1):n])
  cat(" 2 ", cd2, " ", n2 - cd2, " ", n2, "\n")
  cat("--------------------------------------------------\n")
  cat(" Total: ", cd1 + cd2, " ", n - (cd1 + cd2), " ", n, "\n")
  cat("--------------------------------------------------\n")
  cat("Error Rate = ", 100 * (ce / n), "%\n")
  const <- (sum(a * m) + log(p)) / a[2]
  slope <- -a[1] / a[2]
  plot(z[, 1:2], col = class2, pch = 19,
       xlab = colnames(z)[1], ylab = colnames(z)[2],
       main = "FLDA2")
  abline(const, slope, col = "blue")
  points(rbind(xcenter[1:2], ycenter[1:2]), pch = 19, col = 3, cex = 1.5)
  segments(xcenter[1], xcenter[2], ycenter[1], ycenter[2], col = 3)
  invisible(list(xcenter = xcenter[2:1], ycenter = ycenter[2:1], xcov = xcov,
                 ycov = ycov, sp = sp, a = a, slope = slope, const = const,
                 ce = ce, m = m, z = z))
}

loo_multinom <- function(df, vars) {
  pred <- numeric(nrow(df))
  for (i in 1:nrow(df)) {
    train <- df[-i, ]
    test <- df[i, , drop = FALSE]
    xtrain <- scale_xy(train[, vars])
    center <- attr(xtrain, "scaled:center")
    scalev <- attr(xtrain, "scaled:scale")
    xtest <- scale(as.matrix(test[, vars]), center = center, scale = scalev)
    train2 <- data.frame(Type = factor(train$Type), xtrain)
    test2 <- data.frame(xtest)
    names(train2) <- c("Type", vars)
    names(test2) <- vars
    fit <- multinom(Type ~ ., data = train2, trace = FALSE)
    pred[i] <- as.numeric(as.character(predict(fit, test2)))
  }
  pred
}

pairs(glass[, pair_vars], col = glass$Type, pch = 19)

full_bacon <- run_bacon(glass[, vars], glass$Id, "BACON Distances: Full Data")

class_outliers <- rep(FALSE, nrow(glass))
for (g in class_levels) {
  idx <- which(glass$Type == g)
  cat("\nType", g, "-", type_names[as.character(g)], "\n")
  out <- run_bacon(glass[idx, vars], glass$Id[idx], paste("BACON Distances: Type", g))
  if (!is.null(out)) {
    class_outliers[idx[!out$subset]] <- TRUE
  }
}

par(mfrow = c(2, 4))
for (v in vars) {
  qqnorm(glass[[v]], main = paste("Q-Q Plot:", v), pch = 19)
  qqline(glass[[v]], col = 2)
}
par(mfrow = c(1, 1))

x_full <- scale_xy(glass[, vars])
rslt_flda_full <- flda(x_full, glass$Type)
pred_flda_full <- as.numeric(as.character(rslt_flda_full$class))
loo_flda_full <- as.numeric(as.character(lda(x_full, glass$Type, CV = TRUE)$class))

mn_full_data <- data.frame(Type = factor(glass$Type), x_full)
names(mn_full_data) <- c("Type", vars)
mn_full <- multinom(Type ~ ., data = mn_full_data, trace = FALSE)
pred_mn_full <- as.numeric(as.character(predict(mn_full)))
loo_mn_full <- loo_multinom(glass, vars)

glass_clean <- glass[!class_outliers, ]
x_clean <- scale_xy(glass_clean[, vars])
rslt_flda_clean <- flda(x_clean, glass_clean$Type)
pred_flda_clean <- as.numeric(as.character(rslt_flda_clean$class))
loo_flda_clean <- as.numeric(as.character(lda(x_clean, glass_clean$Type, CV = TRUE)$class))

mn_clean_data <- data.frame(Type = factor(glass_clean$Type), x_clean)
names(mn_clean_data) <- c("Type", vars)
mn_clean <- multinom(Type ~ ., data = mn_clean_data, trace = FALSE)
pred_mn_clean <- as.numeric(as.character(predict(mn_clean)))
loo_mn_clean <- loo_multinom(glass_clean, vars)

sub2 <- glass[glass$Type %in% c(1, 5), ]
x2 <- scale_xy(sub2[, c("Mg", "Al")])
colnames(x2) <- c("Scaled Mg", "Scaled Al")
cat("\nFLDA2: Type 1 vs Type 5\n")
rslt_flda2 <- flda2(x2, sub2$Type)
loo_flda2 <- as.numeric(as.character(lda(x2, sub2$Type, CV = TRUE)$class))
fit_flda2 <- lda(x2, sub2$Type)
pred_flda2 <- as.numeric(as.character(predict(fit_flda2, x2)$class))

cat("\nFull data BACON outliers:\n")
if (!is.null(full_bacon)) cat(sum(!full_bacon$subset), "observation(s) flagged\n")
cat("Within-class BACON removed:", sum(class_outliers), "observation(s)\n\n")

cat("FLDA full internal error:", round(err(glass$Type, pred_flda_full), 2), "%\n")
print(cm(glass$Type, pred_flda_full))
cat("\nFLDA full LOO error:", round(err(glass$Type, loo_flda_full), 2), "%\n")
print(cm(glass$Type, loo_flda_full))

cat("\nMultinomial full internal error:", round(err(glass$Type, pred_mn_full), 2), "%\n")
print(cm(glass$Type, pred_mn_full))
cat("\nMultinomial full LOO error:", round(err(glass$Type, loo_mn_full), 2), "%\n")
print(cm(glass$Type, loo_mn_full))

cat("\nFLDA clean internal error:", round(err(glass_clean$Type, pred_flda_clean), 2), "%\n")
print(cm(glass_clean$Type, pred_flda_clean))
cat("\nFLDA clean LOO error:", round(err(glass_clean$Type, loo_flda_clean), 2), "%\n")
print(cm(glass_clean$Type, loo_flda_clean))

cat("\nMultinomial clean internal error:", round(err(glass_clean$Type, pred_mn_clean), 2), "%\n")
print(cm(glass_clean$Type, pred_mn_clean))
cat("\nMultinomial clean LOO error:", round(err(glass_clean$Type, loo_mn_clean), 2), "%\n")
print(cm(glass_clean$Type, loo_mn_clean))

cat("\nFLDA2 internal error:", round(err(sub2$Type, pred_flda2), 2), "%\n")
print(table(True = sub2$Type, Predicted = pred_flda2))
cat("\nFLDA2 LOO error:", round(err(sub2$Type, loo_flda2), 2), "%\n")
print(table(True = sub2$Type, Predicted = loo_flda2))
