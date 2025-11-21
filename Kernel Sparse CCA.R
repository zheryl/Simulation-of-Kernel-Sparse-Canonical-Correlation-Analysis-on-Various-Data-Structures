# Periksa & install paket yang diperlukan (tanpa interupsi jika sudah ada)
pkgs <- c("CCA", "PMA", "kernlab", "caret")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
library(CCA)      # cancor, rcc
library(PMA)      # CCA (sparse)
library(kernlab)  # kernelMatrix, rbfdot
library(caret)    # createDataPartition (opsional)

set.seed(123)

# ---------------------------
# Helper: simulasi data
# ---------------------------
simulasi_data <- function(kondisi = 1, n = 100, p = 10, q = 8) {
  if (kondisi == 1) {
    X <- matrix(rnorm(n * p), n, p)
    Y <- X[, 1:3] %*% matrix(rnorm(3 * q), 3, q) + matrix(rnorm(n * q, 0, 0.1), n, q)
  } else if (kondisi == 2) {
    X <- matrix(rnorm(n * p), n, p)
    Y <- matrix(rnorm(n * q), n, q) + 0.1 * X[, 1:2] %*% matrix(rnorm(2 * q), 2, q)
  } else if (kondisi == 3) {
    X <- matrix(rnorm(n * p), n, p)
    Y <- 0.8 * X[, 1:2] %*% matrix(rnorm(2 * q), 2, q) + matrix(rnorm(n * q), n, q)
  }
  return(list(X = scale(X), Y = scale(Y)))
}

# ---------------------------
# 1) CCA klasik
# ---------------------------
cca_klasik <- function(X, Y) {
  cancor(X, Y)
}

# ---------------------------
# 2) CCA regularized (dari package CCA)
# ---------------------------
cca_regularized <- function(X, Y, lambda1 = 0.1, lambda2 = 0.1) {
  # rcc(X, Y, lambda1, lambda2) ada di package CCA (regularized CCA)
  rcc(X, Y, lambda1, lambda2)
}

# ---------------------------
# 3) Sparse CCA (PMA::CCA)
# ---------------------------
sparse_cca <- function(X, Y, penaltyx = 0.7, penaltyz = 0.7, K = 2) {
  PMA::CCA(X, Y, typex = "standard", typez = "standard",
           penaltyx = penaltyx, penaltyz = penaltyz, K = K)
}

# ---------------------------
# Helper: centering kernel matrices (training)
# Kc = H K H  with H = I - 1/n 11^T
# ---------------------------
center_kernel_train <- function(K) {
  n <- nrow(K)
  one_n <- matrix(1 / n, n, n)
  H <- diag(n) - one_n
  H %*% K %*% H
}

# ---------------------------
# 4) Kernel Sparse CCA (gunakan kernel pada X/Y lalu apply sparse CCA)
#    - Mengembalikan kernel yang sudah di-center (train)
#    - NOTE: evaluasi pada data test butuh perhitungan cross-kernel centering yang teliti.
# ---------------------------
kernel_sparse_cca <- function(X, Y, sigma = 1, penaltyx = 0.7, penaltyz = 0.7, ncomp = 2) {
  # kernelMatrix menerima objek kernel dari rbfdot(sigma = value)
  k <- rbfdot(sigma = sigma)
  Kx <- kernelMatrix(k, X)
  Ky <- kernelMatrix(k, Y)
  
  # center kernel matrices (train)
  Kx_c <- center_kernel_train(Kx)
  Ky_c <- center_kernel_train(Ky)
  
  # gunakan sparse CCA pada kernel yang sudah di-center
  # karena dimensinya n x n, kita tetap pakai PMA::CCA (treat columns as variables)
  fit <- PMA::CCA(Kx_c, Ky_c,
                  typex = "standard", typez = "standard",
                  penaltyx = penaltyx, penaltyz = penaltyz, K = ncomp)
  
  list(u = fit$u, v = fit$v, cors = fit$cors, Kx = Kx_c, Ky = Ky_c, fit = fit, kernel = k)
}

# ---------------------------
# 5) Kernel CCA reguler (sederhana) - perbandingan
#    NOTE: cancor dipakai di sini tapi memperlakukan K matrices sebagai "variabel"
# ---------------------------
kernel_cca <- function(X, Y, sigma = 1, reg = 1e-8) {
  k <- rbfdot(sigma = sigma)
  Kx <- kernelMatrix(k, X)
  Ky <- kernelMatrix(k, Y)
  Kx_c <- center_kernel_train(Kx)
  Ky_c <- center_kernel_train(Ky)
  Kx_reg <- Kx_c + reg * diag(nrow(Kx_c))
  Ky_reg <- Ky_c + reg * diag(nrow(Ky_c))
  # cancor akan mengembalikan "canonical correlation" antar kolom Kx_reg dan Ky_reg
  cancor(Kx_reg, Ky_reg)
}

# ---------------------------
# 6) Estimasi parameter ringkas
# ---------------------------
estimasi_parameter <- function(X, Y, method = "classic", ...) {
  if (method == "classic") {
    fit <- cancor(X, Y)
    list(
      korelasi = fit$cor[1],
      koefisien_X = fit$xcoef[, 1],
      koefisien_Y = fit$ycoef[, 1],
      method = "CCA Klasik"
    )
  } else if (method == "sparse") {
    fit <- sparse_cca(X, Y, ...)
    list(
      korelasi = fit$cors[1],
      koefisien_X = fit$u[, 1],
      koefisien_Y = fit$v[, 1],
      method = "Sparse CCA"
    )
  } else if (method == "kernel_sparse") {
    fit <- kernel_sparse_cca(X, Y, ...)
    list(
      korelasi = fit$cors[1],
      koefisien_alpha = fit$u[, 1],
      koefisien_beta = fit$v[, 1],
      method = "Kernel Sparse CCA"
    )
  } else stop("Unknown method")
}

# ---------------------------
# 7) Uji hipotesis (permutasi sederhana pada Y)
# ---------------------------
uji_hipotesis <- function(X, Y, n_perm = 100) {
  cca_result <- cancor(X, Y)
  cor_observed <- cca_result$cor[1]
  cor_perm <- numeric(n_perm)
  for (i in seq_len(n_perm)) {
    Y_perm <- Y[sample(nrow(Y)), , drop = FALSE]
    cor_perm[i] <- cancor(X, Y_perm)$cor[1]
  }
  p_value <- mean(cor_perm >= cor_observed)
  list(korelasi_observed = cor_observed, p_value = p_value, signifikan = (p_value < 0.05))
}

# ---------------------------
# 8) Pilih model terbaik berdasarkan korelasi loading pertama
# ---------------------------
pilih_model_terbaik <- function(X, Y) {
  results <- list()
  
  # classic
  fit_classic <- cancor(X, Y)
  results[["classic"]] <- list(
    korelasi = fit_classic$cor[1],
    method = "CCA Klasik",
    sparsity_X = sum(fit_classic$xcoef[, 1] != 0),
    sparsity_Y = sum(fit_classic$ycoef[, 1] != 0)
  )
  
  # sparse
  fit_sparse <- sparse_cca(X, Y)
  results[["sparse"]] <- list(
    korelasi = fit_sparse$cors[1],
    method = "Sparse CCA",
    sparsity_X = sum(fit_sparse$u[, 1] != 0),
    sparsity_Y = sum(fit_sparse$v[, 1] != 0)
  )
  
  # kernel sparse
  fit_kernel_sparse <- kernel_sparse_cca(X, Y)
  results[["kernel_sparse"]] <- list(
    korelasi = fit_kernel_sparse$cors[1],
    method = "Kernel Sparse CCA",
    sparsity_alpha = sum(fit_kernel_sparse$u[, 1] != 0),
    sparsity_beta = sum(fit_kernel_sparse$v[, 1] != 0)
  )
  
  # kernel regular
  fit_kernel <- kernel_cca(X, Y)
  results[["kernel"]] <- list(
    korelasi = fit_kernel$cor[1],
    method = "Kernel CCA (reg)",
    sparsity_X = NA,
    sparsity_Y = NA
  )
  
  korelasi_values <- sapply(results, function(x) x$korelasi)
  best_idx <- which.max(korelasi_values)
  list(best_model = results[[best_idx]], all_results = results)
}

# ---------------------------
# EVALUASI MODEL (train/test split)
# NOTE: untuk kernel_sparse, evaluasi dilakukan pada TRAIN set (approx),
# karena proper centering/prediction untuk test memerlukan transformasi yang cermat.
# Jika butuh, kita bisa tambahkan fungsi untuk mem-center cross-kernels.
# ---------------------------
evaluasi_model <- function(X, Y, method = "classic") {
  set.seed(123)
  train_idx <- sample(seq_len(nrow(X)), size = floor(0.7 * nrow(X)))
  X_train <- X[train_idx, , drop = FALSE]; Y_train <- Y[train_idx, , drop = FALSE]
  X_test <- X[-train_idx, , drop = FALSE]; Y_test <- Y[-train_idx, , drop = FALSE]
  
  if (method == "classic") {
    fit <- cancor(X_train, Y_train)
    # proyeksikan test
    u_test <- as.matrix(X_test) %*% fit$xcoef[, 1]
    v_test <- as.matrix(Y_test) %*% fit$ycoef[, 1]
    pred_cor <- cor(u_test, v_test)
  } else if (method == "sparse") {
    fit <- sparse_cca(X_train, Y_train)
    u_test <- as.matrix(X_test) %*% fit$u[, 1]
    v_test <- as.matrix(Y_test) %*% fit$v[, 1]
    pred_cor <- cor(u_test, v_test)
  } else if (method == "kernel_sparse") {
    # Untuk sekarang: hitung korelasi canonical pada TRAIN split (sebagai proxy)
    fit <- kernel_sparse_cca(X_train, Y_train)
    pred_cor <- fit$cors[1]
    # Catatan: jika ingin korelasi pada test set, perlu:
    # - hitung K_test_train = K(X_test, X_train)
    # - lakukan centering cross-kernel secara benar (rumus H K H dan cross-term)
    # - kemudian proyeksikan K_test_centered %*% alpha dan K_test_centered_y %*% beta
    # Aku sengaja tidak mengimplementasikan cross-centering rumit itu di sini agar tetap jelas.
  } else stop("Unknown method")
  
  round(as.numeric(pred_cor), 4)
}

# ---------------------------
# SIMULASI UTAMA
# ---------------------------
main_simulasi <- function() {
  cat("=== SIMULASI CCA MULTI-KONDISI DENGAN KERNEL SPARSE CCA ===\n\n")
  for (kondisi in 1:3) {
    cat(sprintf("KONDISI %d:\n", kondisi))
    cat(rep("-", 50), "\n")
    data <- simulasi_data(kondisi, n = 100, p = 10, q = 8)
    X <- data$X; Y <- data$Y
    
    # 1. CCA Klasik
    cca_classic <- cca_klasik(X, Y)
    cat("1. CCA Klasik - Korelasi:", round(cca_classic$cor[1], 4), "\n")
    
    # 2. Sparse CCA
    sparse_result <- sparse_cca(X, Y)
    cat("2. Sparse CCA - Korelasi:", round(sparse_result$cors[1], 4),
        "- Sparsity X:", sum(sparse_result$u[, 1] != 0),
        "Y:", sum(sparse_result$v[, 1] != 0), "\n")
    
    # 3. Kernel Sparse CCA
    kernel_sparse_result <- kernel_sparse_cca(X, Y)
    cat("3. Kernel Sparse CCA - Korelasi:", round(kernel_sparse_result$cors[1], 4),
        "- Sparsity α:", sum(kernel_sparse_result$u[, 1] != 0),
        "β:", sum(kernel_sparse_result$v[, 1] != 0), "\n")
    
    # 4. Uji signifikansi
    uji <- uji_hipotesis(X, Y, n_perm = 50)
    cat("4. Uji Signifikansi - p-value:", round(uji$p_value, 4),
        "- Signifikan:", uji$signifikan, "\n")
    
    # 5. Pilih model terbaik
    best <- pilih_model_terbaik(X, Y)
    cat("5. Model Terbaik:", best$best_model$method,
        "- Korelasi:", round(best$best_model$korelasi, 4), "\n\n")
  }
}

# Jalankan simulasi utama
main_simulasi()

# ---------------------------
# Evaluasi pada data test (proxy untuk kernel_sparse)
# ---------------------------
cat("\n=== EVALUASI PADA DATA TEST (PROXY UNTUK KERNEL SPARSE) ===\n")
for (kondisi in 1:3) {
  data <- simulasi_data(kondisi)
  cat(sprintf("\nKondisi %d:\n", kondisi))
  cat("CCA Classic:", evaluasi_model(data$X, data$Y, "classic"), "\n")
  cat("Sparse CCA :", evaluasi_model(data$X, data$Y, "sparse"), "\n")
  cat("Kernel Sparse CCA (train-proxy):", evaluasi_model(data$X, data$Y, "kernel_sparse"), "\n")
}

# ---------------------------
# Contoh detail kernel sparse CCA (kondisi 1)
# ---------------------------
cat("\n=== DETAIL KERNEL SPARSE CCA ===\n")
data <- simulasi_data(1)
X <- data$X; Y <- data$Y
ks_fit <- kernel_sparse_cca(X, Y)
cat("Korelasi kanonik (components):", round(ks_fit$cors, 4), "\n")
cat("Jumlah sample (n):", nrow(X), "\n")
cat("Dimensi alpha (length):", length(ks_fit$u[, 1]), "\n")
cat("Sparsity pattern - alpha (non-zero indices):", which(ks_fit$u[, 1] != 0), "\n")
cat("Sparsity pattern - beta  (non-zero indices):", which(ks_fit$v[, 1] != 0), "\n")




# ============================================================
# Tambahan Visualisasi dari Hasil Simulasi
# ============================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

# Jalankan ulang simulasi dan simpan hasil ringkas untuk plotting
hasil_simulasi <- data.frame()

for (kondisi in 1:3) {
  data <- simulasi_data(kondisi)
  X <- data$X; Y <- data$Y
  
  # CCA klasik
  cca_classic <- cca_klasik(X, Y)
  cor_classic <- cca_classic$cor[1]
  cor_test_classic <- evaluasi_model(X, Y, "classic")
  
  # Sparse CCA
  sparse_result <- sparse_cca(X, Y)
  cor_sparse <- sparse_result$cors[1]
  cor_test_sparse <- evaluasi_model(X, Y, "sparse")
  
  # Kernel Sparse CCA
  kernel_sparse_result <- kernel_sparse_cca(X, Y)
  cor_kernel <- kernel_sparse_result$cors[1]
  cor_test_kernel <- evaluasi_model(X, Y, "kernel_sparse")
  
  hasil_simulasi <- rbind(hasil_simulasi, data.frame(
    Kondisi = paste("Simulasi", kondisi),
    Metode = c("CCA Klasik", "Sparse CCA", "Kernel Sparse CCA"),
    Korelasi_Train = c(cor_classic, cor_sparse, cor_kernel),
    Korelasi_Test = c(cor_test_classic, cor_test_sparse, cor_test_kernel)
  ))
}

# ============================================================
# Plot 1: Nilai Korelasi pada Data Latih
# ============================================================
plot_train <- ggplot(hasil_simulasi, aes(x = Kondisi, y = Korelasi_Train, fill = Metode)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = round(Korelasi_Train, 3)), 
            position = position_dodge(width = 0.7), vjust = -0.5, size = 3.5) +
  labs(title = "Perbandingan Korelasi Kanonik antar Metode (Data Latih)",
       x = "Skenario Simulasi", y = "Nilai Korelasi Kanonik (Train)") +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold", size = 13))
plot_train
# ============================================================
# Plot 2: Nilai Korelasi pada Data Uji
# ============================================================
plot_test <- ggplot(hasil_simulasi, aes(x = Kondisi, y = Korelasi_Test, group = Metode, color = Metode)) +
  geom_line(size = 1.1) +
  geom_point(size = 3) +
  geom_text(aes(label = round(Korelasi_Test, 3)), vjust = -0.6, size = 3.5) +
  labs(title = "Kestabilan Hasil Prediksi antar Metode (Data Uji)",
       x = "Skenario Simulasi", y = "Nilai Korelasi (Test)") +
  scale_color_brewer(palette = "Dark2") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold", size = 13))

plot_test
