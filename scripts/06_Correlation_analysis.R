#set working directory
setwd("D:/Thesis/work_package_3/Data_analysis/01_figures/08_correlation_analysis/new_correlation_analysis_04102025")

# ============================
# Temperature monotonicity + smoothed-curve co-occurrence (σ = 8)
# Safe correlations (no NA) and same heatmap style
# ============================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(rlang)
  library(Hmisc); library(RColorBrewer)
})

# ---------- Parameters ----------
gradient_col <- "sample_temp"
taxon_col    <- "GTDB_Tk_Phylum"
sigma_bins   <- 8
n_perm       <- 5000
alpha_lbl    <- 0.05
set.seed(7)

# ---------- Read & tidy ----------
dat <- read.csv("cbb_data.csv", stringsAsFactors = FALSE) %>%
  filter(.data[[gradient_col]] != "Unknown", !is.na(.data[[gradient_col]])) %>%
  mutate(!!taxon_col := trimws(as.character(.data[[taxon_col]])))

# Parse gradient to numeric (supports ranges like "0–200")
parse_num <- function(x) {
  x <- as.character(x)
  val <- suppressWarnings(as.numeric(x))
  if (!is.na(val)) return(val)
  s <- gsub(",", "", x)
  nums <- regmatches(s, gregexpr("\\d+\\.?\\d*", s))[[1]]
  if (length(nums) >= 2 && grepl("-|–|—| to ", s)) {
    return(mean(as.numeric(nums[1:2])))
  } else if (length(nums) >= 1) {
    return(as.numeric(nums[1]))
  }
  NA_real_
}
dat$grad_num <- vapply(dat[[gradient_col]], parse_num, numeric(1))
if (mean(is.na(dat$grad_num)) > 0.5) {
  levs <- unique(dat[[gradient_col]])
  dat$grad_num <- as.numeric(factor(dat[[gradient_col]], levels = levs))
}

# ---------- Presence grid ----------
pres_long <- dat %>%
  distinct(grad_num, !!sym(taxon_col)) %>%
  mutate(Presence = 1L) %>%
  complete(
    grad_num = sort(unique(grad_num)),
    !!sym(taxon_col),
    fill = list(Presence = 0L)
  ) %>%
  arrange(grad_num, !!sym(taxon_col))

A_wide <- pres_long %>%
  tidyr::pivot_wider(names_from = !!sym(taxon_col), values_from = Presence, values_fill = 0) %>%
  arrange(grad_num)

grad_vec <- A_wide$grad_num
A <- as.matrix(A_wide[, -1, drop = FALSE])
phyla <- colnames(A)
L <- nrow(A); P <- ncol(A)

# ---------- Monotonicity per phylum ----------
check_monotone <- function(v, tol = 0) {
  dv <- diff(v)
  inc <- all(dv >= -tol); dec <- all(dv <= tol)
  if (inc && !dec) return(c(TRUE, "non-decreasing"))
  if (dec && !inc) return(c(TRUE, "non-increasing"))
  if (inc &&  dec) return(c(TRUE, "constant"))
  c(FALSE, "non-monotonic")
}
mono_tbl <- lapply(seq_len(P), function(j){
  v  <- A[, j]
  md <- check_monotone(v)
  sp <- suppressWarnings(cor.test(rank(grad_vec), v, method = "spearman"))
  data.frame(
    Phylum        = phyla[j],
    monotone      = as.logical(md[1]),
    direction     = md[2],
    spearman_rho  = unname(sp$estimate),
    spearman_p    = sp$p.value,
    n_levels      = length(v),
    prevalence    = sum(v),
    fraction_pres = mean(v),
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()
write.csv(mono_tbl, "temperature_monotonicity_by_phylum.csv", row.names = FALSE)

# ---------- Smoothing (Gaussian σ = 8) ----------
gaussian_kernel <- function(sigma_bins) {
  half <- max(1L, round(3 * sigma_bins))
  x <- -half:half
  k <- exp(-(x^2) / (2 * sigma_bins^2))
  k / sum(k)
}
smooth_reflect <- function(y, k) {
  pad <- floor((length(k) - 1) / 2)
  ypad <- c(rev(y[2:(pad+1)]), y, rev(y[(length(y)-pad):(length(y)-1)]))
  yy <- stats::filter(ypad, k, sides = 2)
  as.numeric(yy[(pad+1):(pad+length(y))])
}
k <- gaussian_kernel(sigma_bins)
Z <- apply(A, 2, smooth_reflect, k = k)

# Standardize columns; constant columns become all-zero vectors
Z <- scale(Z, center = TRUE, scale = TRUE)
Z[is.na(Z)] <- 0

# Identify constant / near-constant columns post-smoothing
col_sd <- apply(Z, 2, sd)
const_cols <- phyla[col_sd == 0]
if (length(const_cols)) {
  message("Constant (zero-variance) after smoothing: ",
          paste(const_cols, collapse = ", "))
}

# ---------- Safe correlation matrix (no NA) ----------
# Because Z columns are mean-zero and unit-sd (or zero if constant),
# R = Z'Z / L gives correlations; constant columns simply give 0 with others.
R <- (t(Z) %*% Z) / L
diag(R) <- 1
R <- as.matrix(R)
rownames(R) <- colnames(R) <- phyla

# ---------- Circular-shift permutation p-values ----------
perm_pvalue_corr <- function(x, y, L, n_perm) {
  shifts <- sample.int(L, n_perm, replace = TRUE) - 1L
  idx <- (rep(0:(L-1), times = n_perm) - rep(shifts, each = L)) %% L + 1L
  yperm <- matrix(y[idx], nrow = L, ncol = n_perm)
  r_perm <- colSums(x * yperm) / L
  r_obs  <- sum(x * y) / L
  p <- (1 + sum(abs(r_perm) >= abs(r_obs))) / (1 + n_perm)
  c(r_obs = r_obs, p = p)
}

pairs <- list(); cnt <- 1L
for (i in 1:(P-1)) {
  xi <- Z[, i]
  for (j in (i+1):P) {
    out <- perm_pvalue_corr(xi, Z[, j], L, n_perm)
    pairs[[cnt]] <- data.frame(
      Phylum1     = phyla[i],
      Phylum2     = phyla[j],
      curve_corr  = out["r_obs"],
      p_value     = out["p"],
      stringsAsFactors = FALSE
    ); cnt <- cnt + 1L
  }
}
pairs_df <- bind_rows(pairs)

# BH-FDR (optional)
o <- order(pairs_df$p_value)
r <- integer(length(o)); r[o] <- seq_along(o)
q <- pairs_df$p_value * nrow(pairs_df) / r
q[o] <- rev(cummin(rev(q[o])))
pairs_df$q_value <- q
write.csv(pairs_df, "temperature_curve_pairs_sigma8_perm5000.csv", row.names = FALSE)

# ---------- Build heatmap data (same columns as your Spearman code) ----------
# Fill symmetric matrices from pairs and R (already safe)
R_mat <- R
P_mat <- matrix(0, P, P, dimnames = list(phyla, phyla))
apply(pairs_df, 1, function(row){
  i <- row[["Phylum1"]]; j <- row[["Phylum2"]]
  P_mat[i, j] <<- as.numeric(row[["p_value"]])
  P_mat[j, i] <<- as.numeric(row[["p_value"]])
})

cor_df  <- as.data.frame(as.table(R_mat))
pval_df <- as.data.frame(as.table(P_mat))
colnames(cor_df)  <- c("Phylum1","Phylum2","Spearman_r")  # reuse plot mapping
colnames(pval_df) <- c("Phylum1","Phylum2","p_value")

merged_df <- inner_join(cor_df, pval_df, by = c("Phylum1","Phylum2")) %>%
  filter(Phylum1 != Phylum2)

write.csv(merged_df, "CBB_temperature_curvecorr_pvalue.csv", row.names = FALSE)

# ---------- Heatmap (same aesthetics) ----------
div_pal <- brewer.pal(n = 11, name = "RdBu")
heatmap_data <- merged_df %>%
  mutate(label = ifelse(p_value < alpha_lbl, round(Spearman_r, 2), NA))

temp_corr <- ggplot(heatmap_data, aes(x = Phylum1, y = Phylum2, fill = Spearman_r)) +
  geom_tile(color = "grey90", size = 0.3) +
  scale_fill_gradient2(
    low = div_pal[1], mid = div_pal[6], high = div_pal[11],
    midpoint = 0, limits = c(-1, 1), space = "Lab",
    na.value = "white",
    name = "Curve corr\n(σ=8)"
  ) +
  geom_text(aes(label = label), color = "black", size = 3.5, fontface = "bold", na.rm = TRUE) +
  coord_fixed() +
  labs(title = "", subtitle = "", x = NULL, y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x      = element_text(angle = 45, vjust = 1, hjust = 1, face = "bold"),
    axis.text.y      = element_text(face = "bold"),
    legend.title     = element_text(face = "bold", size = 12),
    legend.text      = element_text(size = 10),
    panel.grid       = element_blank(),
    panel.border     = element_rect(fill = NA, color = "grey50"),
    axis.ticks       = element_blank(),
    plot.margin      = margin(t = 20, r = 20, b = 20, l = 20)
  )

print(temp_corr)
ggsave("CBB_temperature_curvecorr_heatmap_sigma8.png", temp_corr, width = 10, height = 8, dpi = 300)


# ============================
# Depth monotonicity + smoothed-curve co-occurrence (σ = 8)
# Safe correlations (no NA) and same heatmap style
# ============================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(rlang)
  library(Hmisc); library(RColorBrewer)
})

# ---------- Parameters ----------
gradient_col <- "sample_depth"
taxon_col    <- "GTDB_Tk_Phylum"
sigma_bins   <- 8
n_perm       <- 5000
alpha_lbl    <- 0.05
set.seed(7)

# ---------- Read & tidy ----------
dat <- read.csv("cbb_data.csv", stringsAsFactors = FALSE) %>%
  filter(.data[[gradient_col]] != "Unknown", !is.na(.data[[gradient_col]])) %>%
  mutate(!!taxon_col := trimws(as.character(.data[[taxon_col]])))

# Parse gradient to numeric (supports ranges like "0–200")
parse_num <- function(x) {
  x <- as.character(x)
  val <- suppressWarnings(as.numeric(x))
  if (!is.na(val)) return(val)
  s <- gsub(",", "", x)
  nums <- regmatches(s, gregexpr("\\d+\\.?\\d*", s))[[1]]
  if (length(nums) >= 2 && grepl("-|–|—| to ", s)) {
    return(mean(as.numeric(nums[1:2])))
  } else if (length(nums) >= 1) {
    return(as.numeric(nums[1]))
  }
  NA_real_
}
dat$grad_num <- vapply(dat[[gradient_col]], parse_num, numeric(1))
if (mean(is.na(dat$grad_num)) > 0.5) {
  levs <- unique(dat[[gradient_col]])
  dat$grad_num <- as.numeric(factor(dat[[gradient_col]], levels = levs))
}

# ---------- Presence grid ----------
pres_long <- dat %>%
  distinct(grad_num, !!sym(taxon_col)) %>%
  mutate(Presence = 1L) %>%
  complete(
    grad_num = sort(unique(grad_num)),
    !!sym(taxon_col),
    fill = list(Presence = 0L)
  ) %>%
  arrange(grad_num, !!sym(taxon_col))

A_wide <- pres_long %>%
  tidyr::pivot_wider(names_from = !!sym(taxon_col), values_from = Presence, values_fill = 0) %>%
  arrange(grad_num)

grad_vec <- A_wide$grad_num
A <- as.matrix(A_wide[, -1, drop = FALSE])
phyla <- colnames(A)
L <- nrow(A); P <- ncol(A)

# ---------- Monotonicity per phylum ----------
check_monotone <- function(v, tol = 0) {
  dv <- diff(v)
  inc <- all(dv >= -tol); dec <- all(dv <= tol)
  if (inc && !dec) return(c(TRUE, "non-decreasing"))
  if (dec && !inc) return(c(TRUE, "non-increasing"))
  if (inc &&  dec) return(c(TRUE, "constant"))
  c(FALSE, "non-monotonic")
}
mono_tbl <- lapply(seq_len(P), function(j){
  v  <- A[, j]
  md <- check_monotone(v)
  sp <- suppressWarnings(cor.test(rank(grad_vec), v, method = "spearman"))
  data.frame(
    Phylum        = phyla[j],
    monotone      = as.logical(md[1]),
    direction     = md[2],
    spearman_rho  = unname(sp$estimate),
    spearman_p    = sp$p.value,
    n_levels      = length(v),
    prevalence    = sum(v),
    fraction_pres = mean(v),
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()
write.csv(mono_tbl, "depth_monotonicity_by_phylum.csv", row.names = FALSE)

# ---------- Smoothing (Gaussian σ = 8) ----------
gaussian_kernel <- function(sigma_bins) {
  half <- max(1L, round(3 * sigma_bins))
  x <- -half:half
  k <- exp(-(x^2) / (2 * sigma_bins^2))
  k / sum(k)
}
smooth_reflect <- function(y, k) {
  pad <- floor((length(k) - 1) / 2)
  ypad <- c(rev(y[2:(pad+1)]), y, rev(y[(length(y)-pad):(length(y)-1)]))
  yy <- stats::filter(ypad, k, sides = 2)
  as.numeric(yy[(pad+1):(pad+length(y))])
}
k <- gaussian_kernel(sigma_bins)
Z <- apply(A, 2, smooth_reflect, k = k)

# Standardize columns; constant columns become all-zero vectors
Z <- scale(Z, center = TRUE, scale = TRUE)
Z[is.na(Z)] <- 0

# Identify constant / near-constant columns post-smoothing
col_sd <- apply(Z, 2, sd)
const_cols <- phyla[col_sd == 0]
if (length(const_cols)) {
  message("Constant (zero-variance) after smoothing: ",
          paste(const_cols, collapse = ", "))
}

# ---------- Safe correlation matrix (no NA) ----------
# Because Z columns are mean-zero and unit-sd (or zero if constant),
# R = Z'Z / L gives correlations; constant columns simply give 0 with others.
R <- (t(Z) %*% Z) / L
diag(R) <- 1
R <- as.matrix(R)
rownames(R) <- colnames(R) <- phyla

# ---------- Circular-shift permutation p-values ----------
perm_pvalue_corr <- function(x, y, L, n_perm) {
  shifts <- sample.int(L, n_perm, replace = TRUE) - 1L
  idx <- (rep(0:(L-1), times = n_perm) - rep(shifts, each = L)) %% L + 1L
  yperm <- matrix(y[idx], nrow = L, ncol = n_perm)
  r_perm <- colSums(x * yperm) / L
  r_obs  <- sum(x * y) / L
  p <- (1 + sum(abs(r_perm) >= abs(r_obs))) / (1 + n_perm)
  c(r_obs = r_obs, p = p)
}

pairs <- list(); cnt <- 1L
for (i in 1:(P-1)) {
  xi <- Z[, i]
  for (j in (i+1):P) {
    out <- perm_pvalue_corr(xi, Z[, j], L, n_perm)
    pairs[[cnt]] <- data.frame(
      Phylum1     = phyla[i],
      Phylum2     = phyla[j],
      curve_corr  = out["r_obs"],
      p_value     = out["p"],
      stringsAsFactors = FALSE
    ); cnt <- cnt + 1L
  }
}
pairs_df <- bind_rows(pairs)

# BH-FDR (optional)
o <- order(pairs_df$p_value)
r <- integer(length(o)); r[o] <- seq_along(o)
q <- pairs_df$p_value * nrow(pairs_df) / r
q[o] <- rev(cummin(rev(q[o])))
pairs_df$q_value <- q
write.csv(pairs_df, "depth_curve_pairs_sigma8_perm5000.csv", row.names = FALSE)

# ---------- Build heatmap data (same columns as your Spearman code) ----------
# Fill symmetric matrices from pairs and R (already safe)
R_mat <- R
P_mat <- matrix(0, P, P, dimnames = list(phyla, phyla))
apply(pairs_df, 1, function(row){
  i <- row[["Phylum1"]]; j <- row[["Phylum2"]]
  P_mat[i, j] <<- as.numeric(row[["p_value"]])
  P_mat[j, i] <<- as.numeric(row[["p_value"]])
})

cor_df  <- as.data.frame(as.table(R_mat))
pval_df <- as.data.frame(as.table(P_mat))
colnames(cor_df)  <- c("Phylum1","Phylum2","Spearman_r")  # reuse plot mapping
colnames(pval_df) <- c("Phylum1","Phylum2","p_value")

merged_df <- inner_join(cor_df, pval_df, by = c("Phylum1","Phylum2")) %>%
  filter(Phylum1 != Phylum2)

write.csv(merged_df, "CBB_depth_curvecorr_pvalue.csv", row.names = FALSE)

# ---------- Heatmap (same aesthetics) ----------
div_pal <- brewer.pal(n = 11, name = "RdBu")
heatmap_data <- merged_df %>%
  mutate(label = ifelse(p_value < alpha_lbl, round(Spearman_r, 2), NA))

depth_corr <- ggplot(heatmap_data, aes(x = Phylum1, y = Phylum2, fill = Spearman_r)) +
  geom_tile(color = "grey90", size = 0.3) +
  scale_fill_gradient2(
    low = div_pal[1], mid = div_pal[6], high = div_pal[11],
    midpoint = 0, limits = c(-1, 1), space = "Lab",
    na.value = "white",
    name = "Curve corr\n(σ=8)"
  ) +
  geom_text(aes(label = label), color = "black", size = 3.5, fontface = "bold", na.rm = TRUE) +
  coord_fixed() +
  labs(title = "", subtitle = "", x = NULL, y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x      = element_text(angle = 45, vjust = 1, hjust = 1, face = "bold"),
    axis.text.y      = element_text(face = "bold"),
    legend.title     = element_text(face = "bold", size = 12),
    legend.text      = element_text(size = 10),
    panel.grid       = element_blank(),
    panel.border     = element_rect(fill = NA, color = "grey50"),
    axis.ticks       = element_blank(),
    plot.margin      = margin(t = 20, r = 20, b = 20, l = 20)
  )

print(depth_corr)
ggsave("CBB_depth_curvecorr_heatmap_sigma8.png", depth_corr, width = 10, height = 8, dpi = 300)


# ============================
# Salinity monotonicity + smoothed-curve co-occurrence (σ = 8)
# Safe correlations (no NA) and same heatmap style
# ============================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(rlang)
  library(Hmisc); library(RColorBrewer)
})

# ---------- Parameters ----------
gradient_col <- "sample_salinity"
taxon_col    <- "GTDB_Tk_Phylum"
sigma_bins   <- 8
n_perm       <- 5000
alpha_lbl    <- 0.05
set.seed(7)

# ---------- Read & tidy ----------
dat <- read.csv("cbb_data.csv", stringsAsFactors = FALSE) %>%
  filter(.data[[gradient_col]] != "Unknown", !is.na(.data[[gradient_col]])) %>%
  mutate(!!taxon_col := trimws(as.character(.data[[taxon_col]])))

# Parse gradient to numeric (supports ranges like "0–200")
parse_num <- function(x) {
  x <- as.character(x)
  val <- suppressWarnings(as.numeric(x))
  if (!is.na(val)) return(val)
  s <- gsub(",", "", x)
  nums <- regmatches(s, gregexpr("\\d+\\.?\\d*", s))[[1]]
  if (length(nums) >= 2 && grepl("-|–|—| to ", s)) {
    return(mean(as.numeric(nums[1:2])))
  } else if (length(nums) >= 1) {
    return(as.numeric(nums[1]))
  }
  NA_real_
}
dat$grad_num <- vapply(dat[[gradient_col]], parse_num, numeric(1))
if (mean(is.na(dat$grad_num)) > 0.5) {
  levs <- unique(dat[[gradient_col]])
  dat$grad_num <- as.numeric(factor(dat[[gradient_col]], levels = levs))
}

# ---------- Presence grid ----------
pres_long <- dat %>%
  distinct(grad_num, !!sym(taxon_col)) %>%
  mutate(Presence = 1L) %>%
  complete(
    grad_num = sort(unique(grad_num)),
    !!sym(taxon_col),
    fill = list(Presence = 0L)
  ) %>%
  arrange(grad_num, !!sym(taxon_col))

A_wide <- pres_long %>%
  tidyr::pivot_wider(names_from = !!sym(taxon_col), values_from = Presence, values_fill = 0) %>%
  arrange(grad_num)

grad_vec <- A_wide$grad_num
A <- as.matrix(A_wide[, -1, drop = FALSE])
phyla <- colnames(A)
L <- nrow(A); P <- ncol(A)

# ---------- Monotonicity per phylum ----------
check_monotone <- function(v, tol = 0) {
  dv <- diff(v)
  inc <- all(dv >= -tol); dec <- all(dv <= tol)
  if (inc && !dec) return(c(TRUE, "non-decreasing"))
  if (dec && !inc) return(c(TRUE, "non-increasing"))
  if (inc &&  dec) return(c(TRUE, "constant"))
  c(FALSE, "non-monotonic")
}
mono_tbl <- lapply(seq_len(P), function(j){
  v  <- A[, j]
  md <- check_monotone(v)
  sp <- suppressWarnings(cor.test(rank(grad_vec), v, method = "spearman"))
  data.frame(
    Phylum        = phyla[j],
    monotone      = as.logical(md[1]),
    direction     = md[2],
    spearman_rho  = unname(sp$estimate),
    spearman_p    = sp$p.value,
    n_levels      = length(v),
    prevalence    = sum(v),
    fraction_pres = mean(v),
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()
write.csv(mono_tbl, "salinity_monotonicity_by_phylum.csv", row.names = FALSE)

# ---------- Smoothing (Gaussian σ = 8) ----------
gaussian_kernel <- function(sigma_bins) {
  half <- max(1L, round(3 * sigma_bins))
  x <- -half:half
  k <- exp(-(x^2) / (2 * sigma_bins^2))
  k / sum(k)
}
smooth_reflect <- function(y, k) {
  pad <- floor((length(k) - 1) / 2)
  ypad <- c(rev(y[2:(pad+1)]), y, rev(y[(length(y)-pad):(length(y)-1)]))
  yy <- stats::filter(ypad, k, sides = 2)
  as.numeric(yy[(pad+1):(pad+length(y))])
}
k <- gaussian_kernel(sigma_bins)
Z <- apply(A, 2, smooth_reflect, k = k)

# Standardize columns; constant columns become all-zero vectors
Z <- scale(Z, center = TRUE, scale = TRUE)
Z[is.na(Z)] <- 0

# Identify constant / near-constant columns post-smoothing
col_sd <- apply(Z, 2, sd)
const_cols <- phyla[col_sd == 0]
if (length(const_cols)) {
  message("Constant (zero-variance) after smoothing: ",
          paste(const_cols, collapse = ", "))
}

# ---------- Safe correlation matrix (no NA) ----------
# Because Z columns are mean-zero and unit-sd (or zero if constant),
# R = Z'Z / L gives correlations; constant columns simply give 0 with others.
R <- (t(Z) %*% Z) / L
diag(R) <- 1
R <- as.matrix(R)
rownames(R) <- colnames(R) <- phyla

# ---------- Circular-shift permutation p-values ----------
perm_pvalue_corr <- function(x, y, L, n_perm) {
  shifts <- sample.int(L, n_perm, replace = TRUE) - 1L
  idx <- (rep(0:(L-1), times = n_perm) - rep(shifts, each = L)) %% L + 1L
  yperm <- matrix(y[idx], nrow = L, ncol = n_perm)
  r_perm <- colSums(x * yperm) / L
  r_obs  <- sum(x * y) / L
  p <- (1 + sum(abs(r_perm) >= abs(r_obs))) / (1 + n_perm)
  c(r_obs = r_obs, p = p)
}

pairs <- list(); cnt <- 1L
for (i in 1:(P-1)) {
  xi <- Z[, i]
  for (j in (i+1):P) {
    out <- perm_pvalue_corr(xi, Z[, j], L, n_perm)
    pairs[[cnt]] <- data.frame(
      Phylum1     = phyla[i],
      Phylum2     = phyla[j],
      curve_corr  = out["r_obs"],
      p_value     = out["p"],
      stringsAsFactors = FALSE
    ); cnt <- cnt + 1L
  }
}
pairs_df <- bind_rows(pairs)

# BH-FDR (optional)
o <- order(pairs_df$p_value)
r <- integer(length(o)); r[o] <- seq_along(o)
q <- pairs_df$p_value * nrow(pairs_df) / r
q[o] <- rev(cummin(rev(q[o])))
pairs_df$q_value <- q
write.csv(pairs_df, "salinity_curve_pairs_sigma8_perm5000.csv", row.names = FALSE)

# ---------- Build heatmap data (same columns as your Spearman code) ----------
# Fill symmetric matrices from pairs and R (already safe)
R_mat <- R
P_mat <- matrix(0, P, P, dimnames = list(phyla, phyla))
apply(pairs_df, 1, function(row){
  i <- row[["Phylum1"]]; j <- row[["Phylum2"]]
  P_mat[i, j] <<- as.numeric(row[["p_value"]])
  P_mat[j, i] <<- as.numeric(row[["p_value"]])
})

cor_df  <- as.data.frame(as.table(R_mat))
pval_df <- as.data.frame(as.table(P_mat))
colnames(cor_df)  <- c("Phylum1","Phylum2","Spearman_r")  # reuse plot mapping
colnames(pval_df) <- c("Phylum1","Phylum2","p_value")

merged_df <- inner_join(cor_df, pval_df, by = c("Phylum1","Phylum2")) %>%
  filter(Phylum1 != Phylum2)

write.csv(merged_df, "CBB_salinity_curvecorr_pvalue.csv", row.names = FALSE)

# ---------- Heatmap (same aesthetics) ----------
div_pal <- brewer.pal(n = 11, name = "RdBu")
heatmap_data <- merged_df %>%
  mutate(label = ifelse(p_value < alpha_lbl, round(Spearman_r, 2), NA))

salinity_corr <- ggplot(heatmap_data, aes(x = Phylum1, y = Phylum2, fill = Spearman_r)) +
  geom_tile(color = "grey90", size = 0.3) +
  scale_fill_gradient2(
    low = div_pal[1], mid = div_pal[6], high = div_pal[11],
    midpoint = 0, limits = c(-1, 1), space = "Lab",
    na.value = "white",
    name = "Curve corr\n(σ=8)"
  ) +
  geom_text(aes(label = label), color = "black", size = 3.5, fontface = "bold", na.rm = TRUE) +
  coord_fixed() +
  labs(title = "", subtitle = "", x = NULL, y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x      = element_text(angle = 45, vjust = 1, hjust = 1, face = "bold"),
    axis.text.y      = element_text(face = "bold"),
    legend.title     = element_text(face = "bold", size = 12),
    legend.text      = element_text(size = 10),
    panel.grid       = element_blank(),
    panel.border     = element_rect(fill = NA, color = "grey50"),
    axis.ticks       = element_blank(),
    plot.margin      = margin(t = 20, r = 20, b = 20, l = 20)
  )

print(salinity_corr)
ggsave("CBB_salinity_curvecorr_heatmap_sigma8.png", salinity_corr, width = 10, height = 8, dpi = 300)

######################
###merge correlation plots

# install.packages(c("gridExtra","cowplot"))
library(ggplot2)
library(gridExtra)
library(cowplot)   # for get_legend()
library(grid)      # for unit()

svg("D:/Thesis/work_package_3/Data_analysis/01_figures/08_correlation_analysis/new_correlation_analysis_04102025/merged_cbb_correlations.svg",width = 24, height = 10)

# 1) Common theme: put legend on the right, vertical orientation
common_theme <- theme_minimal(base_size = 20) +
  theme(
    axis.text.x     = element_text(angle = 45, vjust = 1, hjust = 1, size = 20),
    axis.text.y     = element_text(size = 20),
    legend.position = "right",                # move legend to right
    legend.direction= "vertical",             # force a vertical legend
    legend.title    = element_text(face = "bold", size = 20),
    legend.text     = element_text(size = 20),
    plot.margin     = margin(6, 5, 6, 5)
  )

# 2) Tag theme (unchanged)
tag_theme <- theme(
  plot.tag.position = c(0, 1),
  plot.tag          = element_text(face = "bold", size = 19)
)

# 3) Build each panel with tags + the common fill scale
x2 <- depth_corr + common_theme + labs(tag="(b)") + tag_theme
y2 <- salinity_corr   + common_theme + labs(tag="(c)") + tag_theme
z2 <- temp_corr  + common_theme + labs(tag="(d)") + tag_theme

# 4) Grab the vertical legend from one plot
shared_leg <- get_legend(x2)

# 5) Strip legends off the individual panels
x2 <- x2 + theme(legend.position = "none")
y2 <- y2 + theme(legend.position = "none")
z2 <- z2 + theme(legend.position = "none")

# 6) Arrange your three tagged panels in a row
panel_row <- arrangeGrob(x2, y2, z2, ncol = 3)

# 7) Put panels and legend side by side
grid.arrange(
  panel_row,    # first column: the three plots
  shared_leg,   # second column: the vertical legend
  ncol   = 2,
  widths = c(5, 1)   # adjust as needed: wider left, narrow legend
)

dev.off()


#################################################################

#set working directory
setwd("D:/Thesis/work_package_3/Data_analysis/01_figures/08_correlation_analysis/new_correlation_analysis_04102025/correlation_all")


##################################
### Co-occurenceof phylum by depth
##################################

# ===== FULL smoothed-curve correlations (all phyla) + TOP-50 heatmap =====
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(pheatmap)
})

# ---------------- USER SETTINGS ----------------
infile          <- "merged_data.csv"
gradient_col    <- "sample_depth.x"      # or "sample_temp"
taxon_col       <- "GTDB_Tk_Phylum"

sigma_bins      <- 8                     # Gaussian bandwidth (bins)
n_perm          <- 5000                  # circular-shift permutations (raise for final)
alpha_FDR       <- 0.05                  # FDR threshold for showing tiles
TOP_N           <- 50                    # only for PLOTTING
show_upper_only <- TRUE                  # plot only upper triangle?
fallback_edges  <- 250                   # used only if nothing passes FDR

# Deep blue ↔ white ↔ deep red (scico 'vik' if available)
make_diverging <- function(n = 257) {
  if (requireNamespace("scico", quietly = TRUE))
    return(scico::scico(n, palette = "vik"))
  if (requireNamespace("RColorBrewer", quietly = TRUE))
    return(colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(n))
  colorRampPalette(c("#053061", "#2166AC", "#f7f7f7", "#B2182B", "#67001F"))(n)
}

# ---------------- READ & PRESENCE MATRIX (ALL PHYLA) ----------------
dat <- read.csv(infile, stringsAsFactors = FALSE) %>%
  filter(.data[[gradient_col]] != "Unknown", !is.na(.data[[gradient_col]])) %>%
  mutate(
    !!taxon_col := trimws(as.character(.data[[taxon_col]])),
    grad_str    = as.character(.data[[gradient_col]])
  )

# keep gradient bins in their original order
grad_levels <- unique(dat$grad_str)

presence_matrix <- dat %>%
  select(grad_str, !!sym(taxon_col)) %>%
  distinct() %>%
  mutate(Presence = 1L) %>%
  complete(grad_str = grad_levels, !!sym(taxon_col), fill = list(Presence = 0L)) %>%
  arrange(match(grad_str, grad_levels))

A_all_df <- presence_matrix %>%
  pivot_wider(names_from = !!sym(taxon_col), values_from = Presence, values_fill = 0) %>%
  arrange(match(grad_str, grad_levels))

grad_str <- A_all_df$grad_str           # keep labels for reference
A_all    <- as.matrix(A_all_df[, -1, drop = FALSE])  # L x P
L        <- nrow(A_all); P_all <- ncol(A_all)
phyla_all <- colnames(A_all)

# numeric order (1..L) for monotonicity/Spearman check
grad_num <- seq_len(L)

# ---------------- MONOTONICITY (ALL PHYLA) ----------------
check_monotone <- function(v, tol = 0) {
  dv <- diff(v)
  inc <- all(dv >= -tol); dec <- all(dv <= tol)
  if (inc && !dec) return(c(TRUE, "non-decreasing"))
  if (dec && !inc) return(c(TRUE, "non-increasing"))
  if (inc &&  dec) return(c(TRUE, "constant"))
  c(FALSE, "non-monotonic")
}

mono_tbl <- lapply(seq_len(P_all), function(j){
  v  <- A_all[, j]
  md <- check_monotone(v)
  sp <- suppressWarnings(cor.test(grad_num, v, method = "spearman"))
  data.frame(
    Phylum        = phyla_all[j],
    monotone      = as.logical(md[1]),
    direction     = md[2],
    spearman_rho  = unname(sp$estimate),
    spearman_p    = sp$p.value,
    n_bins        = length(v),
    prevalence    = sum(v),
    fraction_pres = mean(v),
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()
write.csv(mono_tbl, "depth_monotonicity_by_phylum_ALL.csv", row.names = FALSE)

# ---------------- SMOOTH & STANDARDIZE (ALL) ----------------
gaussian_kernel <- function(sigma_bins) {
  half <- max(1L, round(3 * sigma_bins))
  x <- -half:half
  k <- exp(-(x^2) / (2 * sigma_bins^2))
  k / sum(k)
}
smooth_reflect <- function(y, k) {
  pad <- floor((length(k) - 1) / 2)
  ypad <- c(rev(y[2:(pad+1)]), y, rev(y[(length(y)-pad):(length(y)-1)]))
  yy <- stats::filter(ypad, k, sides = 2)
  as.numeric(yy[(pad+1):(pad+length(y))])
}
k  <- gaussian_kernel(sigma_bins)
Z_all <- apply(A_all, 2, smooth_reflect, k = k)   # L x P_all
Z_all <- scale(Z_all, center = TRUE, scale = TRUE)
Z_all[is.na(Z_all)] <- 0

# Safe correlation matrix for ALL phyla (no NA): R = Z'Z / L
R_full <- (t(Z_all) %*% Z_all) / L
R_full <- as.matrix(R_full); diag(R_full) <- 1
rownames(R_full) <- colnames(R_full) <- phyla_all

# ---------------- PERMUTATION p-VALUES (ALL pairs) ----------------
perm_pvalue_corr <- function(x, y, L, n_perm) {
  shifts <- sample.int(L, n_perm, replace = TRUE) - 1L
  idx <- (rep(0:(L-1), times = n_perm) - rep(shifts, each = L)) %% L + 1L
  yperm <- matrix(y[idx], nrow = L, ncol = n_perm)
  r_perm <- colMeans(x * yperm)
  r_obs  <- mean(x * y)
  p <- (1 + sum(abs(r_perm) >= abs(r_obs))) / (1 + n_perm)
  c(r_obs = r_obs, p = p)
}

pairs_all <- vector("list", P_all * (P_all - 1) / 2); ctr <- 1L
for (i in 1:(P_all-1)) {
  xi <- Z_all[, i]
  for (j in (i+1):P_all) {
    out <- perm_pvalue_corr(xi, Z_all[, j], L, n_perm)
    pairs_all[[ctr]] <- data.frame(
      Phylum1    = phyla_all[i],
      Phylum2    = phyla_all[j],
      curve_corr = out["r_obs"],
      p_value    = out["p"],
      Sigma_bins = sigma_bins,
      Permutations = n_perm,
      stringsAsFactors = FALSE
    ); ctr <- ctr + 1L
  }
}
pairs_all <- bind_rows(pairs_all)

# BH-FDR across ALL pairs
o <- order(pairs_all$p_value)
r <- integer(length(o)); r[o] <- seq_along(o)
q <- pairs_all$p_value * nrow(pairs_all) / r
q[o] <- rev(cummin(rev(q[o])))
pairs_all$q_value <- q

# ---- SAVE FULL TABLES (ALL PHYLA) ----
write.csv(pairs_all, "depth_curve_correlations_pairs_ALL.csv", row.names = FALSE)
write.csv(R_full,   "depth_curve_corr_matrix_ALL.csv", row.names = TRUE)

# Build symmetric q-matrix (for later subsetting)
Q_full <- matrix(NA_real_, P_all, P_all, dimnames = list(phyla_all, phyla_all))
apply(pairs_all, 1, function(row){
  i <- row[["Phylum1"]]; j <- row[["Phylum2"]]
  Q_full[i, j] <<- Q_full[j, i] <<- as.numeric(row[["q_value"]])
})
diag(Q_full) <- NA_real_
write.csv(Q_full, "depth_curve_q_matrix_ALL.csv", row.names = TRUE)

# ---------------- PLOT: TOP 50 only ----------------
# Choose top 50 by occupancy, tie-break by variance (computed on raw presence)
occ_all  <- colSums(A_all, na.rm = TRUE)
varv_all <- apply(A_all, 2, var)
ord_all  <- order(-occ_all, -varv_all)
keep     <- phyla_all[ord_all][seq_len(min(TOP_N, length(ord_all)))]

R_top <- R_full[keep, keep, drop = FALSE]
Q_top <- Q_full[keep, keep, drop = FALSE]

# Mask non-significant tiles (FDR)
R_sig <- R_top
R_sig[is.na(Q_top) | Q_top >= alpha_FDR] <- NA_real_
diag(R_sig) <- NA_real_

# Fallback: if everything is NA after masking, show top |r| edges
if (all(is.na(R_sig))) {
  message("No cells passed FDR<", alpha_FDR,
          ". Showing top ", fallback_edges, " edges by |r| as fallback.")
  Rt <- R_top; diag(Rt) <- NA
  ord_edges <- order(abs(Rt), decreasing = TRUE, na.last = NA)
  show_idx <- ord_edges[seq_len(min(fallback_edges, length(ord_edges)))]
  R_sig[,] <- NA_real_
  R_sig[show_idx] <- R_top[show_idx]
}

# Drop rows/cols with no visible links and cluster
keep_idx <- rowSums(!is.na(R_sig)) > 0
R_sig <- R_sig[keep_idx, keep_idx, drop = FALSE]

R_for_order <- R_sig; R_for_order[is.na(R_for_order)] <- 0
hc <- hclust(as.dist(1 - R_for_order), method = "average")

R_plot <- R_sig
if (show_upper_only) R_plot[lower.tri(R_plot, diag = TRUE)] <- NA_real_

pal  <- make_diverging(257)
brks <- seq(-1, 1, length.out = length(pal))

#####################
### Figure S7

# SVG
svg("heatmap_phylum_depth_TOP50_sigma8_FDR.svg", width = 10, height = 10)
pheatmap(
  mat = R_plot,
  color = pal, breaks = brks, na_col = "grey95",
  cluster_rows = hc, cluster_cols = hc,
  border_color = NA, legend = TRUE,
  show_rownames = TRUE, show_colnames = TRUE,
  fontsize_row = 11, fontsize_col = 11, angle_col = 90,
  treeheight_row = 50, treeheight_col = 50, fontsize = 12
)
dev.off()

# High-res PNG
png("heatmap_phylum_depth_TOP50_sigma8_FDR.png", width = 5200, height = 5000, res = 550)
pheatmap(
  mat = R_plot,
  color = pal, breaks = brks, na_col = "grey95",
  cluster_rows = hc, cluster_cols = hc,
  border_color = NA, legend = TRUE,
  show_rownames = TRUE, show_colnames = TRUE,
  fontsize_row = 11, fontsize_col = 11, angle_col = 90,
  treeheight_row = 50, treeheight_col = 50, fontsize = 12
)
dev.off()


##################################
### Co-occurenceof phylum by salinity
##################################

# ===== FULL smoothed-curve correlations (all phyla) + TOP-50 heatmap =====
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(pheatmap)
})

# ---------------- USER SETTINGS ----------------
infile          <- "merged_data.csv"
gradient_col    <- "sample_salinity.x"      # or "sample_temp"
taxon_col       <- "GTDB_Tk_Phylum"

sigma_bins      <- 8                     # Gaussian bandwidth (bins)
n_perm          <- 5000                  # circular-shift permutations (raise for final)
alpha_FDR       <- 0.05                  # FDR threshold for showing tiles
TOP_N           <- 50                    # only for PLOTTING
show_upper_only <- TRUE                  # plot only upper triangle?
fallback_edges  <- 250                   # used only if nothing passes FDR

# Deep blue ↔ white ↔ deep red (scico 'vik' if available)
make_diverging <- function(n = 257) {
  if (requireNamespace("scico", quietly = TRUE))
    return(scico::scico(n, palette = "vik"))
  if (requireNamespace("RColorBrewer", quietly = TRUE))
    return(colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(n))
  colorRampPalette(c("#053061", "#2166AC", "#f7f7f7", "#B2182B", "#67001F"))(n)
}

# ---------------- READ & PRESENCE MATRIX (ALL PHYLA) ----------------
dat <- read.csv(infile, stringsAsFactors = FALSE) %>%
  filter(.data[[gradient_col]] != "Unknown", !is.na(.data[[gradient_col]])) %>%
  mutate(
    !!taxon_col := trimws(as.character(.data[[taxon_col]])),
    grad_str    = as.character(.data[[gradient_col]])
  )

# keep gradient bins in their original order
grad_levels <- unique(dat$grad_str)

presence_matrix <- dat %>%
  select(grad_str, !!sym(taxon_col)) %>%
  distinct() %>%
  mutate(Presence = 1L) %>%
  complete(grad_str = grad_levels, !!sym(taxon_col), fill = list(Presence = 0L)) %>%
  arrange(match(grad_str, grad_levels))

A_all_df <- presence_matrix %>%
  pivot_wider(names_from = !!sym(taxon_col), values_from = Presence, values_fill = 0) %>%
  arrange(match(grad_str, grad_levels))

grad_str <- A_all_df$grad_str           # keep labels for reference
A_all    <- as.matrix(A_all_df[, -1, drop = FALSE])  # L x P
L        <- nrow(A_all); P_all <- ncol(A_all)
phyla_all <- colnames(A_all)

# numeric order (1..L) for monotonicity/Spearman check
grad_num <- seq_len(L)

# ---------------- MONOTONICITY (ALL PHYLA) ----------------
check_monotone <- function(v, tol = 0) {
  dv <- diff(v)
  inc <- all(dv >= -tol); dec <- all(dv <= tol)
  if (inc && !dec) return(c(TRUE, "non-decreasing"))
  if (dec && !inc) return(c(TRUE, "non-increasing"))
  if (inc &&  dec) return(c(TRUE, "constant"))
  c(FALSE, "non-monotonic")
}

mono_tbl <- lapply(seq_len(P_all), function(j){
  v  <- A_all[, j]
  md <- check_monotone(v)
  sp <- suppressWarnings(cor.test(grad_num, v, method = "spearman"))
  data.frame(
    Phylum        = phyla_all[j],
    monotone      = as.logical(md[1]),
    direction     = md[2],
    spearman_rho  = unname(sp$estimate),
    spearman_p    = sp$p.value,
    n_bins        = length(v),
    prevalence    = sum(v),
    fraction_pres = mean(v),
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()
write.csv(mono_tbl, "depth_monotonicity_by_phylum_ALL.csv", row.names = FALSE)

# ---------------- SMOOTH & STANDARDIZE (ALL) ----------------
gaussian_kernel <- function(sigma_bins) {
  half <- max(1L, round(3 * sigma_bins))
  x <- -half:half
  k <- exp(-(x^2) / (2 * sigma_bins^2))
  k / sum(k)
}
smooth_reflect <- function(y, k) {
  pad <- floor((length(k) - 1) / 2)
  ypad <- c(rev(y[2:(pad+1)]), y, rev(y[(length(y)-pad):(length(y)-1)]))
  yy <- stats::filter(ypad, k, sides = 2)
  as.numeric(yy[(pad+1):(pad+length(y))])
}
k  <- gaussian_kernel(sigma_bins)
Z_all <- apply(A_all, 2, smooth_reflect, k = k)   # L x P_all
Z_all <- scale(Z_all, center = TRUE, scale = TRUE)
Z_all[is.na(Z_all)] <- 0

# Safe correlation matrix for ALL phyla (no NA): R = Z'Z / L
R_full <- (t(Z_all) %*% Z_all) / L
R_full <- as.matrix(R_full); diag(R_full) <- 1
rownames(R_full) <- colnames(R_full) <- phyla_all

# ---------------- PERMUTATION p-VALUES (ALL pairs) ----------------
perm_pvalue_corr <- function(x, y, L, n_perm) {
  shifts <- sample.int(L, n_perm, replace = TRUE) - 1L
  idx <- (rep(0:(L-1), times = n_perm) - rep(shifts, each = L)) %% L + 1L
  yperm <- matrix(y[idx], nrow = L, ncol = n_perm)
  r_perm <- colMeans(x * yperm)
  r_obs  <- mean(x * y)
  p <- (1 + sum(abs(r_perm) >= abs(r_obs))) / (1 + n_perm)
  c(r_obs = r_obs, p = p)
}

pairs_all <- vector("list", P_all * (P_all - 1) / 2); ctr <- 1L
for (i in 1:(P_all-1)) {
  xi <- Z_all[, i]
  for (j in (i+1):P_all) {
    out <- perm_pvalue_corr(xi, Z_all[, j], L, n_perm)
    pairs_all[[ctr]] <- data.frame(
      Phylum1    = phyla_all[i],
      Phylum2    = phyla_all[j],
      curve_corr = out["r_obs"],
      p_value    = out["p"],
      Sigma_bins = sigma_bins,
      Permutations = n_perm,
      stringsAsFactors = FALSE
    ); ctr <- ctr + 1L
  }
}
pairs_all <- bind_rows(pairs_all)

# BH-FDR across ALL pairs
o <- order(pairs_all$p_value)
r <- integer(length(o)); r[o] <- seq_along(o)
q <- pairs_all$p_value * nrow(pairs_all) / r
q[o] <- rev(cummin(rev(q[o])))
pairs_all$q_value <- q

# ---- SAVE FULL TABLES (ALL PHYLA) ----
write.csv(pairs_all, "salinity_curve_correlations_pairs_ALL.csv", row.names = FALSE)
write.csv(R_full,   "salinity_curve_corr_matrix_ALL.csv", row.names = TRUE)

# Build symmetric q-matrix (for later subsetting)
Q_full <- matrix(NA_real_, P_all, P_all, dimnames = list(phyla_all, phyla_all))
apply(pairs_all, 1, function(row){
  i <- row[["Phylum1"]]; j <- row[["Phylum2"]]
  Q_full[i, j] <<- Q_full[j, i] <<- as.numeric(row[["q_value"]])
})
diag(Q_full) <- NA_real_
write.csv(Q_full, "salinity_curve_q_matrix_ALL.csv", row.names = TRUE)

# ---------------- PLOT: TOP 50 only ----------------
# Choose top 50 by occupancy, tie-break by variance (computed on raw presence)
occ_all  <- colSums(A_all, na.rm = TRUE)
varv_all <- apply(A_all, 2, var)
ord_all  <- order(-occ_all, -varv_all)
keep     <- phyla_all[ord_all][seq_len(min(TOP_N, length(ord_all)))]

R_top <- R_full[keep, keep, drop = FALSE]
Q_top <- Q_full[keep, keep, drop = FALSE]

# Mask non-significant tiles (FDR)
R_sig <- R_top
R_sig[is.na(Q_top) | Q_top >= alpha_FDR] <- NA_real_
diag(R_sig) <- NA_real_

# Fallback: if everything is NA after masking, show top |r| edges
if (all(is.na(R_sig))) {
  message("No cells passed FDR<", alpha_FDR,
          ". Showing top ", fallback_edges, " edges by |r| as fallback.")
  Rt <- R_top; diag(Rt) <- NA
  ord_edges <- order(abs(Rt), decreasing = TRUE, na.last = NA)
  show_idx <- ord_edges[seq_len(min(fallback_edges, length(ord_edges)))]
  R_sig[,] <- NA_real_
  R_sig[show_idx] <- R_top[show_idx]
}

# Drop rows/cols with no visible links and cluster
keep_idx <- rowSums(!is.na(R_sig)) > 0
R_sig <- R_sig[keep_idx, keep_idx, drop = FALSE]

R_for_order <- R_sig; R_for_order[is.na(R_for_order)] <- 0
hc <- hclust(as.dist(1 - R_for_order), method = "average")

R_plot <- R_sig
if (show_upper_only) R_plot[lower.tri(R_plot, diag = TRUE)] <- NA_real_

pal  <- make_diverging(257)
brks <- seq(-1, 1, length.out = length(pal))

#####################
### Figure S8

# SVG
svg("heatmap_phylum_salinity_TOP50_sigma8_FDR.svg", width = 10, height = 10)
pheatmap(
  mat = R_plot,
  color = pal, breaks = brks, na_col = "grey95",
  cluster_rows = hc, cluster_cols = hc,
  border_color = NA, legend = TRUE,
  show_rownames = TRUE, show_colnames = TRUE,
  fontsize_row = 11, fontsize_col = 11, angle_col = 90,
  treeheight_row = 50, treeheight_col = 50, fontsize = 12
)
dev.off()

# High-res PNG
png("heatmap_phylum_salinity_TOP50_sigma8_FDR.png", width = 5200, height = 5000, res = 550)
pheatmap(
  mat = R_plot,
  color = pal, breaks = brks, na_col = "grey95",
  cluster_rows = hc, cluster_cols = hc,
  border_color = NA, legend = TRUE,
  show_rownames = TRUE, show_colnames = TRUE,
  fontsize_row = 11, fontsize_col = 11, angle_col = 90,
  treeheight_row = 50, treeheight_col = 50, fontsize = 12
)
dev.off()


##################################
### Co-occurenceof phylum by temperature
##################################

# ===== FULL smoothed-curve correlations (all phyla) + TOP-50 heatmap =====
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(pheatmap)
})

# ---------------- USER SETTINGS ----------------
infile          <- "merged_data.csv"
gradient_col    <- "sample_temp.x"      # or "sample_temp"
taxon_col       <- "GTDB_Tk_Phylum"

sigma_bins      <- 8                     # Gaussian bandwidth (bins)
n_perm          <- 5000                  # circular-shift permutations (raise for final)
alpha_FDR       <- 0.05                  # FDR threshold for showing tiles
TOP_N           <- 50                    # only for PLOTTING
show_upper_only <- TRUE                  # plot only upper triangle?
fallback_edges  <- 250                   # used only if nothing passes FDR

# Deep blue ↔ white ↔ deep red (scico 'vik' if available)
make_diverging <- function(n = 257) {
  if (requireNamespace("scico", quietly = TRUE))
    return(scico::scico(n, palette = "vik"))
  if (requireNamespace("RColorBrewer", quietly = TRUE))
    return(colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(n))
  colorRampPalette(c("#053061", "#2166AC", "#f7f7f7", "#B2182B", "#67001F"))(n)
}

# ---------------- READ & PRESENCE MATRIX (ALL PHYLA) ----------------
dat <- read.csv(infile, stringsAsFactors = FALSE) %>%
  filter(.data[[gradient_col]] != "Unknown", !is.na(.data[[gradient_col]])) %>%
  mutate(
    !!taxon_col := trimws(as.character(.data[[taxon_col]])),
    grad_str    = as.character(.data[[gradient_col]])
  )

# keep gradient bins in their original order
grad_levels <- unique(dat$grad_str)

presence_matrix <- dat %>%
  select(grad_str, !!sym(taxon_col)) %>%
  distinct() %>%
  mutate(Presence = 1L) %>%
  complete(grad_str = grad_levels, !!sym(taxon_col), fill = list(Presence = 0L)) %>%
  arrange(match(grad_str, grad_levels))

A_all_df <- presence_matrix %>%
  pivot_wider(names_from = !!sym(taxon_col), values_from = Presence, values_fill = 0) %>%
  arrange(match(grad_str, grad_levels))

grad_str <- A_all_df$grad_str           # keep labels for reference
A_all    <- as.matrix(A_all_df[, -1, drop = FALSE])  # L x P
L        <- nrow(A_all); P_all <- ncol(A_all)
phyla_all <- colnames(A_all)

# numeric order (1..L) for monotonicity/Spearman check
grad_num <- seq_len(L)

# ---------------- MONOTONICITY (ALL PHYLA) ----------------
check_monotone <- function(v, tol = 0) {
  dv <- diff(v)
  inc <- all(dv >= -tol); dec <- all(dv <= tol)
  if (inc && !dec) return(c(TRUE, "non-decreasing"))
  if (dec && !inc) return(c(TRUE, "non-increasing"))
  if (inc &&  dec) return(c(TRUE, "constant"))
  c(FALSE, "non-monotonic")
}

mono_tbl <- lapply(seq_len(P_all), function(j){
  v  <- A_all[, j]
  md <- check_monotone(v)
  sp <- suppressWarnings(cor.test(grad_num, v, method = "spearman"))
  data.frame(
    Phylum        = phyla_all[j],
    monotone      = as.logical(md[1]),
    direction     = md[2],
    spearman_rho  = unname(sp$estimate),
    spearman_p    = sp$p.value,
    n_bins        = length(v),
    prevalence    = sum(v),
    fraction_pres = mean(v),
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()
write.csv(mono_tbl, "temp_monotonicity_by_phylum_ALL.csv", row.names = FALSE)

# ---------------- SMOOTH & STANDARDIZE (ALL) ----------------
gaussian_kernel <- function(sigma_bins) {
  half <- max(1L, round(3 * sigma_bins))
  x <- -half:half
  k <- exp(-(x^2) / (2 * sigma_bins^2))
  k / sum(k)
}
smooth_reflect <- function(y, k) {
  pad <- floor((length(k) - 1) / 2)
  ypad <- c(rev(y[2:(pad+1)]), y, rev(y[(length(y)-pad):(length(y)-1)]))
  yy <- stats::filter(ypad, k, sides = 2)
  as.numeric(yy[(pad+1):(pad+length(y))])
}
k  <- gaussian_kernel(sigma_bins)
Z_all <- apply(A_all, 2, smooth_reflect, k = k)   # L x P_all
Z_all <- scale(Z_all, center = TRUE, scale = TRUE)
Z_all[is.na(Z_all)] <- 0

# Safe correlation matrix for ALL phyla (no NA): R = Z'Z / L
R_full <- (t(Z_all) %*% Z_all) / L
R_full <- as.matrix(R_full); diag(R_full) <- 1
rownames(R_full) <- colnames(R_full) <- phyla_all

# ---------------- PERMUTATION p-VALUES (ALL pairs) ----------------
perm_pvalue_corr <- function(x, y, L, n_perm) {
  shifts <- sample.int(L, n_perm, replace = TRUE) - 1L
  idx <- (rep(0:(L-1), times = n_perm) - rep(shifts, each = L)) %% L + 1L
  yperm <- matrix(y[idx], nrow = L, ncol = n_perm)
  r_perm <- colMeans(x * yperm)
  r_obs  <- mean(x * y)
  p <- (1 + sum(abs(r_perm) >= abs(r_obs))) / (1 + n_perm)
  c(r_obs = r_obs, p = p)
}

pairs_all <- vector("list", P_all * (P_all - 1) / 2); ctr <- 1L
for (i in 1:(P_all-1)) {
  xi <- Z_all[, i]
  for (j in (i+1):P_all) {
    out <- perm_pvalue_corr(xi, Z_all[, j], L, n_perm)
    pairs_all[[ctr]] <- data.frame(
      Phylum1    = phyla_all[i],
      Phylum2    = phyla_all[j],
      curve_corr = out["r_obs"],
      p_value    = out["p"],
      Sigma_bins = sigma_bins,
      Permutations = n_perm,
      stringsAsFactors = FALSE
    ); ctr <- ctr + 1L
  }
}
pairs_all <- bind_rows(pairs_all)

# BH-FDR across ALL pairs
o <- order(pairs_all$p_value)
r <- integer(length(o)); r[o] <- seq_along(o)
q <- pairs_all$p_value * nrow(pairs_all) / r
q[o] <- rev(cummin(rev(q[o])))
pairs_all$q_value <- q

# ---- SAVE FULL TABLES (ALL PHYLA) ----
write.csv(pairs_all, "temp_curve_correlations_pairs_ALL.csv", row.names = FALSE)
write.csv(R_full,   "temp_curve_corr_matrix_ALL.csv", row.names = TRUE)

# Build symmetric q-matrix (for later subsetting)
Q_full <- matrix(NA_real_, P_all, P_all, dimnames = list(phyla_all, phyla_all))
apply(pairs_all, 1, function(row){
  i <- row[["Phylum1"]]; j <- row[["Phylum2"]]
  Q_full[i, j] <<- Q_full[j, i] <<- as.numeric(row[["q_value"]])
})
diag(Q_full) <- NA_real_
write.csv(Q_full, "temp_curve_q_matrix_ALL.csv", row.names = TRUE)

# ---------------- PLOT: TOP 50 only ----------------
# Choose top 50 by occupancy, tie-break by variance (computed on raw presence)
occ_all  <- colSums(A_all, na.rm = TRUE)
varv_all <- apply(A_all, 2, var)
ord_all  <- order(-occ_all, -varv_all)
keep     <- phyla_all[ord_all][seq_len(min(TOP_N, length(ord_all)))]

R_top <- R_full[keep, keep, drop = FALSE]
Q_top <- Q_full[keep, keep, drop = FALSE]

# Mask non-significant tiles (FDR)
R_sig <- R_top
R_sig[is.na(Q_top) | Q_top >= alpha_FDR] <- NA_real_
diag(R_sig) <- NA_real_

# Fallback: if everything is NA after masking, show top |r| edges
if (all(is.na(R_sig))) {
  message("No cells passed FDR<", alpha_FDR,
          ". Showing top ", fallback_edges, " edges by |r| as fallback.")
  Rt <- R_top; diag(Rt) <- NA
  ord_edges <- order(abs(Rt), decreasing = TRUE, na.last = NA)
  show_idx <- ord_edges[seq_len(min(fallback_edges, length(ord_edges)))]
  R_sig[,] <- NA_real_
  R_sig[show_idx] <- R_top[show_idx]
}

# Drop rows/cols with no visible links and cluster
keep_idx <- rowSums(!is.na(R_sig)) > 0
R_sig <- R_sig[keep_idx, keep_idx, drop = FALSE]

R_for_order <- R_sig; R_for_order[is.na(R_for_order)] <- 0
hc <- hclust(as.dist(1 - R_for_order), method = "average")

R_plot <- R_sig
if (show_upper_only) R_plot[lower.tri(R_plot, diag = TRUE)] <- NA_real_

pal  <- make_diverging(257)
brks <- seq(-1, 1, length.out = length(pal))


#####################
### Figure S9

# SVG
svg("heatmap_phylum_temp_TOP50_sigma8_FDR.svg", width = 10, height = 10)
pheatmap(
  mat = R_plot,
  color = pal, breaks = brks, na_col = "grey95",
  cluster_rows = hc, cluster_cols = hc,
  border_color = NA, legend = TRUE,
  show_rownames = TRUE, show_colnames = TRUE,
  fontsize_row = 11, fontsize_col = 11, angle_col = 90,
  treeheight_row = 50, treeheight_col = 50, fontsize = 12
)
dev.off()

# High-res PNG
png("heatmap_phylum_temp_TOP50_sigma8_FDR.png", width = 5200, height = 5000, res = 550)
pheatmap(
  mat = R_plot,
  color = pal, breaks = brks, na_col = "grey95",
  cluster_rows = hc, cluster_cols = hc,
  border_color = NA, legend = TRUE,
  show_rownames = TRUE, show_colnames = TRUE,
  fontsize_row = 11, fontsize_col = 11, angle_col = 90,
  treeheight_row = 50, treeheight_col = 50, fontsize = 12
)
dev.off()






















































