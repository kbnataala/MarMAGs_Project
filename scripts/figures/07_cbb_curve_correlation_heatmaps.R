################################################################################
# Script: 06_cbb_curve_correlation_heatmaps.R
#
# Purpose:
#   Documentation template for generating smoothed-curve co-occurrence
#   heatmaps among CBB-capable phyla across environmental gradients.
#
# Repository role:
#   This script is provided as a documentation template and may require
#   adaptation for local execution. It is not presented as a fully reproducible
#   pipeline component.
#
# Expected input:
#   - path/to/input/table
#
# Expected outputs:
#   - results/figures/CBB_temperature_curvecorr_heatmap_sigma8.png
#   - results/figures/CBB_depth_curvecorr_heatmap_sigma8.png
#   - results/figures/CBB_salinity_curvecorr_heatmap_sigma8.png
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(RColorBrewer)
  library(cowplot)
  library(gridExtra)
})

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------

base_dir <- "."
input_file <- file.path(base_dir, "results", "tables", "cbb_data.csv")
table_dir <- file.path(base_dir, "results", "tables")
figure_dir <- file.path(base_dir, "results", "figures")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) {
  stop("Missing input file: ", input_file)
}

# ------------------------------------------------------------------------------
# Parameters
# ------------------------------------------------------------------------------

taxon_col <- "GTDB_Tk_Phylum"
sigma_bins <- 8
n_perm <- 5000
alpha_lbl <- 0.05
set.seed(7)

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

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

check_monotone <- function(v, tol = 0) {
  dv <- diff(v)
  inc <- all(dv >= -tol)
  dec <- all(dv <= tol)

  if (inc && !dec) return(c(TRUE, "non-decreasing"))
  if (dec && !inc) return(c(TRUE, "non-increasing"))
  if (inc && dec)  return(c(TRUE, "constant"))
  c(FALSE, "non-monotonic")
}

gaussian_kernel <- function(sigma_bins) {
  half <- max(1L, round(3 * sigma_bins))
  x <- -half:half
  k <- exp(-(x^2) / (2 * sigma_bins^2))
  k / sum(k)
}

smooth_reflect <- function(y, k) {
  pad <- floor((length(k) - 1) / 2)
  ypad <- c(rev(y[2:(pad + 1)]), y, rev(y[(length(y) - pad):(length(y) - 1)]))
  yy <- stats::filter(ypad, k, sides = 2)
  as.numeric(yy[(pad + 1):(pad + length(y))])
}

perm_pvalue_corr <- function(x, y, L, n_perm) {
  shifts <- sample.int(L, n_perm, replace = TRUE) - 1L
  idx <- (rep(0:(L - 1), times = n_perm) - rep(shifts, each = L)) %% L + 1L
  yperm <- matrix(y[idx], nrow = L, ncol = n_perm)
  r_perm <- colSums(x * yperm) / L
  r_obs  <- sum(x * y) / L
  p <- (1 + sum(abs(r_perm) >= abs(r_obs))) / (1 + n_perm)
  c(r_obs = r_obs, p = p)
}

run_curve_corr_analysis <- function(data, gradient_col, prefix) {

  dat <- data %>%
    filter(.data[[gradient_col]] != "Unknown", !is.na(.data[[gradient_col]])) %>%
    mutate(!!taxon_col := trimws(as.character(.data[[taxon_col]])))

  dat$grad_num <- vapply(dat[[gradient_col]], parse_num, numeric(1))

  if (mean(is.na(dat$grad_num)) > 0.5) {
    levs <- unique(dat[[gradient_col]])
    dat$grad_num <- as.numeric(factor(dat[[gradient_col]], levels = levs))
  }

  pres_long <- dat %>%
    distinct(grad_num, !!rlang::sym(taxon_col)) %>%
    mutate(Presence = 1L) %>%
    complete(
      grad_num = sort(unique(grad_num)),
      !!rlang::sym(taxon_col),
      fill = list(Presence = 0L)
    ) %>%
    arrange(grad_num, !!rlang::sym(taxon_col))

  A_wide <- pres_long %>%
    pivot_wider(
      names_from = !!rlang::sym(taxon_col),
      values_from = Presence,
      values_fill = 0
    ) %>%
    arrange(grad_num)

  grad_vec <- A_wide$grad_num
  A <- as.matrix(A_wide[, -1, drop = FALSE])
  phyla <- colnames(A)
  L <- nrow(A)
  P <- ncol(A)

  # Monotonicity table
  mono_tbl <- lapply(seq_len(P), function(j) {
    v <- A[, j]
    md <- check_monotone(v)
    sp <- suppressWarnings(cor.test(rank(grad_vec), v, method = "spearman"))

    data.frame(
      Phylum = phyla[j],
      monotone = as.logical(md[1]),
      direction = md[2],
      spearman_rho = unname(sp$estimate),
      spearman_p = sp$p.value,
      n_levels = length(v),
      prevalence = sum(v),
      fraction_pres = mean(v),
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows()

  write.csv(
    mono_tbl,
    file.path(table_dir, paste0(prefix, "_monotonicity_by_phylum.csv")),
    row.names = FALSE
  )

  # Smoothed standardized abundance
  k <- gaussian_kernel(sigma_bins)
  Z <- apply(A, 2, smooth_reflect, k = k)
  Z <- scale(Z, center = TRUE, scale = TRUE)
  Z[is.na(Z)] <- 0

  R <- (t(Z) %*% Z) / L
  diag(R) <- 1
  R <- as.matrix(R)
  rownames(R) <- colnames(R) <- phyla

  # Permutation p-values
  pairs <- list()
  cnt <- 1L
  for (i in 1:(P - 1)) {
    xi <- Z[, i]
    for (j in (i + 1):P) {
      out <- perm_pvalue_corr(xi, Z[, j], L, n_perm)
      pairs[[cnt]] <- data.frame(
        Phylum1 = phyla[i],
        Phylum2 = phyla[j],
        curve_corr = out["r_obs"],
        p_value = out["p"],
        stringsAsFactors = FALSE
      )
      cnt <- cnt + 1L
    }
  }

  pairs_df <- bind_rows(pairs)

  # BH-FDR
  o <- order(pairs_df$p_value)
  r <- integer(length(o))
  r[o] <- seq_along(o)
  q <- pairs_df$p_value * nrow(pairs_df) / r
  q[o] <- rev(cummin(rev(q[o])))
  pairs_df$q_value <- q

  write.csv(
    pairs_df,
    file.path(table_dir, paste0(prefix, "_curve_pairs_sigma8_perm5000.csv")),
    row.names = FALSE
  )

  # Build heatmap data
  P_mat <- matrix(0, P, P, dimnames = list(phyla, phyla))
  apply(pairs_df, 1, function(row) {
    i <- row[["Phylum1"]]
    j <- row[["Phylum2"]]
    P_mat[i, j] <<- as.numeric(row[["p_value"]])
    P_mat[j, i] <<- as.numeric(row[["p_value"]])
  })

  cor_df <- as.data.frame(as.table(R))
  pval_df <- as.data.frame(as.table(P_mat))
  colnames(cor_df) <- c("Phylum1", "Phylum2", "Spearman_r")
  colnames(pval_df) <- c("Phylum1", "Phylum2", "p_value")

  merged_df <- inner_join(cor_df, pval_df, by = c("Phylum1", "Phylum2")) %>%
    filter(Phylum1 != Phylum2)

  write.csv(
    merged_df,
    file.path(table_dir, paste0("CBB_", prefix, "_curvecorr_pvalue.csv")),
    row.names = FALSE
  )

  div_pal <- brewer.pal(n = 11, name = "RdBu")

  heatmap_data <- merged_df %>%
    mutate(label = ifelse(p_value < alpha_lbl, round(Spearman_r, 2), NA))

  p <- ggplot(heatmap_data, aes(x = Phylum1, y = Phylum2, fill = Spearman_r)) +
    geom_tile(color = "grey90", linewidth = 0.3) +
    scale_fill_gradient2(
      low = div_pal[1],
      mid = div_pal[6],
      high = div_pal[11],
      midpoint = 0,
      limits = c(-1, 1),
      space = "Lab",
      na.value = "white",
      name = "Curve corr\n(σ=8)"
    ) +
    geom_text(
      aes(label = label),
      color = "black",
      size = 3.5,
      fontface = "bold",
      na.rm = TRUE
    ) +
    coord_fixed() +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "bold"),
      axis.text.y = element_text(face = "bold"),
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 10),
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, color = "grey50"),
      axis.ticks = element_blank(),
      plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
    )

  ggsave(
    file.path(figure_dir, paste0("CBB_", prefix, "_curvecorr_heatmap_sigma8.png")),
    p,
    width = 10,
    height = 8,
    dpi = 300
  )

  return(p)
}

# ------------------------------------------------------------------------------
# Run analyses
# ------------------------------------------------------------------------------

cbb_data <- read.csv(input_file, stringsAsFactors = FALSE)

temp_corr <- run_curve_corr_analysis(cbb_data, "sample_temp", "temperature")
depth_corr <- run_curve_corr_analysis(cbb_data, "sample_depth", "depth")
salinity_corr <- run_curve_corr_analysis(cbb_data, "sample_salinity", "salinity")

# ------------------------------------------------------------------------------
# Combined figure
# ------------------------------------------------------------------------------

common_theme <- theme_minimal(base_size = 20) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    legend.position = "right",
    legend.direction = "vertical",
    legend.title = element_text(face = "bold", size = 20),
    legend.text = element_text(size = 20),
    plot.margin = margin(6, 5, 6, 5)
  )

tag_theme <- theme(
  plot.tag.position = c(0, 1),
  plot.tag = element_text(face = "bold", size = 19)
)

x2 <- depth_corr + common_theme + labs(tag = "(b)") + tag_theme
y2 <- salinity_corr + common_theme + labs(tag = "(c)") + tag_theme
z2 <- temp_corr + common_theme + labs(tag = "(d)") + tag_theme

shared_leg <- get_legend(x2)

x2 <- x2 + theme(legend.position = "none")
y2 <- y2 + theme(legend.position = "none")
z2 <- z2 + theme(legend.position = "none")

panel_row <- arrangeGrob(x2, y2, z2, ncol = 3)

svg(file.path(figure_dir, "merged_cbb_correlations.svg"), width = 24, height = 10)
grid.arrange(
  panel_row,
  shared_leg,
  ncol = 2,
  widths = c(5, 1)
)
dev.off()
