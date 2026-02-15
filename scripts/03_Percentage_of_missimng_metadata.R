##########################
# set working directory
setwd("D:/Thesis/work_package_3/Data_analysis/01_figures/10_missingness")

# =======================
# Missingness (flipped)
# =======================
suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

# ---- Inputs ----
in_csv     <- "MarMAGs_Additional_File_01.csv"  # your file
out_prefix <- "metadata_missingness_flipped"

# What counts as missing?
is_missing <- function(x) {
  y <- if (is.character(x)) str_trim(x) else x
  is.na(y) | (is.character(y) & y == "")
}

# ---- Load ----
raw <- suppressMessages(readr::read_csv(in_csv, show_col_types = FALSE))

# Best-effort sample labels from the first non-all-missing of the first 2 cols
sample_id <- {
  cands <- names(raw)[seq_len(min(2, ncol(raw)))]
  lab <- NULL
  for (nm in cands) {
    v <- raw[[nm]]
    if (!all(is_missing(v))) { lab <- as.character(v); break }
  }
  if (is.null(lab)) sprintf("Sample_%03d", seq_len(nrow(raw))) else lab
}
# Make labels unique to avoid "factor level duplicated" error
sample_id_unique <- make.unique(sample_id, sep = "_dup")

# Exclude first two columns from analysis
dat <- raw[, -(1:min(2, ncol(raw))), drop = FALSE]

# 0/1 missingness matrix (1 = missing)
miss_mat <- as.data.frame(lapply(dat, is_missing))
miss_mat[] <- lapply(miss_mat, as.integer)

n_samples <- nrow(miss_mat)
n_fields  <- ncol(miss_mat)

# =========================
# Panel A — field-wise bars (FLIPPED)
# =========================
field_summary <- tibble(
  Field   = names(miss_mat),
  Missing = colSums(miss_mat, na.rm = TRUE)
) |>
  mutate(Total = n_samples,
         Pct   = Missing / Total,
         Label = sprintf("%d (%.1f%%)", Missing, 100*Pct)) |>
  arrange(desc(Pct)) |>
  mutate(Field = factor(Field, levels = Field))   # keep sorted order

panelA <- ggplot(field_summary, aes(y = Field, x = Pct)) +
  geom_col(fill = "#2C7FB8", width = 0.72) +
  geom_text(aes(label = Label), hjust = -0.05, size = 4.2, fontface = "bold") +
  scale_x_continuous(labels = percent_format(accuracy = 1),
                     expand  = expansion(mult = c(0, 0.05))) +
  labs(x = "Missing values", y = "Metadata field", title = "") +
  theme_minimal(base_size = 20) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.title = element_text(face = "bold")
  ) +
  coord_cartesian(xlim = c(0, max(field_summary$Pct) + 0.10))

# autosize (height now depends on #fields)
wA <- 14
hA <- max(6, 0.24 * n_fields + 2)
ggsave(paste0(out_prefix, "_panelA_all.png"), panelA, width = wA, height = hA, dpi = 600)
ggsave(paste0(out_prefix, "_panelA_all.svg"), panelA, width = wA, height = hA)

# =========================
# Panel B — missingness matrix (FLIPPED)
#   y = Field, x = Sample
# =========================
# Columns ordered like Panel A; rows (samples) by total missingness
col_order <- levels(field_summary$Field)
row_missing <- rowSums(miss_mat, na.rm = TRUE)
row_order   <- order(row_missing, decreasing = TRUE)

miss_long <- miss_mat |>
  mutate(.sample = sample_id_unique) |>
  mutate(.sample = factor(.sample, levels = sample_id_unique[row_order])) |>
  pivot_longer(-.sample, names_to = "Field", values_to = "Missing") |>
  mutate(
    Field   = factor(Field, levels = col_order),
    Missing = factor(Missing, levels = c(0,1), labels = c("Present","Missing"))
  )

panelB <- ggplot(miss_long, aes(x = .sample, y = Field, fill = Missing)) +
  geom_tile(width = 0.98, height = 0.98) +
  scale_fill_manual(values = c("Present" = "#F1F1F1", "Missing" = "#2C7FB8")) +
  guides(fill = guide_legend(override.aes = list(color = NA))) +
  labs(x = "Sample", y = "Metadata field", fill = NULL, title = "Sample × field missingness") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid  = element_blank(),
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )

# autosize (width ~ #samples, height ~ #fields; cap to 45 in)
wB <- min(max(10, 0.10 * n_samples + 2), 45)
hB <- min(max(10, 0.26 * n_fields  + 2), 45)

ggsave(paste0(out_prefix, "_panelB_all.png"), panelB, width = wB, height = hB, dpi = 300)
ggsave(paste0(out_prefix, "_panelB_all.svg"), panelB, width = wB, height = hB)

