# ============================================================
# Redox / bioenergetics bubbleplot
#
# Raw data (xlsx or csv) -> per-sample normalisation -> statistics -> bubbleplots
#
#   x axis : assay          y axis : species (adults/juveniles, then larvae)
#   facets : time (rows) x treatment (columns)
#   colour : mean log2 fold change  OR  Hedges' g
#   size   : FDR-adjusted p-value class (black outline = p.adj <= 0.05)
#
# Usage
#   Rscript redox_bubbleplot.R <input.xlsx|input.csv> <output_folder>
#   or run interactively (e.g. in RStudio): you will be asked for both paths.
#
# See README.md for the analysis and the input format.
# ============================================================

# ---- 0) Settings ----

# Treatment codes in the input file -> labels used in tables and figures
TREATMENT_RENAME <- c("DB-S" = "Vessel-DNS", "MSFD-S" = "MSFD-NB")
CONTROL_LABEL    <- "Control"
TIME_POINTS_H    <- c(1, 24)

# Species codes in the input file -> labels on the y axis (top to bottom)
SPECIES_LABELS <- c(
  "Blanceolatum"  = "B. lanceolatum",
  "Bschlosseri"   = "B. schlosseri",
  "Clepadiformis" = "C. lepadiformis",
  "Crobusta"      = "C. robusta",
  "Cintestinalis" = "C. intestinalis"
)

# Assays: order on the x axis (three groups) and labels.
# "P/O" is not read from the input: it is computed as ATP synthesis / OCR.
ASSAY_LABELS <- c(
  "ATP synthesis" = "ATP synthesis",
  "OCR"           = "Oxygen consumption rate",
  "P/O"           = "OxPhos efficiency",
  "GR"            = "Glutathione reductase",
  "CAT"           = "Catalase",
  "NADPHOX"       = "NADPH oxidase",
  "MDA"           = "Lipid peroxidation",
  "Nitrotyrosine" = "Protein oxidation",
  "8-OHdG"        = "DNA damage"
)
ASSAY_GROUP_BOUNDS <- c(3.5, 6.5)   # dashed lines between the three assay groups

# Rows whose `notes` contain "larva" are analysed as a separate unit
# (own controls, own tests). Set to TRUE to pool larvae with the adults
# of the same species (proof-of-concept analysis).
POOL_LIFE_STAGES <- FALSE

# Figure format (Science, two columns)
FIG_WIDTH_CM  <- 12.1
FIG_HEIGHT_CM <- 8
FIG_DPI       <- 600
FONT_PT_AXIS   <- 6
FONT_PT_STRIP  <- 7
FONT_PT_LEGEND <- 6
SYMBOL_MIN_PT  <- 6     # smallest bubble; ggplot sizes are in mm (1 pt = 1/2.845 mm)
SYMBOL_STEP_MM <- 0.8

# ---- 1) Packages ----
packages <- c("tidyverse", "readxl", "emmeans", "broom", "scales")
for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(emmeans)
  library(broom)
  library(scales)
})

# ---- 2) Input and output paths ----
choose_paths <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 2) return(list(input = args[1], output = args[2]))

  if (!interactive()) {
    stop("Usage: Rscript redox_bubbleplot.R <input.xlsx|input.csv> <output_folder>")
  }
  has_rstudio <- requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()

  message("Select the input data file (.xlsx or .csv)")
  input <- if (has_rstudio) {
    rstudioapi::selectFile(caption = "Input data (.xlsx or .csv)",
                           filter = "Data files (*.xlsx *.csv)")
  } else {
    file.choose()
  }
  if (is.null(input) || !nzchar(input)) stop("No input file selected.")

  default_out <- file.path(dirname(input), "bubbleplot_output")
  message("Select the output folder")
  output <- if (has_rstudio) {
    rstudioapi::selectDirectory(caption = "Output folder", path = dirname(input))
  } else {
    ans <- readline(sprintf("Output folder [%s]: ", default_out))
    if (nzchar(ans)) ans else default_out
  }
  if (is.null(output) || !nzchar(output)) stop("No output folder selected.")
  list(input = input, output = output)
}

paths      <- choose_paths()
out_dir    <- paths$output
fig_dir    <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
message("Input : ", paths$input)
message("Output: ", out_dir)

# ---- 3) Read and validate input ----
read_input <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx", "xls")) {
    read_excel(path, sheet = 1)
  } else if (ext == "csv") {
    read_csv(path, col_types = cols(.default = col_character()), show_col_types = FALSE)
  } else {
    stop("Unsupported file type: ", ext, " (use .xlsx or .csv)")
  }
}

data_raw <- read_input(paths$input)

required_cols <- c("sample_id", "species", "treatment", "time_h", "assay",
                   "value", "unit", "replicate_id", "batch", "notes")
missing_cols <- setdiff(required_cols, names(data_raw))
if (length(missing_cols) > 0) stop("Missing columns: ", paste(missing_cols, collapse = ", "))

with_rows <- data_raw %>% mutate(excel_row = row_number() + 1)

bad_time <- with_rows %>% filter(!(as.character(time_h) %in% as.character(TIME_POINTS_H)))
if (nrow(bad_time) > 0) {
  print(as.data.frame(bad_time %>% select(excel_row, sample_id, species, treatment, time_h, assay)))
  stop("Unexpected time_h values in the rows above: fix the input file and run again.")
}

known_trt <- c(CONTROL_LABEL, names(TREATMENT_RENAME))
bad_trt <- with_rows %>% filter(!(treatment %in% known_trt))
if (nrow(bad_trt) > 0) {
  print(as.data.frame(bad_trt %>% distinct(treatment)))
  stop("Unknown treatment codes (allowed: ", paste(known_trt, collapse = ", "), ").")
}

if (any(data_raw$assay == "OxPhos efficiency", na.rm = TRUE)) {
  message("Ignoring 'OxPhos efficiency' rows: P/O is computed from ATP synthesis / OCR.")
}

LARVAE <- "larvae"
ADULTS <- "adults/juveniles"
TRT_LEVELS <- unname(TREATMENT_RENAME)

is_larva <- str_detect(coalesce(data_raw$notes, ""), regex("larva", ignore_case = TRUE))
species_with_larvae <- unique(data_raw$species[is_larva])

data <- data_raw %>%
  filter(assay != "OxPhos efficiency") %>%
  mutate(
    value      = as.numeric(value),
    treatment  = recode(as.character(treatment), !!!TREATMENT_RENAME),
    treatment  = factor(treatment, levels = c(CONTROL_LABEL, TRT_LEVELS)),
    time_h     = factor(as.character(time_h), levels = as.character(TIME_POINTS_H)),
    species    = as.character(species),
    assay      = as.character(assay),
    batch      = as.factor(batch),
    life_stage = if_else(!POOL_LIFE_STAGES &
                           str_detect(coalesce(notes, ""), regex("larva", ignore_case = TRUE)),
                         LARVAE, ADULTS)
  ) %>%
  filter(!is.na(value)) %>%
  # larvae become a separate analysis unit ("<species>_larva")
  mutate(species = if_else(life_stage == LARVAE & !str_detect(species, "_larva$"),
                           paste0(species, "_larva"), species))

stage_check <- data %>% distinct(species, life_stage) %>% count(species) %>% filter(n > 1)
if (nrow(stage_check) > 0) {
  stop("Analysis units with mixed life stages: ", paste(stage_check$species, collapse = ", "))
}
message("Analysis units:")
print(data %>% distinct(species, life_stage) %>% arrange(life_stage, species))

# ---- 4) Average technical replicates -> one value per biological sample ----
data_avg <- data %>%
  group_by(species, life_stage, assay, time_h, treatment, batch, sample_id) %>%
  summarise(value = mean(value, na.rm = TRUE), n_technical = n(), .groups = "drop")

# ---- 5) P/O = ATP synthesis / OCR (per biological sample) ----
po <- data_avg %>%
  filter(assay %in% c("ATP synthesis", "OCR")) %>%
  select(species, life_stage, time_h, treatment, batch, sample_id, assay, value) %>%
  pivot_wider(names_from = assay, values_from = value)

po_incomplete <- po %>% filter(is.na(`ATP synthesis`) | is.na(OCR))
if (nrow(po_incomplete) > 0) {
  message("P/O not computed (ATP or OCR missing) for ", nrow(po_incomplete), " samples:")
  print(as.data.frame(po_incomplete %>% select(species, time_h, treatment, sample_id)))
}

po <- po %>%
  filter(!is.na(`ATP synthesis`), !is.na(OCR), OCR > 0) %>%
  transmute(species, life_stage, time_h, treatment, batch, sample_id,
            assay = "P/O", value = `ATP synthesis` / OCR, n_technical = 1L)

data_avg <- bind_rows(data_avg, po)
write.csv(po, file.path(out_dir, "01_PO_calculated.csv"), row.names = FALSE)

# ---- 6) Normalisation to controls (batch-matched, pooled as fallback) ----
controls_batch <- data_avg %>%
  filter(treatment == CONTROL_LABEL) %>%
  group_by(species, assay, time_h, batch) %>%
  summarise(ctrl_mean_batch = mean(value, na.rm = TRUE), .groups = "drop")

controls_pooled <- data_avg %>%
  filter(treatment == CONTROL_LABEL) %>%
  group_by(species, assay, time_h) %>%
  summarise(ctrl_mean_pooled = mean(value, na.rm = TRUE), .groups = "drop")

data_norm <- data_avg %>%
  left_join(controls_batch, by = c("species", "assay", "time_h", "batch")) %>%
  left_join(controls_pooled, by = c("species", "assay", "time_h")) %>%
  mutate(
    ctrl_ref = if_else(!is.na(ctrl_mean_batch), ctrl_mean_batch, ctrl_mean_pooled),
    fc       = value / ctrl_ref,
    log2_fc  = log2(fc)
  ) %>%
  filter(!is.na(ctrl_ref), is.finite(fc), is.finite(log2_fc))

write.csv(data_norm, file.path(out_dir, "02_data_normalized_log2FC.csv"), row.names = FALSE)

# ---- 7) Statistics: each treatment vs control (unadjusted) ----
stats_raw <- data_norm %>%
  group_by(species, life_stage, assay, time_h) %>%
  group_modify(~ {
    df <- .x
    if (n_distinct(df$treatment) < 2 || !(CONTROL_LABEL %in% df$treatment)) return(tibble())
    model <- lm(log2_fc ~ treatment, data = df)
    emmeans(model, ~ treatment) %>%
      contrast(method = "trt.vs.ctrl", ref = CONTROL_LABEL, adjust = "none") %>%
      broom::tidy()
  }) %>%
  ungroup() %>%
  mutate(treatment = map_chr(contrast, function(ct) {
    hit <- TRT_LEVELS[str_detect(ct, fixed(TRT_LEVELS))]
    if (length(hit) == 1) hit else NA_character_
  })) %>%
  filter(!is.na(treatment)) %>%
  select(-term, -null.value)

# ---- 8) Effect size: Cohen's d and Hedges' g (treatment vs its controls) ----
effect_sizes <- data_norm %>%
  group_by(species, life_stage, assay, time_h) %>%
  group_modify(~ {
    df <- .x
    if (!(CONTROL_LABEL %in% df$treatment) || n_distinct(df$treatment) < 2) return(tibble())
    ctrl_vals <- df %>% filter(treatment == CONTROL_LABEL) %>% pull(log2_fc)
    trts <- df %>% filter(treatment != CONTROL_LABEL) %>% pull(treatment) %>% unique() %>% as.character()
    map_dfr(trts, function(trt) {
      trt_vals <- df %>% filter(treatment == trt) %>% pull(log2_fc)
      n1 <- length(trt_vals); n2 <- length(ctrl_vals)
      if (n1 < 2 || n2 < 2) return(tibble())
      df_pooled <- n1 + n2 - 2
      pooled_sd <- sqrt(((n1 - 1) * sd(trt_vals)^2 + (n2 - 1) * sd(ctrl_vals)^2) / df_pooled)
      if (!is.finite(pooled_sd) || pooled_sd == 0) return(tibble())
      d <- (mean(trt_vals) - mean(ctrl_vals)) / pooled_sd
      J <- 1 - 3 / (4 * df_pooled - 1)   # Hedges' small-sample correction
      tibble(treatment = trt, cohens_d = d, hedges_g = d * J, n_treatment = n1, n_control = n2)
    })
  }) %>%
  ungroup()

# ---- 9) Labels ----
# y-axis labels as plotmath: species name in italics, "larvae" upright
species_plotmath <- map_chr(names(SPECIES_LABELS), function(sp) {
  lab <- sprintf("italic('%s')", SPECIES_LABELS[[sp]])
  if (POOL_LIFE_STAGES && sp %in% species_with_larvae) lab <- paste0(lab, "~'(adults + larvae)'")
  lab
})
names(species_plotmath) <- names(SPECIES_LABELS)
larva_units <- paste0(names(SPECIES_LABELS), "_larva")
species_plotmath <- c(species_plotmath,
                      setNames(sprintf("italic('%s')~larvae", SPECIES_LABELS), larva_units))
species_order_top_down <- c(names(SPECIES_LABELS), larva_units)

unknown_species <- setdiff(unique(data_norm$species), species_order_top_down)
if (length(unknown_species) > 0) {
  stop("Species not listed in SPECIES_LABELS: ", paste(unknown_species, collapse = ", "))
}

assay_order <- names(ASSAY_LABELS)

sig_levels <- c("p > 0.05", "p ≤ 0.05", "p ≤ 0.01", "p ≤ 0.001", "p ≤ 0.0001")
make_sig_bin <- function(p) {
  factor(
    case_when(
      p > 0.05    ~ sig_levels[1],
      p <= 0.0001 ~ sig_levels[5],
      p <= 0.001  ~ sig_levels[4],
      p <= 0.01   ~ sig_levels[3],
      p <= 0.05   ~ sig_levels[2],
      TRUE ~ NA_character_
    ),
    levels = sig_levels
  )
}

symbol_min_mm   <- SYMBOL_MIN_PT / 2.845276
sig_bin_sizes   <- setNames(symbol_min_mm + (0:4) * SYMBOL_STEP_MM, sig_levels)
sig_bin_borders <- setNames(c("grey60", rep("black", 4)), sig_levels)

# ---- 10) Plotting ----
build_bubble_plot <- function(df, color_var, limits, legend_name) {
  sp_present <- species_order_top_down[species_order_top_down %in% unique(df$species)]

  p <- ggplot(df, aes(x = assay, y = species)) +
    geom_vline(xintercept = ASSAY_GROUP_BOUNDS, linetype = "dashed",
               color = "grey55", linewidth = 0.25)

  # separator between adults/juveniles and larvae
  n_larvae_rows <- sum(sp_present %in% (df %>% filter(life_stage == LARVAE) %>% pull(species)))
  if (n_larvae_rows > 0 && n_larvae_rows < length(sp_present)) {
    p <- p + geom_hline(yintercept = n_larvae_rows + 0.5, color = "grey40", linewidth = 0.3)
  }

  p +
    geom_point(aes(size = sig_bin_adj, fill = .data[[color_var]], color = sig_bin_adj),
               shape = 21, stroke = 0.3) +
    scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426",
                         midpoint = 0, limits = limits, oob = scales::squish,
                         name = legend_name, guide = guide_colorbar(order = 1)) +
    scale_color_manual(values = sig_bin_borders, drop = TRUE, name = "FDR p-adj",
                       guide = guide_legend(order = 2)) +
    scale_size_manual(values = sig_bin_sizes, drop = TRUE, name = "FDR p-adj",
                      guide = guide_legend(order = 2)) +
    scale_y_discrete(limits = rev(sp_present),
                     labels = function(x) parse(text = species_plotmath[x])) +
    scale_x_discrete(labels = ASSAY_LABELS) +
    facet_grid(time_label ~ treatment) +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = FONT_PT_AXIS) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = FONT_PT_AXIS),
      axis.text.y = element_text(color = "black", size = FONT_PT_AXIS),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = FONT_PT_STRIP),
      legend.text = element_text(size = FONT_PT_LEGEND),
      legend.title = element_text(size = FONT_PT_LEGEND + 1),
      legend.key.size = unit(3, "mm"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey93", linewidth = 0.2),
      panel.spacing.y = unit(0.5, "lines"),
      panel.spacing.x = unit(2/3, "lines"),
      plot.margin = margin(t = 2, r = 2, b = 2, l = 8, unit = "mm")
    )
}

# PDF device that can write "≤": cairo if it works, Quartz on macOS
# (cairo needs XQuartz there), base pdf otherwise.
pick_pdf_device <- function() {
  cairo_ok <- tryCatch({
    f <- tempfile(fileext = ".pdf")
    grDevices::cairo_pdf(f, width = 1, height = 1)
    grid::grid.text("≤")
    grDevices::dev.off()
    file.exists(f) && file.size(f) > 0
  }, warning = function(w) FALSE, error = function(e) FALSE)
  while (grDevices::dev.cur() > 1 && names(grDevices::dev.cur()) %in% c("cairo_pdf", "pdf")) grDevices::dev.off()

  if (cairo_ok) return(grDevices::cairo_pdf)
  if (Sys.info()[["sysname"]] == "Darwin") {
    return(function(filename, width, height, ...) {
      grDevices::quartz(type = "pdf", file = filename, width = width, height = height)
    })
  }
  message("Neither cairo nor Quartz available: using the base pdf device ('≤' may not render).")
  "pdf"
}
pdf_device <- pick_pdf_device()

save_bubble <- function(plot, stub) {
  w <- FIG_WIDTH_CM / 2.54; h <- FIG_HEIGHT_CM / 2.54
  ggsave(file.path(fig_dir, paste0(stub, ".pdf")), plot, width = w, height = h, device = pdf_device)
  ggsave(file.path(fig_dir, paste0(stub, ".png")), plot, width = w, height = h, dpi = FIG_DPI)
}

# ---- 11) One analysis set = FDR + table + 3 figures ----
# FDR (Benjamini-Hochberg) is computed over all tests of the set, and only those.
run_set <- function(stages, tag) {
  st <- stats_raw %>%
    filter(life_stage %in% stages) %>%
    mutate(p.adj = p.adjust(p.value, method = "BH"))
  es <- effect_sizes %>% filter(life_stage %in% stages)

  tab <- st %>%
    left_join(es, by = c("species", "life_stage", "assay", "time_h", "treatment")) %>%
    arrange(life_stage, species, assay, time_h, treatment)
  write.csv(tab, file.path(out_dir, paste0("03_stats_FDR_effectsizes_", tag, ".csv")), row.names = FALSE)

  bubble_df <- data_norm %>%
    filter(life_stage %in% stages, treatment != CONTROL_LABEL) %>%
    group_by(species, life_stage, assay, time_h, treatment) %>%
    summarise(mean_log2fc = mean(log2_fc, na.rm = TRUE), n_bio = n(), .groups = "drop") %>%
    mutate(treatment = as.character(treatment)) %>%
    left_join(tab %>% mutate(treatment = as.character(treatment)) %>%
                select(species, assay, time_h, treatment, p.value, p.adj, hedges_g),
              by = c("species", "assay", "time_h", "treatment")) %>%
    mutate(
      time_label  = factor(paste0(as.character(time_h), " h"), levels = paste0(TIME_POINTS_H, " h")),
      treatment   = factor(treatment, levels = TRT_LEVELS),
      assay       = factor(assay, levels = assay_order),
      sig_bin_adj = make_sig_bin(p.adj)
    ) %>%
    filter(!is.na(assay), !is.na(sig_bin_adj))

  max_abs_log2fc <- max(abs(bubble_df$mean_log2fc), na.rm = TRUE)
  g_cap   <- unname(quantile(abs(bubble_df$hedges_g), 0.9, na.rm = TRUE))   # 90th percentile
  g_nocap <- max(abs(bubble_df$hedges_g), na.rm = TRUE)

  save_bubble(build_bubble_plot(bubble_df, "mean_log2fc", c(-max_abs_log2fc, max_abs_log2fc),
                                expression(log[2]~"FC")),
              paste0("bubbleplot_log2FC_padj_", tag))
  save_bubble(build_bubble_plot(bubble_df, "hedges_g", c(-g_cap, g_cap), "Hedges' g"),
              paste0("bubbleplot_hedgesG_padj_cap_", tag))
  save_bubble(build_bubble_plot(bubble_df, "hedges_g", c(-g_nocap, g_nocap), "Hedges' g"),
              paste0("bubbleplot_hedgesG_padj_nocap_", tag))

  message(sprintf("[%s] tests: %d | max|log2FC| = %.2f | g cap (90th pct) = %.2f | max|g| = %.2f",
                  tag, nrow(st), max_abs_log2fc, g_cap, g_nocap))
  invisible(tab)
}

run_set(ADULTS, if (POOL_LIFE_STAGES) "pooledStages" else "noLarvae")
if (any(data_norm$life_stage == LARVAE)) run_set(c(ADULTS, LARVAE), "withLarvae")

message("Done. Tables in ", out_dir, ", figures in ", fig_dir)
