library(dplyr)
library(ggplot2)
library(patchwork)

# ----------------------------------------------------
# 1. Load Data
# ----------------------------------------------------
main_df <- read.csv("MR_diseases_handedness_main.csv", stringsAsFactors = FALSE)

# ----------------------------------------------------
# 2. Process Data
# ----------------------------------------------------
df_plot <- main_df %>%
  filter(
    direction    == "disease_to_handedness",
    method_label %in% c("IVW", "MR-Egger", "Weighted_Median", "Weighted_Mode")
  ) %>%
  select(id.exposure, method_label, pval, or, or_lci95, or_uci95) %>%
  mutate(
    across(c(pval, or, or_lci95, or_uci95), as.numeric),
    # Add a flag for significance
    is_sig      = pval < 0.05,
    # If pval < 0.001 use scientific, else 3 decimals fixed
    pval_str    = ifelse(pval < 0.001, 
                         formatC(pval, format = "e", digits = 2), 
                         sprintf("%.3f", pval)),
    # Stack OR and P-value vertically without text labels
    annot_text  = paste0(sprintf("%.3f (%.3f, %.3f)", or, or_lci95, or_uci95), "\n", pval_str)
  )

# Factor order: top row = IVW, bottom row = Weighted_Mode
method_order <- rev(c("IVW", "MR-Egger", "Weighted_Median", "Weighted_Mode"))
df_plot$method_label <- factor(df_plot$method_label, levels = method_order)

disease_colors <- c(
  "ADHD"          = "#D92B03",
  "PTSD"          = "#038766",
  "Schizophrenia" = "#A77D9A"
)

# Shared y-scale expansion: top space for column headers
y_exp <- expansion(add = c(0.5, 1.2))

# ----------------------------------------------------
# 3. Build each panel as a standalone ggplot object
# ----------------------------------------------------

# --- Forest plot panel ---
make_forest <- function(dis, title_text, col) {
  df <- df_plot %>% filter(id.exposure == dis)
  max_dev <- max(abs(df$or_lci95 - 1), abs(df$or_uci95 - 1), na.rm = TRUE)
  max_dev <- max(max_dev * 1.2, 0.06)
  xlim    <- c(1 - max_dev, 1 + max_dev)

  ggplot(df, aes(y = method_label, alpha = is_sig)) +
    geom_errorbarh(
      aes(xmin = or_lci95, xmax = or_uci95),
      height = 0, linewidth = 0.7, color = col
    ) +
    geom_point(
      aes(x = or), size = 3, shape = 21,
      fill = col, color = "black", stroke = 0.8
    ) +
    geom_vline(
      xintercept = 1, linetype = "dashed",
      color = "gray50", linewidth = 0.5
    ) +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3), guide = "none") +
    scale_x_continuous(limits = xlim, n.breaks = 4) +
    scale_y_discrete(expand = y_exp) +
    labs(title = title_text, x = "Odds Ratio", y = NULL) +
    theme_classic() +
    theme(
      plot.title      = element_text(face = "bold", size = 12, hjust = 0.5),
      axis.line.y     = element_blank(),
      axis.ticks.y    = element_blank(),
      axis.text.y     = element_text(face = "bold", size = 10, color = "black"),
      axis.text.x     = element_text(size = 9),
      axis.title.x    = element_text(size = 10),
      panel.grid      = element_blank(),
      # right margin = 0 so table panel sits flush
      plot.margin     = margin(t = 8, r = 0, b = 8, l = 10)
    )
}

# --- Table annotation panel ---
# Headers placed once at the top via annotate(); rows contain bare values only.
make_table <- function(dis) {
  df    <- df_plot %>% filter(id.exposure == dis)
  n_lev <- nlevels(df$method_label)   # 4 levels, mapped to y = 1..4

  # y-position for column headers: just inside the top expansion zone
  header_y <- n_lev + 0.85

  ggplot(df, aes(y = method_label)) +
    # Stacked bare OR and P-values
    geom_text(
      aes(x = 0.5, label = annot_text),
      hjust = 0.5, size = 3.5, color = "black", lineheight = 1.2
    ) +
    # Single column header for the stacked values
    annotate(
      "text", x = 0.5, y = header_y,
      label = "OR (95% CI)\nP-value", fontface = "bold",
      hjust = 0.5, size = 3.8, lineheight = 1.1
    ) +
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_discrete(expand = y_exp) +
    # Invisible title & x-axis label to match the height of the forest panel
    labs(title = " ", x = " ") +
    theme_classic() +
    theme(
      plot.title    = element_text(size = 12),          # invisible but reserves space
      axis.line     = element_blank(),
      axis.ticks    = element_blank(),
      axis.text     = element_blank(),
      axis.title    = element_text(size = 10, color = "white"),  # hidden but takes space
      panel.grid    = element_blank(),
      # left margin = 0 so it sits flush against the forest panel
      plot.margin   = margin(t = 8, r = 10, b = 8, l = 0)
    )
}

# ----------------------------------------------------
# 4. Build all 6 panels and combine in ONE flat layout
# ----------------------------------------------------
diseases <- c("ADHD",             "PTSD",             "Schizophrenia")
titles   <- c("ADHD to Left-handedness",
              "PTSD to Left-handedness",
              "SCZ to Left-handedness")

forests <- mapply(make_forest, diseases, titles, disease_colors[diseases],
                  SIMPLIFY = FALSE)
tables  <- lapply(diseases, make_table)

# One plot_layout call - no nesting, avoids all overlap
final_plot <-
  forests[[1]] + tables[[1]] +
  forests[[2]] + tables[[2]] +
  forests[[3]] + tables[[3]] +
  plot_layout(nrow = 1, widths = c(2, 1.2, 2, 1.2, 2, 1.2)) +
  plot_annotation(title = "Mendelian Randomization: Disease to Left-handedness")

# Save
ggsave("MR_disease_to_handedness_methods_forest.pdf",
       final_plot, width = 18, height = 4.5)
ggsave("MR_disease_to_handedness_methods_forest.png",
       final_plot, width = 18, height = 4.5, dpi = 200)



