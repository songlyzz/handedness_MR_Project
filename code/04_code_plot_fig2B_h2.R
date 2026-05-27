library(dplyr)
library(ggplot2)
library(patchwork)

# ----------------------------------------------------
# 1. Load/Construct Data
# ----------------------------------------------------
# Data derived from LDSC log & rg output, previously in plot_figure1AB.py
df <- data.frame(
  trait = c("LH", "SCZ", "ADHD", "PTSD"),
  h2    = c(0.0106, 0.1008, 0.0947, 0.0282),
  se    = c(0.0010, 0.0035, 0.0045, 0.0011)
)

df <- df %>%
  mutate(
    # Truncate lower CI at 0 if needed
    lci = ifelse(h2 - 1.96 * se < 0, 0, h2 - 1.96 * se),
    uci = h2 + 1.96 * se,
    # Create text annotation for the table panel
    annot_text = sprintf("%.4f (%.4f, %.4f)", h2, lci, uci)
  )

# Factor order: order determines bottom-to-top rendering in ggplot
trait_levels <- c("LH", "PTSD", "ADHD", "SCZ")
df$trait <- factor(df$trait, levels = trait_levels)

# ----------------------------------------------------
# 2. Styling
# ----------------------------------------------------
trait_colors <- c(
  "LH"    = "#4DA6D0",
  "SCZ"   = "#8B4DAB",
  "ADHD"  = "#D94040",
  "PTSD"  = "#34A86B"
)

y_exp <- expansion(add = c(0.5, 0.8))

# ----------------------------------------------------
# 3. Build Panels
# ----------------------------------------------------
# --- Forest Plot (Left Panel) ---
p_forest <- ggplot(df, aes(y = trait)) +
  geom_errorbarh(
    aes(xmin = lci, xmax = uci, color = trait),
    height = 0.25, linewidth = 1.2
  ) +
  geom_point(
    aes(x = h2, fill = trait), size = 3.5, shape = 21,
    color = "white", stroke = 1.0
  ) +
  geom_vline(
    xintercept = 0, linetype = "dashed", 
    color = "gray50", linewidth = 0.5
  ) +
  scale_color_manual(values = trait_colors, guide = "none") +
  scale_fill_manual(values = trait_colors, guide = "none") +
  scale_x_continuous(limits = c(-0.005, 0.12), n.breaks = 5) +
  scale_y_discrete(expand = y_exp) +
  labs(
    title = "B", 
    x = expression(h[SNP]^2~"(observed scale)"), 
    y = NULL
  ) +
  theme_classic() +
  theme(
    plot.title      = element_text(face = "bold", size = 16, hjust = 0),
    axis.line.y     = element_blank(),
    axis.ticks.y    = element_blank(),
    axis.text.y     = element_text(face = "bold", size = 11, color = "black"),
    axis.text.x     = element_text(size = 10),
    axis.title.x    = element_text(size = 12),
    panel.grid      = element_blank(),
    plot.margin     = margin(t = 10, r = 10, b = 10, l = 10)
  )

# ----------------------------------------------------
# 4. Save
# ----------------------------------------------------
ggsave("figure1B_h2_forest.pdf", p_forest, width = 5.0, height = 3.2)
ggsave("figure1B_h2_forest.png", p_forest, width = 5.0, height = 3.2, dpi = 300)
print("Saved: figure1B_h2_forest.pdf and figure1B_h2_forest.png")
