library(ggplot2)
library(dplyr)

setwd("./MR_Mediation_Confounder")

# Path 1: Memory -> LH
beta_1 <- -0.1295
se_1 <- 0.0440
pval_1 <- 0.00328

# Path 2: Memory -> PTSD
beta_2 <- -0.1146
se_2 <- 0.0442
pval_2 <- 0.00952

format_p <- function(p) {
  if (p < 0.001) return(formatC(p, format="e", digits=2))
  return(sprintf("%.3f", p))
}

p <- ggplot() +
  annotate("label", x = 1.0, y = 1.0, label = "       LH      ", size = 5.5, fontface = "bold", fill = "#AED6F1", label.padding = unit(0.6, "lines"), label.r = unit(0.3, "lines"), color="black") +
  annotate("label", x = 2.0, y = 2.0, label = "Pairs Matching\n(Confounder)", size = 5.5, fontface = "bold", fill = "#F5CBA7", label.padding = unit(0.6, "lines"), label.r = unit(0.3, "lines"), color="black") +
  annotate("label", x = 3.0, y = 1.0, label = "     PTSD      ", size = 5.5, fontface = "bold", fill = "#D7BDE2", label.padding = unit(0.6, "lines"), label.r = unit(0.3, "lines"), color="black") +

  # Arrows: Both downward, plus horizontal
  annotate("segment", x = 1.7, xend = 1.3, y = 1.7, yend = 1.3, arrow = arrow(length=unit(0.35, "cm"), type="closed"), linewidth = 1, color="gray30") +
  annotate("segment", x = 2.3, xend = 2.7, y = 1.7, yend = 1.3, arrow = arrow(length=unit(0.35, "cm"), type="closed"), linewidth = 1, color="gray30") +
  annotate("segment", x = 1.45, xend = 2.55, y = 1.0, yend = 1.0, arrow = arrow(length=unit(0.35, "cm"), type="closed"), linewidth = 1, color="gray30", linetype="dashed") +

  # Labels for the Confounder Paths Only
  annotate("text", x = 1.35, y = 1.65, label = sprintf("\u03B21 = %.3f\nP = %s", beta_1, format_p(pval_1)), fontface="bold", size = 4.5, color="black") +
  annotate("text", x = 2.65, y = 1.65, label = sprintf("\u03B22 = %.3f\nP = %s", beta_2, format_p(pval_2)), fontface="bold", size = 4.5, color="black") +
  
  scale_x_continuous(limits = c(0.4, 3.6)) +
  scale_y_continuous(limits = c(0.5, 2.3)) +
  theme_void() +
  labs(title="Mendelian Randomization Confounder Path") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18, margin = margin(t=10, b=10)))

ggsave("MR_Confounder_Memory_LH_PTSD.png", plot=p, width=10, height=6, dpi=300, bg="white")
ggsave("MR_Confounder_Memory_LH_PTSD.pdf", plot=p, width=10, height=6, dpi=300, bg="white")