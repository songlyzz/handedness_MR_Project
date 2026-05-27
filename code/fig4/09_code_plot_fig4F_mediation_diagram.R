library(ggplot2)
library(dplyr)

setwd("./MR_Mediation_Confounder")

# Path A: SCZ -> Memory
beta_a <- -0.0304
se_a <- 0.0041
pval_a <- 1.59e-13

# Path B: Memory -> LH
beta_b <- -0.1295
se_b <- 0.0440
pval_b <- 0.00328

# Total Effect: SCZ -> LH
beta_total <- 0.0145921183243691
se_total <- 0.00423633585149393
pval_total <- 0.000572086920598469

# Mediation Effect (A * B)
beta_m <- beta_a * beta_b
se_m <- sqrt( (beta_a^2 * se_b^2) + (beta_b^2 * se_a^2) )
z_m <- beta_m / se_m
pval_m <- 2 * pnorm(-abs(z_m))

# Direct Effect (Total - Mediation)
beta_direct <- beta_total - beta_m
se_direct <- sqrt(se_total^2 + se_m^2)
pval_direct <- 2 * pnorm(-abs(beta_direct / se_direct))

format_p <- function(p) {
  if (p < 0.001) return(formatC(p, format="e", digits=2))
  return(sprintf("%.3f", p))
}

p <- ggplot() +
  # Nodes
  annotate("label", x = 1.0, y = 1.0, label = "      SCZ      ", size = 5.5, fontface = "bold", fill = "#D7BDE2", label.padding = unit(0.6, "lines"), label.r = unit(0.3, "lines"), color="black") +
  annotate("label", x = 2.0, y = 2.0, label = "Pairs Matching\n(Mediator)", size = 5.5, fontface = "bold", fill = "#F5CBA7", label.padding = unit(0.6, "lines"), label.r = unit(0.3, "lines"), color="black") +
  annotate("label", x = 3.0, y = 1.0, label = "       LH      ", size = 5.5, fontface = "bold", fill = "#AED6F1", label.padding = unit(0.6, "lines"), label.r = unit(0.3, "lines"), color="black") +

  # Arrows
  annotate("segment", x = 1.3, xend = 1.7, y = 1.3, yend = 1.7, arrow = arrow(length=unit(0.35, "cm"), type="closed"), linewidth = 1, color="gray30") +
  annotate("segment", x = 2.3, xend = 2.7, y = 1.7, yend = 1.3, arrow = arrow(length=unit(0.35, "cm"), type="closed"), linewidth = 1, color="gray30") +
  annotate("segment", x = 1.45, xend = 2.55, y = 1.0, yend = 1.0, arrow = arrow(length=unit(0.35, "cm"), type="closed"), linewidth = 1, color="gray30") +

  # Path A and Path B Labels
  annotate("text", x = 1.35, y = 1.65, label = sprintf("\u03B21 = %.3f\nP = %s", beta_a, format_p(pval_a)), fontface="bold", size = 4.5, color="black") +
  annotate("text", x = 2.65, y = 1.65, label = sprintf("\u03B22 = %.3f\nP = %s", beta_b, format_p(pval_b)), fontface="bold", size = 4.5, color="black") +
  
  # Direct Beta annotation
  annotate("text", x = 2.0, y = 0.85, label = sprintf("Direct \u03B2 = %.5f, P = %s", beta_direct, format_p(pval_direct)), fontface="bold", size = 4.8, color="black") +

  # Mediation Beta annotation
  annotate("text", x = 2.0, y = 1.25, label = sprintf("\u03B21 \u00D7 \u03B22 = %.5f, P = %s", beta_m, format_p(pval_m)), fontface="bold", size=4.8, color="black") +
  
  scale_x_continuous(limits = c(0.4, 3.6)) +
  scale_y_continuous(limits = c(0.5, 2.3)) +
  theme_void() +
  labs(title="Mendelian Randomization Mediation Path") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18, margin = margin(t=10, b=10)))

ggsave("MR_Mediation_SCZ_Memory_Handedness.png", plot=p, width=10, height=6, dpi=300, bg="white")
ggsave("MR_Mediation_SCZ_Memory_Handedness.pdf", plot=p, width=10, height=6, dpi=300, bg="white")