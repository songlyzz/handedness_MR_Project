library(dplyr)
library(ggplot2)
library(patchwork)

setwd("./MR_Cog_Handedness_Bidirectional")

main_df <- read.csv("bidirectional_mr_cog_handedness.csv", stringsAsFactors = FALSE)

target_tasks <- c("Memory", "SymbolDigit", "Tower")

df_cog_to_handedness <- main_df %>%
  filter(
    method == "Inverse variance weighted",
    exposure_label %in% target_tasks,
    outcome_label == "Handedness"
  ) %>%
  mutate(
    across(c(pval, or, or_lci95, or_uci95), as.numeric),
    is_sig = pval < 0.05,
    pval_str = ifelse(pval < 0.001,
                      formatC(pval, format = "e", digits = 2),
                      sprintf("%.3f", pval)),
    annot_text = paste0(sprintf("%.3f (%.3f, %.3f)", or, or_lci95, or_uci95), "\n", pval_str),
    exposure_label = factor(exposure_label, levels = rev(c("Memory", "SymbolDigit", "Tower"))),
    # Labels corresponding to user prompts
    exposure_name = case_when(
        exposure_label == "Memory" ~ "Pairs Matching",
        exposure_label == "SymbolDigit" ~ "Symbol Digit",
        exposure_label == "Tower" ~ "Tower Rearranging",
        TRUE ~ as.character(exposure_label)
    ),
    exposure_name = factor(exposure_name, levels = rev(c("Pairs Matching", "Symbol Digit", "Tower Rearranging")))
  )

max_dev_lh <- max(abs(df_cog_to_handedness$or_lci95 - 1), abs(df_cog_to_handedness$or_uci95 - 1), na.rm = TRUE)
xlim_lh <- c(1 - max_dev_lh * 1.1, 1 + max_dev_lh * 1.1)
y_exp <- expansion(add = c(0.5, 1.2))

p1_left <- ggplot(df_cog_to_handedness, aes(y = exposure_name, alpha = is_sig)) +
  geom_errorbarh(aes(xmin = or_lci95, xmax = or_uci95), height = 0, linewidth = 0.7, color = "#4DA6D0") +
  geom_point(aes(x = or), size = 3, shape = 21, fill = "#4DA6D0", color = "black", stroke = 0.8) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3), guide = "none") +
  scale_x_continuous(limits = xlim_lh, n.breaks = 4) +
  scale_y_discrete(expand = y_exp) +
  labs(title = "Cognition to Left-handedness", x = "Odds Ratio", y = NULL) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(face = "bold", size = 10, color = "black"),
    plot.margin = margin(t = 20, r = 0, b = 10, l = 10)
  )

p1_right <- ggplot(df_cog_to_handedness, aes(y = exposure_name)) +
  geom_text(aes(x = 0.5, label = annot_text), hjust = 0.5, size = 3.5, color = "black", lineheight = 1.2) +
  annotate("text", x = 0.5, y = nlevels(df_cog_to_handedness$exposure_name) + 0.8,
           label = "OR (95% CI)\nP-value", fontface = "bold", size = 3.8, lineheight = 1.2) +
  scale_y_discrete(expand = y_exp) +
  theme_void() +
  theme(plot.margin = margin(t = 20, r = 10, b = 10, l = 0))

panel_cog_to_lh <- p1_left + p1_right + plot_layout(widths = c(2, 1))

df_handedness_to_cog <- main_df %>%
  filter(
    method == "Inverse variance weighted",
    exposure_label == "Handedness",
    outcome_label %in% target_tasks
  ) %>%
  mutate(
    across(c(pval, b, lo_ci, up_ci), as.numeric),
    is_sig = pval < 0.05,
    pval_str = ifelse(pval < 0.001,
                      formatC(pval, format = "e", digits = 2),
                      sprintf("%.3f", pval)),
    annot_text = paste0(sprintf("%.3f (%.3f, %.3f)", b, lo_ci, up_ci), "\n", pval_str),
    outcome_label = factor(outcome_label, levels = rev(c("Memory", "SymbolDigit", "Tower"))),
    outcome_name = case_when(
        outcome_label == "Memory" ~ "Pairs Matching",
        outcome_label == "SymbolDigit" ~ "Symbol Digit",
        outcome_label == "Tower" ~ "Tower Rearranging",
        TRUE ~ as.character(outcome_label)
    ),
    outcome_name = factor(outcome_name, levels = rev(c("Pairs Matching", "Symbol Digit", "Tower Rearranging")))

  )

max_dev_cog <- max(abs(df_handedness_to_cog$lo_ci), abs(df_handedness_to_cog$up_ci), na.rm = TRUE)
xlim_cog <- c(-max_dev_cog * 1.1, max_dev_cog * 1.1)

p2_left <- ggplot(df_handedness_to_cog, aes(y = outcome_name, alpha = is_sig)) +
  geom_errorbarh(aes(xmin = lo_ci, xmax = up_ci), height = 0, linewidth = 0.7, color = "#D05B43") +
  geom_point(aes(x = b), size = 3, shape = 21, fill = "#D05B43", color = "black", stroke = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3), guide = "none") +
  scale_x_continuous(limits = xlim_cog, n.breaks = 4) +
  scale_y_discrete(expand = y_exp) +
  labs(title = "Left-handedness to Cognition", x = "Beta", y = NULL) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(face = "bold", size = 10, color = "black"),
    plot.margin = margin(t = 20, r = 0, b = 10, l = 10)
  )

p2_right <- ggplot(df_handedness_to_cog, aes(y = outcome_name)) +
  geom_text(aes(x = 0.5, label = annot_text), hjust = 0.5, size = 3.5, color = "black", lineheight = 1.2) +
  annotate("text", x = 0.5, y = nlevels(df_handedness_to_cog$outcome_name) + 0.8,
           label = "Beta (95% CI)\nP-value", fontface = "bold", size = 3.8, lineheight = 1.2) +
  scale_y_discrete(expand = y_exp) +
  theme_void() +
  theme(plot.margin = margin(t = 20, r = 10, b = 10, l = 0))

panel_lh_to_cog <- p2_left + p2_right + plot_layout(widths = c(2, 1))

final_plot <- panel_cog_to_lh | panel_lh_to_cog
final_plot <- final_plot + plot_annotation(
  title = "Bidirectional Mendelian Randomization: Left-handedness and Cognition (IVW only)",
  theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 20)))
)

ggsave("MR_Cog_Handedness_Bidirectional_IVW.png", final_plot, width = 12, height = 4, dpi = 300, bg = "white")
ggsave("MR_Cog_Handedness_Bidirectional_IVW.pdf", final_plot, width = 12, height = 4, dpi = 300, bg = "white")