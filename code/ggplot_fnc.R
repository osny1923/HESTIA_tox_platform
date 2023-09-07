library(ggplot2)
plot_group_uncertainty <- function(grp, df, Tag) {
  plot_state <- {{grp}} == "Pesticides"
  max_row <- which.max(df$GStDev)
  sec_max_row <- which(df$GStDev == df$GStDev[order(df$GStDev, decreasing = TRUE)][2])
  
  plot <- df[df$Group == paste({{grp}}),] %>% 
    arrange(GStDev) %>% 
    ggplot(aes(x = n_recs, y = GStDev, color = status)) +
    geom_point() +
    scale_x_log10() +
    xlim(1, 200) +
    theme_bw() +
    guides(color = guide_legend(title = "Output status"))+
    theme(
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 10),
      legend.title.align = 0.5,
      plot.title = element_text(size = 10, face = "bold"),
      axis.text = element_text(color = "black", size = 10),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.line = element_line(colour = "black")
    ) +
    ggtitle(paste(Tag, grp, sep = " "))
  
  if (!plot_state) {
    plot <- plot + guides(color = "none")
  }
  return(plot)
}
