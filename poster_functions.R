plot_dpcrtest <- function(test_stats) {
  test_df <- data.frame(t(rbind(rep(c(1, 1:5*10, 2:10*50), 20), 
                                as.vector(sapply(c(0:5*10, 2:10*50, 6:10*100), rep, 15)), 
                                test_stats)))
  colnames(test_df) <- c("added_molecules", "base_number", "mean", "sd")
  
  test_df[["added_molecules"]] <- as.factor(test_df[["added_molecules"]])
  test_df[["base_number"]] <- as.factor(test_df[["base_number"]])
  test_df[["mean"]] <- cut(test_df[["mean"]], c(0, 0.001, 0.01, 0.05, 0.25, 0.5, 1), include.lowest = TRUE)
  ggplot(test_df, aes(x = base_number, y = added_molecules)) +
    scale_fill_manual(name = "mean p-value", 
                      values = c("#0072B2", "#56B4E9", "#009E73", 
                                 "#F0E442", "#E69F00", "#D55E00")) +
    geom_tile(aes(x = base_number, y = added_molecules, fill = mean)) +
    theme(legend.background = element_rect(fill="NA")) +
    geom_point(aes(x = base_number, y = added_molecules, size = sd), range = c(3, 9)) +
    scale_size_continuous(name = "standard deviation\nof p-value", 
                          breaks = c(0, 0.005, 0.05, 0.5)) +
    scale_x_discrete("Base number of molecules") +
    scale_y_discrete("Added number of molecules") +
    theme(plot.background=element_rect(fill = "transparent",colour = "transparent"),
          axis.text.x = element_text(size=16), axis.text.y = element_text(size=16),
          axis.title.x = element_text(size=19), axis.title.y = element_text(size=19))
}

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}