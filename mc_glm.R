library(dpcR)
library(multcomp)
library(pbapply)

test1 <- pblapply(1L:4000, function(dummy_variable) {
  sapply(c(0:5*10, 2:10*50, 6:10*100), function(base_number) {
    sapply(c(1, 1:5*10, 2:10*50), function(added_molecules) {
      dat <- sim_adpcr(m = base_number, n = 1000, times = 10000, pos_sums = FALSE, n_panels = 1) 
      dat2 <- dat
      ids <- sample(1L:1000, added_molecules)
      dat2[ids, ] <- dat2[ids, ] + 1
      #sample added below to reshuffle positive wells
      counts_data <- data.frame(experiment = c(rep("baseline", 1000), rep("modyfied", 1000)), 
                                counts = as.numeric(c(as.vector(dat) > 0, 
                                                      sample(as.vector(dat2)) > 0)))
      fit <- glm(counts ~ experiment + 0, data = counts_data, family = quasipoisson())
      summary(glht(fit, linfct = mcp(experiment = "Tukey")))[["test"]][["pvalues"]]
    })
  })
})

save(test1, file = "dpcrposter_data.RData")


test2 <- pblapply(1L:4000, function(dummy_variable) {
  sapply(c(0:5*10, 2:10*50, 6:10*100), function(base_number) {
    sapply(c(1, 1:5*10, 2:10*50), function(added_molecules) {
      dat <- sim_adpcr(m = base_number, n = 2000, times = 10000, pos_sums = FALSE, n_panels = 1) 
      dat2 <- dat
      ids <- sample(1L:2000, added_molecules)
      dat2[ids, ] <- dat2[ids, ] + 1
      #sample added below to reshuffle positive wells
      counts_data <- data.frame(experiment = c(rep("baseline", 1000), rep("modyfied", 1000)), 
                                counts = as.numeric(c(as.vector(dat) > 0, 
                                                      sample(as.vector(dat2)) > 0)))
      fit <- glm(counts ~ experiment + 0, data = counts_data, family = quasipoisson)
      summary(glht(fit, linfct = mcp(experiment = "Tukey")))[["test"]][["pvalues"]]
    })
  })
})

save(test1, test2, file = "dpcrposter_data.RData")












plot_dpcrtest <- function(test_data) {
  test1_stats <- sapply(1L:length(test_data[[1]]), function(position) {
    tmp <- sapply(test_data, function(test) test[position])
    c(mean(tmp), sd(tmp))
  })
  
  test_df <- data.frame(t(rbind(rep(c(1, 1:5*10, 2:10*50), 20), 
                                as.vector(sapply(c(0:5*10, 2:10*50, 6:10*100), rep, 15)), 
                                test1_stats)))
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
    scale_size_continuous(name = "standard deviation\nof p-value") +
    scale_x_discrete("Base number of molecules") +
    scale_y_discrete("Added number of molecules")
}


test_dpcr <- function(size, times = 4000, cl) 
  pblapply(1L:times, function(dummy_variable) {
    parSapply(cl, c(0:5*10, 2:10*50, 6:10*100), function(base_number) {
      sapply(c(1, 1:5*10, 2:10*50), function(added_molecules) {
        dat <- dpcR::sim_adpcr(m = base_number, n = size, times = 10000, pos_sums = FALSE, 
                               n_panels = 1) 
        dat2 <- dat
        ids <- sample(1L:size, added_molecules)
        dat2[ids, ] <- dat2[ids, ] + 1
        #sample added below to reshuffle positive wells
        counts_data <- data.frame(experiment = c(rep("baseline", size), rep("modyfied", size)), 
                                  counts = as.numeric(c(as.vector(dat) > 0, 
                                                        sample(as.vector(dat2)) > 0)))
        fit <- glm(counts ~ experiment + 0, data = counts_data, family = quasipoisson)
        multcomp:::summary.glht(multcomp::glht(fit, linfct = mcp(experiment = "Tukey")))[["test"]][["pvalues"]]
      })
    })
  })

test_dpcr_slow <- function(size, times = 4000) 
  pblapply(1L:times, function(dummy_variable) {
    sapply(c(0:5*10, 2:10*50, 6:10*100), function(base_number) {
      sapply(c(1, 1:5*10, 2:10*50), function(added_molecules) {
        dat <- sim_adpcr(m = base_number, n = size, times = 10000, pos_sums = FALSE, 
                         n_panels = 1) 
        dat2 <- dat
        ids <- sample(1L:size, added_molecules)
        dat2[ids, ] <- dat2[ids, ] + 1
        #sample added below to reshuffle positive wells
        counts_data <- data.frame(experiment = c(rep("baseline", size), rep("modyfied", size)), 
                                  counts = as.numeric(c(as.vector(dat) > 0, 
                                                        sample(as.vector(dat2)) > 0)))
        fit <- glm(counts ~ experiment + 0, data = counts_data, family = quasipoisson)
        summary(glht(fit, linfct = mcp(experiment = "Tukey")))[["test"]][["pvalues"]]
      })
    })
  })


test2 <- test_dpcr(2000, 40)
test3 <- test_dpcr(5000, 40)


cl <- makeCluster(6)
clusterEvalQ(cl, {
  library(dpcR)
  library(multcomp)
  NULL
})
#system.time(tmp <- test_dpcr_slow(1000, 4))
system.time(tmp <- test_dpcr(5000, 4, cl))
stopCluster(cl)

plot_dpcrtest(test1)
plot_dpcrtest(test2)
plot_dpcrtest(test3)

save(test1, test2, test3, file = "dpcrposter_data.RData")

load("dpcrposter_data.RData")
