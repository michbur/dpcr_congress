test_dpcr <- function(size, times = 4000, cl) 
  pblapply(1L:times, function(dummy_variable) {
    parSapply(cl, c(0:5*10, 2:10*50, 6:10*100), function(base_number) {
      sapply(c(1, 5, 1:9*10, 2:5*50), function(added_molecules) {
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

library(pbapply)
library(parallel)
cl <- makeCluster(6)
clusterEvalQ(cl, {
  library(dpcR)
  library(multcomp)
  NULL
})
#system.time(tmp <- test_dpcr_slow(1000, 4))
mc1000 <- test_dpcr(1000, 1000, cl)
mc5000 <- test_dpcr(5000, 1000, cl)
mc10000 <- test_dpcr(10000, 1000, cl)
stopCluster(cl)

save(mc1000, mc5000, mc10000, file = "fin_dpcrposter_rawdata.RData")

