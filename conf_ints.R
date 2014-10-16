library(pbapply)
library(dpcR)

c(1L:5*3, 4L:10*5, 6L:10*10)

get_ints <- function(number_of_exps) {
  adpcr1 <- sim_adpcr(m = 10, n = 765, times = 10000, pos_sums = FALSE, 
                      n_panels = number_of_exps)
  adpcr2 <- sim_adpcr(m = 60, n = 765, times = 10000, pos_sums = FALSE, 
                      n_panels = number_of_exps)
  dat <- bind_dpcr(adpcr1, adpcr2)
  modeltab <- slot(test_counts(dat), "group_coef")
  dpcrtab <- summary(dat, print = FALSE)[["summary"]]
  list(ints = cbind(dpcrtab[seq(1, nrow(dpcrtab), by = 2), 3:5], 
        dpcrtab[seq(1, nrow(dpcrtab), by = 2) + 1, 3:5],
        modeltab[, -1]),
       real_m = colSums(dat))
}

ints5 <- get_ints(5)

save(ints5, file = "ints_plot.RData")
