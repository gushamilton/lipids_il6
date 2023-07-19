pacman::p_load(tidyverse, vroom, TwoSampleMR, ggforestplot)

dat <- vroom("harmonised_exposure_outcome.dat")
res <-mr(dat, method = "mr_ivw")
res %>%
  mutate(outcome =  sub("\\|\\|.*", "", outcome)) %>%
  arrange(b) %>%
  forestplot(name = outcome,
             colour = exposure,
             estimate = b,
             se = se) +
    theme_bw() +
    theme(legend.position = "bottom") +
  xlab("SD change in each outcome for unit change in CRP")
