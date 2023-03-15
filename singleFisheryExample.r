library(dplyr)
library(furrr)
library(ggplot2)
library(cowplot)


source("fishery_simulate.R")
source("inst_mort_mfa.r")

singleFisheryRun <- function(catch) {
  cohort_abundance <- 10000L
  cohort_encounter_rate <- 0.3
  cohort_adclip_rate <- 0.5
  cohort_pnl <- 0.5

  release_mort_rate <- 0.15
  drop_off_rate <- 0.05

  #background cohorts within the fishery
  fishery_adclip_rate <- 0.4
  fishery_pnl <- 0.5
  fishery_unclip_release_rate <- 0.8
  fishery_clip_release_rate <- 0.2

  sim_result <-
    createCohort(cohort_size = cohort_abundance,
                 adclip_rate = cohort_adclip_rate,
                 pnl = cohort_pnl) |>
    sequencialFisherySim(fishery_catch = catch,
                         encounter_rate = cohort_encounter_rate,
                         drop_off_rate = drop_off_rate,
                         unclipped_release_rate = fishery_unclip_release_rate,
                         clipped_release_rate = fishery_clip_release_rate,
                         drop_mort_rate = 1.0,
                         legal_release_mort_rate = release_mort_rate,
                         nonlegal_release_mort_rate = release_mort_rate,
                         fishery_adclip_rate = fishery_adclip_rate,
                         fishery_pnl = fishery_pnl)

  fishery_summary <- summarizeFishery(sim_result$fishery_df)

  kept_mark <-
    fishery_summary |>
    filter(is_clipped == TRUE, outcome == Kept) |>
    pull(events) |>
    sum()

  kept_mark_cohort <-
    fishery_summary |>
    filter(is_clipped == TRUE, outcome == Kept, is_cohort == TRUE) |>
    pull(events) |>
    sum()

  legal_release_mark <-
    fishery_summary |>
    filter(is_clipped == TRUE, outcome == Released, is_legal == TRUE) |>
    pull(events) |>
    sum()

  kept_unmark <-
    fishery_summary |>
    filter(is_clipped == FALSE, outcome == Kept) |>
    pull(events) |>
    sum()

  legal_release_unmark <-
    fishery_summary |>
    filter(is_clipped == FALSE, outcome == Released, is_legal == TRUE) |>
    pull(events) |>
    sum()

  post_mark_cohort <-
    sim_result$cohort_df |>
    filter(mortality == FALSE, is_clipped == TRUE) |>
    nrow()

  post_unmark_cohort <-
    sim_result$cohort_df |>
    filter(mortality == FALSE, is_clipped == FALSE) |>
    nrow()

  terminal_unmark_cohort <-
    sim_result$cohort_df |>
    filter(is_clipped == FALSE) |>
    nrow()


  non_pnv_est_fishery <-
    instantaneousMfa(kept_mark,
                     legal_release_mark,
                     kept_unmark,
                     legal_release_unmark,
                     kept_mark_cohort,
                     post_mark_cohort,
                     release_mort_rate,
                     drop_off_rate,
                     terminal_unmark_cohort = terminal_unmark_cohort) |>
    rename_with(.fn = ~ paste0("nopnl_", .x))

  pnv_est_fishery <-
    instMfaNonLegal(kept_mark,
                    legal_release_mark,
                    kept_unmark,
                    legal_release_unmark,
                    kept_mark_cohort,
                    post_mark_cohort,
                    release_mort_rate,
                    release_mort_rate,
                    cohort_pnl,
                    drop_off_rate,
                    terminal_unmark_cohort = terminal_unmark_cohort) |>
    rename_with(.fn = ~ paste0("pnl_", .x))

  sim_pre_fishery_mark_cohort <-
    sim_result$cohort_df |>
    filter(is_clipped == TRUE) |>
    nrow()

  sim_pre_fishery_unmark_cohort <-
    sim_result$cohort_df |>
    filter(is_clipped == FALSE) |>
    nrow()

  sim_post_fishery_mark_cohort <-
    sim_result$cohort_df |>
    filter(is_clipped == TRUE, mortality == FALSE) |>
    nrow()

  sim_post_fishery_unmark_cohort <-
    sim_result$cohort_df |>
    filter(mortality == FALSE, is_clipped == FALSE) |>
    nrow()

  sim_result <-
    tibble(sim_post_fishery_mark_cohort = sim_post_fishery_mark_cohort,
           sim_pre_fishery_mark_cohort = sim_pre_fishery_mark_cohort,
           sim_post_fishery_unmark_cohort = sim_post_fishery_unmark_cohort,
           sim_pre_fishery_unmark_cohort = sim_pre_fishery_unmark_cohort) |>
    bind_cols(non_pnv_est_fishery) |>
    bind_cols(pnv_est_fishery)


  return(sim_result)
}
singleFisheryRun(1000)
plan(multisession, workers = 5)

sim_result <-
  furrr::future_map_dfr(as.integer(runif(100, 100, 5000)),singleFisheryRun, .options = furrr_options(seed = T)) |>
  mutate(nopnl_mfa_unmark_mort = nopnl_pre_fishery_unmark_cohort - nopnl_post_fishery_unmark_cohort,
         pnl_mfa_unmark_mort = pnl_pre_fishery_unmark_cohort - pnl_post_fishery_unmark_cohort,
         nopnl_mfa_mark_mort = nopnl_pre_fishery_mark_cohort - nopnl_post_fishery_mark_cohort,
         pnl_mfa_mark_mort = pnl_pre_fishery_mark_cohort - pnl_post_fishery_mark_cohort,
         sim_unmark_mort = sim_pre_fishery_unmark_cohort - sim_post_fishery_unmark_cohort,
         sim_mark_mort = sim_pre_fishery_mark_cohort - sim_post_fishery_mark_cohort)

nonpnl_mark_plot <-
  ggplot(sim_result, aes(sim_mark_mort,
                         nopnl_mfa_mark_mort)) +
  geom_point() +
  geom_smooth(method='lm', formula= y~x) +
  geom_abline(intercept = 0, slope = 1) +
  labs(
    x = "Simulated Marked Cohort Mortalities",
    y = "MFA Marked Cohort Mortalities",
    title = "PNL Not Included"
  )


pnl_mark_plot <-
  ggplot(sim_result, aes(sim_mark_mort,
                         pnl_mfa_mark_mort)) +
  geom_point() +
  geom_smooth(method='lm', formula= y~x) +
  geom_abline(intercept = 0, slope = 1) +
  labs(
    x = "Simulated Marked Cohort Mortalities",
    y = "MFA Marked Cohort Mortalities",
    title = "PNL Included"
  )


nonpnl_unmark_plot <-
  ggplot(sim_result, aes(sim_unmark_mort,
                         nopnl_mfa_unmark_mort)) +
  geom_point() +
  geom_smooth(method='lm', formula= y~x) +
  geom_abline(intercept = 0, slope = 1) +
  labs(
    x = "Simulated Unmarked Cohort Mortalities",
    y = "MFA Unmarked Cohort Mortalities",
    title = "PNL Not Included"
  )


pnl_unmark_plot <-
  ggplot(sim_result, aes(sim_unmark_mort,
                         pnl_mfa_unmark_mort)) +
  geom_point() +
  geom_smooth(method='lm', formula= y~x) +
  geom_abline(intercept = 0, slope = 1) +
  labs(
    x = "Simulated Unmarked Cohort Mortalities",
    y = "MFA Unmarked Cohort Mortalities",
    title = "PNL Included"
  )

plot_grid(nonpnl_mark_plot, pnl_mark_plot,
          nonpnl_unmark_plot, pnl_unmark_plot, labels = "AUTO")

lm(nopnl_mfa_mark_mort ~ sim_mark_mort, sim_result)
lm(pnl_mfa_mark_mort ~ sim_mark_mort, sim_result)



