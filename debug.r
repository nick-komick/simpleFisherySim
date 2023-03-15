cohort_abundance <- 10000L
cohort_encounter_rate <- 0.5
cohort_adclip_rate <- 0.5
cohort_pnl <- 0.5

release_mort_rate <- 0.15
drop_off_rate <- 0.05

#background cohorts within the fishery
fishery_catch <- 10000
fishery_adclip_rate <- 0.4
fishery_pnl <- 0.5
fishery_unclip_release_rate <- 0.8
fishery_clip_release_rate <- 0.2

sim_result <-
  createCohort(cohort_size = cohort_abundance,
               adclip_rate = cohort_adclip_rate,
               pnl = cohort_pnl) |>
  sequencialFisherySim(fishery_catch = fishery_catch,
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
                  terminal_unmark_cohort = terminal_unmark_cohort)