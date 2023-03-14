library(dplyr)


## Changing cohort size
lambda_low_df <-
  cohort_df <- createCohort(cohort_size = 10000L) |>
  fisherySim(encounter_rate = .5) |>
  mutate(time_step = row_number(),
         clipped_cohort = sum(is_clipped == TRUE),
         clipped_mort_total = cumsum(if_else(is_clipped == TRUE & mortality, 1, 0)),
         unclipped_cohort = sum(is_clipped == FALSE),
         unclipped_mort_total = cumsum(if_else(is_clipped == FALSE & mortality, 1, 0)),
         lambda = (unclipped_cohort - unclipped_mort_total)/(clipped_cohort - clipped_mort_total)) |>
  filter(mortality) |>
  slice_sample(n=1000) |>
  mutate(encounter_rate = 0.5) |>
  arrange(time_step) |>
  mutate(time_step = row_number())

lambda_med_df <-
  cohort_df <- createCohort(cohort_size = 100000L) |>
  fisherySim(encounter_rate = .5) |>
  mutate(time_step = row_number(),
         clipped_cohort = sum(is_clipped == TRUE),
         clipped_mort_total = cumsum(if_else(is_clipped == TRUE & mortality, 1, 0)),
         unclipped_cohort = sum(is_clipped == FALSE),
         unclipped_mort_total = cumsum(if_else(is_clipped == FALSE & mortality, 1, 0)),
         lambda = (unclipped_cohort - unclipped_mort_total)/(clipped_cohort - clipped_mort_total)) |>
  filter(mortality) |>
  slice_sample(n=1000) |>
  mutate(encounter_rate = 0.5) |>
  arrange(time_step) |>
  mutate(time_step = row_number())


lambda_high_df <-
  cohort_df <- createCohort(cohort_size = 1000000L) |>
  fisherySim(encounter_rate = .5) |>
  mutate(time_step = row_number(),
         clipped_cohort = sum(is_clipped == TRUE),
         clipped_mort_total = cumsum(if_else(is_clipped == TRUE & mortality, 1, 0)),
         unclipped_cohort = sum(is_clipped == FALSE),
         unclipped_mort_total = cumsum(if_else(is_clipped == FALSE & mortality, 1, 0)),
         lambda = (unclipped_cohort - unclipped_mort_total)/(clipped_cohort - clipped_mort_total)) |>
  filter(mortality) |>
  slice_sample(n=1000) |>
  mutate(encounter_rate = 0.5) |>
  arrange(time_step) |>
  mutate(time_step = row_number())

lambda_df <-
  lambda_low_df |>
  bind_rows(lambda_med_df) |>
  bind_rows(lambda_high_df)


ggplot(lambda_df, aes(time_step, lambda, colour = clipped_cohort)) +
  geom_point()

## Changing encounter rate
lambda_low_df <-
  cohort_df <- createCohort(cohort_size = 100000L) |>
  fisherySim(encounter_rate = .3) |>
  mutate(time_step = row_number(),
         clipped_cohort = sum(is_clipped == TRUE),
         clipped_mort_total = cumsum(if_else(is_clipped == TRUE & mortality, 1, 0)),
         unclipped_cohort = sum(is_clipped == FALSE),
         unclipped_mort_total = cumsum(if_else(is_clipped == FALSE & mortality, 1, 0)),
         lambda = (unclipped_cohort - unclipped_mort_total)/(clipped_cohort - clipped_mort_total)) |>
  filter(mortality) |>
  slice_sample(n=10000) |>
  mutate(encounter_rate = 0.3) |>
  arrange(time_step) |>
  mutate(time_step = row_number())

lambda_med_df <-
  cohort_df <- createCohort(cohort_size = 100000L) |>
  fisherySim(encounter_rate = .5) |>
  mutate(time_step = row_number(),
         clipped_cohort = sum(is_clipped == TRUE),
         clipped_mort_total = cumsum(if_else(is_clipped == TRUE & mortality, 1, 0)),
         unclipped_cohort = sum(is_clipped == FALSE),
         unclipped_mort_total = cumsum(if_else(is_clipped == FALSE & mortality, 1, 0)),
         lambda = (unclipped_cohort - unclipped_mort_total)/(clipped_cohort - clipped_mort_total)) |>
  filter(mortality) |>
  slice_sample(n=10000) |>
  mutate(encounter_rate = 0.5) |>
  arrange(time_step) |>
  mutate(time_step = row_number())


lambda_high_df <-
  cohort_df <- createCohort(cohort_size = 100000L) |>
  fisherySim(encounter_rate = .7) |>
  mutate(time_step = row_number(),
         clipped_cohort = sum(is_clipped == TRUE),
         clipped_mort_total = cumsum(if_else(is_clipped == TRUE & mortality, 1, 0)),
         unclipped_cohort = sum(is_clipped == FALSE),
         unclipped_mort_total = cumsum(if_else(is_clipped == FALSE & mortality, 1, 0)),
         lambda = (unclipped_cohort - unclipped_mort_total)/(clipped_cohort - clipped_mort_total)) |>
  filter(mortality) |>
  slice_sample(n=10000) |>
  mutate(encounter_rate = 0.7) |>
  arrange(time_step) |>
  mutate(time_step = row_number())

lambda_df <-
  lambda_low_df |>
  bind_rows(lambda_med_df) |>
  bind_rows(lambda_high_df) #|>
#  bind_rows(tibble(time_step=rep(5000, 3), encouter_rate= "medium",
#                   lambda=c(log((exp(max(lambda_low_df$lambda) - min(lambda_low_df$lambda))/2) + min(lambda_low_df$lambda)),
#                            log((exp(max(lambda_med_df$lambda) - min(lambda_med_df$lambda))/2) + min(lambda_med_df$lambda)),
#                            log((exp(max(lambda_high_df$lambda) - min(lambda_high_df$lambda))/2) + min(lambda_high_df$lambda)))))


ggplot(lambda_df, aes(time_step, lambda, colour = encounter_rate)) +
  geom_point()
