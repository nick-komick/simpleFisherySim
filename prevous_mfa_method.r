

mfaLambda <- function(abundance, pre_lambda, release_mort_rate, dropoff_mort_rate, fishery_summary, mid_lambda = FALSE) {
  if (mid_lambda) {
    fishery_summary <-
      fishery_summary |>
      mutate(events = events / 2,
             mortality= mortality / 2)
  }

  catch_df <-
    fishery_summary |>
    filter(outcome %in% c(Kept, Released)) |>
    group_by(is_clipped) |>
    summarize(encounters = sum(events) * (1 + dropoff_mort_rate), .groups = "drop")

  marked_encounters <-
    catch_df |>
    filter(is_clipped) |>
    pull(encounters)

  mortality_rate_df <-
    fishery_summary |>
    filter(outcome %in% c(Kept, Released)) |>
    group_by(outcome, is_clipped) |>
    summarize(catch = sum(events), .groups="drop") |>
    mutate(drop_off_mortality = catch * dropoff_mort_rate,
           catch_mortality = if_else(outcome == Released, catch*release_mort_rate, as.double(catch)),
           mortality = drop_off_mortality + catch_mortality) |>
    group_by(is_clipped) |>
    summarize(mortality = sum(mortality), .groups="drop") |>
    inner_join(catch_df, by=c("is_clipped")) |>
    mutate(mortality_rate = mortality/encounters)

  marked_mort_rate <-
    mortality_rate_df |>
    filter(is_clipped) |>
    pull(mortality_rate)

  #marked_mort_rate <- 1 - exp(-1 * marked_mort_rate)

  unmarked_mort_rate <-
    mortality_rate_df |>
    filter(is_clipped == FALSE) |>
    pull(mortality_rate)

  #unmarked_mort_rate <- 1 - exp(-1 * unmarked_mort_rate)

  marked_cohort_catch <-
    fishery_summary |>
    filter(is_clipped, is_cohort, outcome %in% c(Released,Kept)) |>
    pull(events) |>
    sum()


  abundance_marked <- abundance /(pre_lambda + 1)

  mid_lambda <- (abundance_marked * pre_lambda - pre_lambda * (marked_cohort_catch * 1.05) * unmarked_mort_rate)/
    (abundance_marked - (marked_cohort_catch * 1.05) * marked_mort_rate)

  return(mid_lambda)
}


mfaMethod <- function(lambda_pre,
                      release_mort_rate,
                      drop_mort_rate,
                      fishery_summary) {
  total_catch <-
    fishery_summary |>
    filter(outcome %in% c(Kept, Released)) |>
    group_by(outcome, is_clipped) |>
    summarize(events = sum(events, na.rm=TRUE), .groups="drop")


  kept_mark_cohort <-
    fishery_summary |>
    filter(outcome %in% c(Kept),
           is_clipped == TRUE,
           is_cohort == TRUE) |>
    pull(events)

  kept_mark <-
    total_catch |>
    filter(outcome %in% c(Kept),
           is_clipped == TRUE) |>
    pull(events)


  release_mark <-
    total_catch |>
    filter(outcome %in% c(Released),
           is_clipped == TRUE) |>
    pull(events)


  kept_unmark <-
    total_catch |>
    filter(outcome %in% c(Kept),
           is_clipped == FALSE) |>
    pull(events)


  release_unmark <-
    total_catch |>
    filter(outcome %in% c(Released),
           is_clipped == FALSE) |>
    pull(events)

  release_ratio <- release_mark/ kept_mark


  mark_cohort_mortality <- kept_mark_cohort * (1 + release_mort_rate * release_ratio +
                                                 drop_mort_rate * (1 + release_ratio) )

  unmark_cohort_mortality <- lambda_pre * (kept_mark_cohort + release_ratio * kept_mark_cohort ) *
    (((kept_unmark + release_mort_rate * release_unmark) / (kept_unmark + release_unmark) ) + drop_mort_rate)


  return(tibble(is_clipped = c(TRUE, FALSE),
                mfa_mortality = c(mark_cohort_mortality, unmark_cohort_mortality)))
}

