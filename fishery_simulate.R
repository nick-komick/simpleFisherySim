library(dplyr)
library(furrr)
library(ggplot2)

NoEncounter <- "NO ENCOUNTER"
DropOff <- "DROP-OFF"
Released <- "RELEASED"
Kept <- "KEPT"


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

#' Create Cohort Data Frame
#'
#' Create a data frame of random fish based on the
#' cohort size and the adipose clip rate
#'
#' @param adclip_rate Adipose clip rate of the cohort (related to lambda)
#' @param cohort_size Size of the cohort
#'
#' @return A data frame with `cohort_size` rows
#' @export
#'
createCohort <- function(adclip_rate = 0.5,
                         cohort_size = 10000L) {
  cohort_df <-
    tibble(fish_number = seq_len(cohort_size),
           adclip_event = runif(cohort_size),
           mortality = FALSE) |>
    mutate(is_clipped = adclip_event <= adclip_rate) |>
    select(-adclip_event)

  return(cohort_df)
}

#' (Re-)Select a random cohort fish
#'
#' Select a random cohort fish, if it is dead select a different one
#' Do this until you find a fish or you have made as many tries as there are
#' fish in the cohort.
#'
#' @param cohort_df Cohort data frame
#' @param fish_event A random number to try as first fish
#'
#' @return Fish/row number of a fish that is encountered in a fishery
#'
selectFish <- function(cohort_df, potential_fish) {
  find <- 1
  fish_number <- potential_fish

  while(cohort_df$mortality[fish_number] == TRUE) {
    fish_number <- as.integer(runif(1, min = 1, max = nrow(cohort_df)))

    if(find > nrow(cohort_df)) {
      stop("Ran out of fish")
    }
    find <- find + 1
  }

  return(fish_number)
}

#' Sequential Fishery Simulation
#'
#' Sequentially simulate a fishery up to approximately a certain number of catch
#' (ie. kept or released) events.  Events may be slightly more or less to accommodate
#' the random probability of drop-off events.  If a fish from the cohort dies, it is
#' no longer available for selection in the fishery
#'
#' @param cohort_df Cohort data frame
#' @param fish_event A random number to try as first fish
#'
#' @return Fish/row number of a fish that is encountered in a fishery
#'
sequencialFisherySim <- function(cohort_df,
                                 fishery_catch = 1000L,
                                 encounter_rate = 0.7,
                                 drop_off_rate = .05,
                                 unclipped_release_rate = .9,
                                 clipped_release_rate = .1,
                                 drop_mort_rate = 1.0,
                                 release_mort_rate = 0.15,
                                 fishery_adclip_rate = 0.6) {

  event_seq <- seq_len(fishery_catch * 2)
  event_len <- length(event_seq)


  total_kept_catch <- 0L
  outcomes <- rep(NA_character_, event_len)
  mortalities <- rep(NA, event_len)
  fish_numbers <- rep(NA_integer_, event_len)

  potential_fish <- as.integer(runif(event_len, min = 1, max = nrow(cohort_df)))
  encouter_event <- runif(event_len)
  drop_event <- runif(event_len)
  rel_event <- runif(event_len)
  mort_event <- runif(event_len)
  clip_event <- runif(event_len)
  is_clipped <- rep(NA, event_len)

  for(event_id in event_seq) {
    if(total_kept_catch  >= fishery_catch) {
      break
    }

    if(encouter_event[event_id] <= encounter_rate) {
      #Encountered fish is from the provided cohort
      cohort_fish_number <- selectFish(cohort_df,
                                       potential_fish[event_id])

      if(cohort_fish_number > nrow(cohort_df)) {
        stop()
      }

    } else {
      #Encountered fish is not from the provided cohort
      cohort_fish_number <- NA_integer_
    }


    fish_numbers[event_id] <- cohort_fish_number

    #The fish was caught, we need to figure out if it was kept or released
    if(!is.na(cohort_fish_number)) {
      #If the fish was from our model cohort, use the fish adclip status
      is_clipped[event_id] <- cohort_df$is_clipped[cohort_fish_number]
    } else {
      #If the fish was NOT from our model cohort, use randomly assign clip
      # based on the default fishery adclip rate
      is_clipped[event_id]<- clip_event[event_id] <= fishery_adclip_rate
    }


    if(drop_event[event_id] <= drop_off_rate) {
      #fish dropped off
      outcomes[event_id] <- DropOff
      if(mort_event[event_id] <= drop_mort_rate) {
        #Fish died after dropping off
        mortalities[event_id] <- TRUE
        if(!is.na(cohort_fish_number)) {
          #Kill the fish in the cohort data frame
          cohort_df$mortality[cohort_fish_number] <- TRUE
        }
      } else {
        mortalities[event_id] <- FALSE
      }
    }else {
      total_kept_catch <- total_kept_catch + 1L
      #Select the release rate based on clip status
      fish_rel_rate <- NA_real_
      if(is_clipped[event_id]) {
        fish_rel_rate <- clipped_release_rate
      } else {
        fish_rel_rate <- unclipped_release_rate
      }

      if(rel_event[event_id] <= fish_rel_rate) {
        #The fish was released
        outcomes[event_id] <- Released
        if(mort_event[event_id] <= release_mort_rate) {
          #The fish died immediately after being released
          mortalities[event_id] <- TRUE
          if(!is.na(cohort_fish_number)) {
            #Kill the fish in the cohort data frame
            cohort_df$mortality[cohort_fish_number] <- TRUE
          }
        } else {
          mortalities[event_id] <- FALSE
        }
      } else {
        #The fish was kept and died
        outcomes[event_id] <- Kept
        mortalities[event_id] <- TRUE
        if(!is.na(cohort_fish_number)) {
          #Kill the fish in the cohort data frame
          cohort_df$mortality[cohort_fish_number] <- TRUE
        }
      }
    }
  }

  fishery_df <-
    tibble(event_id = event_seq,
           potential_fish = potential_fish,
           encouter_event = encouter_event,
           drop_event = drop_event,
           rel_event = rel_event,
           mort_event = mort_event,
           clip_event = clip_event,
           is_clipped  = is_clipped,
           outcome = outcomes,
           mortality = mortalities,
           fish_number = fish_numbers) |>
    filter(!is.na(outcome))

  return(list(fishery_df = fishery_df, cohort_df = cohort_df))
}

#' Summarize Fishery Simulation
#'
#' Summarize the individual fishing events from a simulated fishery
#'
#' @param cohort_df Cohort data frame
#' @param fish_event A random number to try as first fish
#'
#' @return Fish/row number of a fish that is encountered in a fishery
#'
summarizeFishery <- function(fishery_df) {
  total_catch <-
    fishery_df |>
    mutate(is_cohort = if_else(is.na(fish_number), FALSE, TRUE)) |>
    group_by(outcome, is_cohort, is_clipped) |>
    summarize(mortality = sum(mortality),
              events = n(),
              .groups = "drop")

  return(total_catch)
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


singleFisheryRun <- function(catch) {
  cohort_abundance <- 2000L
  encounter_rate <- 0.1
  release_mort_rate <- 0.15
  drop_mort_rate <- 0.05
  lambda <- 1
  adclip_rate <- 0.5

  sim_result <-
    createCohort(cohort_size = cohort_abundance,
                 adclip_rate = adclip_rate) |>
    sequencialFisherySim(fishery_catch = catch,
                         encounter_rate = encounter_rate,
                         drop_off_rate = drop_mort_rate,
                         unclipped_release_rate = .9,
                         clipped_release_rate = .1,
                         drop_mort_rate = 1.0,
                         release_mort_rate = release_mort_rate,
                         fishery_adclip_rate = 0.2)

  fishery_summary <- summarizeFishery(sim_result$fishery_df)

  model_mortality_df <-
    fishery_summary |>
    filter(is_cohort == TRUE) |>
    group_by(is_clipped) |>
    summarize(model_mortality = sum(mortality), .groups = "drop")


  post_cohort <-
    sim_result$cohort_df |>
    filter(mortality == FALSE) |>
    count(is_clipped)

  sim_post_lambda <- post_cohort$n[post_cohort$is_clipped == FALSE] / post_cohort$n[post_cohort$is_clipped == TRUE]

  mfa_post_lambda <- mfaLambda(cohort_abundance,
                               lambda,
                               release_mort_rate,
                               drop_mort_rate,
                               fishery_summary)

  mid_lambda <- mfaLambda(cohort_abundance,
                          lambda,
                          release_mort_rate,
                          drop_mort_rate,
                          fishery_summary,
                          TRUE)


  mfa_mortality_df <-
    mfaMethod(mid_lambda, release_mort_rate, drop_mort_rate, fishery_summary) |>
    full_join(model_mortality_df, by="is_clipped") |>
    mutate(catch = catch,
           bias = (abs(model_mortality) - abs(mfa_mortality)) / abs(model_mortality),
           sim_post_lambda = sim_post_lambda,
           mfa_post_lambda = mfa_post_lambda)

  return(mfa_mortality_df)
}

plan(multisession, workers = 5)


bias_estimate <-
  furrr::future_map_dfr(runif(2000, 1000, 10000),singleFisheryRun, .options = furrr_options(seed = T))

ggplot(bias_estimate, aes(catch, bias, colour = is_clipped)) +
  geom_point()

lambda_bias <-
  bias_estimate |>
  filter(is_clipped == FALSE) |>
  mutate(lambda_bias = (mfa_post_lambda- sim_post_lambda)/ sim_post_lambda)


bias_estimate  |> filter(!is_clipped) |> pull(bias) |> mean()

ggplot(bias_estimate, aes(catch, bias, color = is_clipped)) +
  geom_point() +
  geom_smooth(method = "loess")


ggplot(lambda_bias, aes(catch, lambda_bias)) +
  geom_point() +
  geom_smooth(method = "loess")

#fishery_summary <-
#  lapply(seq_len(10), function(.) {
#    sim_result <-
#      createCohort(cohort_size = 100000L) |>
#      sequencialFisherySim(fishery_catch = fishery_catch,
#                           encounter_rate = encounter_rate)
#    return(summarizeFishery(sim_result$fishery_df))
#  }) %>%
#  bind_rows() %>%
#  group_by(outcome, is_cohort, is_clipped) %>%
#  summarize(mortality = mean(mortality), .groups="drop")
#

#
#cohort_df <- sim_result$cohort_df
#fishery_df <- sim_result$fishery_df
#
#cohort_fishery_df <-
#  fishery_df |>
#  select(-mortality) |>
#  inner_join(cohort_df, by="fish_number") |>
#  mutate(clipped_cohort = sum(cohort_df$is_clipped == TRUE),
#         clipped_mort_total = cumsum(if_else(is_clipped == TRUE & mortality, 1, 0)),
#         unclipped_cohort = sum(cohort_df$is_clipped == FALSE),
#         unclipped_mort_total = cumsum(if_else(is_clipped == FALSE & mortality, 1, 0)),
#         lambda = (unclipped_cohort - unclipped_mort_total)/(clipped_cohort - clipped_mort_total))
#
#
#ggplot(cohort_fishery_df, aes(event_id, lambda)) +
#  geom_point() +
#  expand_limits(x = 0, y = 1)#