library(dplyr)
library(ggplot2)

NoEncounter <- "NO ENCOUNTER"
DropOff <- "DROP-OFF"
Released <- "RELEASED"
Kept <- "KEPT"

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
                         cohort_size = 1000000L) {
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
                                 fishery_adclip_rate = 0.5) {

  event_seq <- seq_len(fishery_catch * (1 + drop_off_rate))
  event_len <- length(event_seq)

  outcomes <- NULL
  mortalities <- NULL
  fish_numbers <- NULL

  potential_fish <- as.integer(runif(event_len, min = 1, max = nrow(cohort_df)))
  encouter_event <- runif(event_len)
  drop_event <- runif(event_len)
  rel_event <- runif(event_len)
  mort_event <- runif(event_len)
  clip_event <- runif(event_len)

  for(event_id in event_seq) {
    if(encouter_event[event_id] <= encounter_rate) {
      #Encountered fish is from the provided cohort
      cohort_fish_number <- selectFish(cohort_df,
                                       potential_fish[event_id])
    } else {
      #Encountered fish is not from the provided cohort
      cohort_fish_number <- NA_integer_
    }

    fish_numbers <- c(fish_numbers, cohort_fish_number)

    if(drop_event[event_id] <= drop_off_rate) {
      #fish dropped off
      outcomes <- c(outcomes, DropOff)
      if(mort_event[event_id] <= drop_mort_rate) {
        #Fish died after dropping off
        mortalities <- c(mortalities, TRUE)
        if(!is.na(cohort_fish_number)) {
          #Kill the fish in the cohort data frame
          cohort_df$mortality[cohort_fish_number] <- TRUE
        }
      } else {
        mortalities <- c(mortalities, FALSE)
      }
    }else {
      #The fish was caught, we need to figure out if it was kept or released

      fish_rel_rate <- NA_real_
      is_clipped <- NA
      if(!is.na(cohort_fish_number)) {
        #If the fish was from our model cohort, use the fish adclip status
        is_clipped <- cohort_df$is_clipped[cohort_fish_number]
      } else {
        #If the fish was NOT from our model cohort, use randomly assign clip
        # based on the default fishery adclip rate
        is_clipped <- clip_event[event_id] <= fishery_adclip_rate
      }

      #Select the release rate based on clip status
      if(is_clipped) {
        fish_rel_rate <- clipped_release_rate
      } else {
        fish_rel_rate <- unclipped_release_rate
      }

      if(rel_event[event_id] <= fish_rel_rate) {
        #The fish was released
        outcomes <- c(outcomes, Released)
        if(mort_event[event_id] <= release_mort_rate) {
          #The fish died immediately after being released
          mortalities <- c(mortalities, TRUE)
          if(!is.na(cohort_fish_number)) {
            #Kill the fish in the cohort data frame
            cohort_df$mortality[cohort_fish_number] <- TRUE
          }
        } else {
          mortalities <- c(mortalities, FALSE)
        }
      } else {
        #The fish was kept and died
        outcomes <- c(outcomes, Kept)
        mortalities <- c(mortalities, TRUE)
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
           outcome = outcomes,
           mortality = mortalities,
           fish_number = fish_numbers)

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
    mutate(is_cohort = if_else(is.na(fish_number), "NO", "YES")) |>
    group_by(outcome, is_cohort) |>
    summarize(mortality = sum(mortality), .groups = "drop")

  return(total_catch)
}


encounter_rate <- 0.3


sim_result <-
  createCohort(cohort_size = 100000L) |>
  sequencialFisherySim(fishery_catch = fishery_catch,
                       encounter_rate = encounter_rate)


sequencialFisherySimNew(createCohort(cohort_size = 100000L),
                        fishery_catch = fishery_catch,
                        encounter_rate = encounter_rate)


fishery_sim_result <-
  lapply(seq_len(10), function(.) {
    sim_result <-
      createCohort(cohort_size = 100000L) |>
      sequencialFisherySim(fishery_catch = fishery_catch,
                           encounter_rate = encounter_rate)
    return(summarizeFishery(sim_result$fishery_df))
  }) %>%
  bind_rows() %>%
  group_by(outcome, is_cohort) %>%
  summarize(mortality = mean(mortality), .groups="drop")
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