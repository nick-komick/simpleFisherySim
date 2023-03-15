library(dplyr)
library(furrr)
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
#' @param pnl Proportion not legal size
#' @param cohort_size Size of the cohort
#'
#' @return A data frame with `cohort_size` rows
#' @export
#'
createCohort <- function(adclip_rate = 0.5,
                         pnl = 0.2,
                         cohort_size = 10000L) {
  cohort_df <-
    tibble(fish_number = seq_len(cohort_size),
           adclip_event = runif(cohort_size),
           legal_event = runif(cohort_size),
           mortality = FALSE) |>
    mutate(is_clipped = adclip_event <= adclip_rate,
           is_legal = legal_event > pnl) |>
    select(-adclip_event,
           -legal_event)

  return(cohort_df)
}

#' (Re-)Select a random cohort fish
#'
#' Select a random cohort fish, if it is dead select a different one
#' Do this until you find a fish or you have made as many tries as there are
#' fish in the cohort.
#'
#' @param cohort_df Cohort data frame
#' @param potential_fish A random number to try as first fish
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
                                 legal_release_mort_rate = 0.15,
                                 nonlegal_release_mort_rate = 0.15,
                                 fishery_adclip_rate = 0.6,
                                 fishery_pnl = 0.2) {

  event_seq <- seq_len(fishery_catch * 10)
  event_len <- length(event_seq)

  total_kept_catch <- 0L
  outcomes <- rep(NA_character_, event_len)
  mortalities <- rep(NA, event_len)
  fish_numbers <- rep(NA_integer_, event_len)

  potential_fish <- as.integer(runif(event_len, min = 1, max = nrow(cohort_df)))
  encounter_event <- runif(event_len)
  drop_event <- runif(event_len)
  rel_event <- runif(event_len)
  mort_event <- runif(event_len)
  clip_event <- runif(event_len)
  legal_event <- runif(event_len)
  is_clipped <- rep(NA, event_len)
  is_legal <- rep(NA, event_len)

  for(event_id in event_seq) {
    if(total_kept_catch  >= fishery_catch) {
      #Total allowable catch is reach, so no more fishing events to process
      break
    }

    if(encounter_event[event_id] <= encounter_rate) {
      #Encountered fish is from the provided cohort
      cohort_fish_number <- selectFish(cohort_df,
                                       potential_fish[event_id])

      if(cohort_fish_number > nrow(cohort_df)) {
        stop("Bad fish number")
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

      #If the fish was from our model cohort, use the fish legal status
      is_legal[event_id] <- cohort_df$is_legal[cohort_fish_number]
    } else {
      #If the fish was NOT from our model cohort, use randomly assign clip
      # based on the default fishery adclip rate
      is_clipped[event_id]<- clip_event[event_id] <= fishery_adclip_rate


      #If the fish was NOT from our model cohort, use randomly assign
      # vulnerability status based on the default fishery vulnerability rate
      is_legal[event_id]<- legal_event[event_id] >= fishery_pnl
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
    } else {
      fish_rel_rate <- 0
      if(is_legal[event_id] == FALSE) {
        #If the fish is not legal, then 100% chance of release
        outcomes[event_id] <- Released

        if(mort_event[event_id] <= nonlegal_release_mort_rate) {
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
        #For legal fish, select the release rate based on clip status
        fish_rel_rate <- ifelse(is_clipped[event_id],
                                clipped_release_rate,
                                unclipped_release_rate)

        if(rel_event[event_id] <= fish_rel_rate) {
          #The fish was released
          outcomes[event_id] <- Released
          if(mort_event[event_id] <= legal_release_mort_rate) {
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
          total_kept_catch <- total_kept_catch + 1L
          outcomes[event_id] <- Kept
          mortalities[event_id] <- TRUE
          if(!is.na(cohort_fish_number)) {
            #Kill the fish in the cohort data frame
            cohort_df$mortality[cohort_fish_number] <- TRUE
          }
        }
      }
    }
  }

  if(!is.na(outcomes[length(outcomes)])) {
    stop("Ran out of fishing events :-(")
  }

  fishery_df <-
    tibble(event_id = event_seq,
           potential_fish = potential_fish,
           encounter_event = encounter_event,
           drop_event = drop_event,
           rel_event = rel_event,
           mort_event = mort_event,
           clip_event = clip_event,
           is_clipped  = is_clipped,
           is_legal = is_legal,
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
    group_by(outcome, is_cohort, is_clipped, is_legal) |>
    summarize(mortality = sum(mortality),
              events = n(),
              .groups = "drop")

  return(total_catch)
}

