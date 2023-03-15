
instantaneousMfa <- function(kept_mark,
                             release_mark,
                             kept_unmark,
                             release_unmark,
                             kept_mark_cohort,
                             escapement_mark_cohort,
                             release_mort_rate,
                             drop_off_rate,
                             escapement_unmark_cohort = NULL,
                             terminal_unmark_cohort = NULL) {

  if(is.null(escapement_unmark_cohort) && is.null(terminal_unmark_cohort) ) {
    stop("Either the unmarked cohort escapement or terminal unmarked cohort must be provided")
  }

  if(!is.null(escapement_unmark_cohort) && !is.null(terminal_unmark_cohort) ) {
    stop("Only the the unmarked cohort escapement or terminal unmarked cohort can be provided")
  }

  drop_off_ratio <- drop_off_rate / (1-drop_off_rate)

  catch_mark <- kept_mark + release_mark

  release_mark_cohort <- release_mark * kept_mark_cohort / kept_mark

  catch_mark_cohort <- release_mark_cohort + kept_mark_cohort

  drop_mark_cohort <- catch_mark_cohort * drop_off_ratio

  terminal_mark_cohort <- kept_mark_cohort +
    release_mort_rate * release_mark_cohort +
    escapement_mark_cohort + drop_mark_cohort

  A_m <- 1 - escapement_mark_cohort / terminal_mark_cohort

  Z_m <- -1 * log(1 - A_m)

  mu_u_k <- kept_unmark / (kept_unmark + release_unmark) *
    catch_mark_cohort / terminal_mark_cohort

  F_u_k <- mu_u_k * Z_m / A_m

  mu_u_r <- release_mort_rate * release_unmark / (kept_unmark + release_unmark) *
    catch_mark_cohort / terminal_mark_cohort

  F_u_r <- mu_u_r * Z_m / A_m

  mu_m_d <- drop_mark_cohort / terminal_mark_cohort

  F_u_d <- mu_m_d * Z_m / A_m

  Z_u <- F_u_r + F_u_k + F_u_d

  if(is.null(terminal_unmark_cohort)) {
    terminal_unmark_cohort <- escapement_unmark_cohort / exp(-Z_u)
  } else {
    escapement_unmark_cohort <- terminal_unmark_cohort * exp(-Z_u)
  }
  return(tibble(post_fishery_mark_cohort = escapement_mark_cohort,
                pre_fishery_mark_cohort = terminal_mark_cohort,
                post_fishery_unmark_cohort = escapement_unmark_cohort,
                pre_fishery_unmark_cohort = terminal_unmark_cohort))
}





instMfaNonLegal <- function(kept_mark,
                            legal_release_mark,
                            kept_unmark,
                            legal_release_unmark,
                            kept_mark_cohort,
                            escapement_mark_cohort,
                            legal_release_mort_rate,
                            nonlegal_release_mort_rate,
                            prop_nonlegal,
                            drop_off_rate,
                            escapement_unmark_cohort = NULL,
                            terminal_unmark_cohort = NULL) {

  if(is.null(escapement_unmark_cohort) && is.null(terminal_unmark_cohort) ) {
    stop("Either the unmarked cohort escapement or terminal unmarked cohort must be provided")
  }

  if(!is.null(escapement_unmark_cohort) && !is.null(terminal_unmark_cohort) ) {
    stop("Only the the unmarked cohort escapement or terminal unmarked cohort can be provided")
  }

  drop_off_ratio <- drop_off_rate / (1-drop_off_rate)

  legal_catch_mark <- kept_mark + legal_release_mark
  legal_catch_unmark <- kept_unmark + legal_release_unmark

  legal_release_mark_cohort <- legal_release_mark * kept_mark_cohort / kept_mark
  legal_catch_mark_cohort <- legal_release_mark_cohort + kept_mark_cohort

  nonlegal_rel_mark_cohort <- legal_catch_mark_cohort * prop_nonlegal / (1-prop_nonlegal)



  drop_mark_cohort <- (legal_catch_mark_cohort + nonlegal_rel_mark_cohort) * drop_off_ratio

  terminal_mark_cohort <- kept_mark_cohort +
    legal_release_mort_rate * legal_release_mark_cohort +
    nonlegal_release_mort_rate * nonlegal_rel_mark_cohort +
    escapement_mark_cohort + drop_mark_cohort

  A_m <- 1 - escapement_mark_cohort / terminal_mark_cohort

  Z_m <- -1 * log(1 - A_m)

  mu_u_k <- kept_unmark / legal_catch_unmark * legal_catch_mark_cohort / terminal_mark_cohort

  F_u_k <- mu_u_k * Z_m / A_m

  #Legal Release of Unmarked Cohort
  mu_u_lr <- legal_release_mort_rate * legal_release_unmark / legal_catch_unmark *
    legal_catch_mark_cohort / terminal_mark_cohort

  F_u_lr <- mu_u_lr * Z_m / A_m

  #Non-Legal Release of Unmarked Cohort
  mu_u_nlr <- nonlegal_release_mort_rate * legal_catch_unmark * prop_nonlegal / (1- prop_nonlegal) / legal_catch_unmark *
    legal_catch_mark_cohort / terminal_mark_cohort

  F_u_nlr <- mu_u_nlr * Z_m / A_m

  #Dropoff Mortality Rate of Cohort
  mu_m_d <- drop_mark_cohort / terminal_mark_cohort

  F_u_d <- mu_m_d * Z_m / A_m

  Z_u <- F_u_nlr + F_u_lr + F_u_k + F_u_d

  if(is.null(terminal_unmark_cohort)) {
    terminal_unmark_cohort <- escapement_unmark_cohort / exp(-Z_u)
  } else {
    escapement_unmark_cohort <- terminal_unmark_cohort * exp(-Z_u)
  }
  return(tibble(post_fishery_mark_cohort = escapement_mark_cohort,
                pre_fishery_mark_cohort = terminal_mark_cohort,
                post_fishery_unmark_cohort = escapement_unmark_cohort,
                pre_fishery_unmark_cohort = terminal_unmark_cohort))
}



kept_mark <- 100
legal_release_mark <- 10
kept_unmark <- 1000
legal_release_unmark <- 1000
kept_mark_cohort <- 36.363636
escapement_mark_cohort <- 209.09091
legal_release_mort_rate <- 0.2
drop_off_rate <- 0.05
escapement_unmark_cohort <- 535
terminal_unmark_cohort <- 606


instantaneousMfa(kept_mark,
                 legal_release_mark,
                 kept_unmark,
                 legal_release_unmark,
                 kept_mark_cohort,
                 escapement_mark_cohort,
                 legal_release_mort_rate,
                 drop_off_rate,
                 escapement_unmark_cohort = escapement_unmark_cohort)

instantaneousMfa(kept_mark,
                 legal_release_mark,
                 kept_unmark,
                 legal_release_unmark,
                 kept_mark_cohort,
                 escapement_mark_cohort,
                 legal_release_mort_rate,
                 drop_off_rate,
                 terminal_unmark_cohort = terminal_unmark_cohort)

nonlegal_release_mort_rate <- .2
prop_nonlegal <- .2

instMfaNonLegal(kept_mark,
                legal_release_mark,
                kept_unmark,
                legal_release_unmark,
                kept_mark_cohort,
                escapement_mark_cohort,
                legal_release_mort_rate,
                nonlegal_release_mort_rate,
                prop_nonlegal,
                drop_off_rate,
                escapement_unmark_cohort = escapement_unmark_cohort)

instMfaNonLegal(kept_mark,
                legal_release_mark,
                kept_unmark,
                legal_release_unmark,
                kept_mark_cohort,
                escapement_mark_cohort,
                legal_release_mort_rate,
                nonlegal_release_mort_rate,
                prop_nonlegal,
                drop_off_rate,
                terminal_unmark_cohort = terminal_unmark_cohort)

