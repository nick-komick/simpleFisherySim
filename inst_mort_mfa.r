
instantaiousMfa <- function(kept_mark,
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

  catch_mark <- kept_mark + release_mark

  release_mark_cohort <- release_mark * kept_mark_cohort / kept_mark

  catch_mark_cohort <- release_mark_cohort + kept_mark_cohort

  drop_mark_cohort <- catch_mark_cohort * drop_off_rate

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
  return(c(escapement_mark_cohort, terminal_mark_cohort, escapement_unmark_cohort, terminal_unmark_cohort))
}

instUnclip(100,
           10,
           1000,
           1000,
           36.363636,
           20.909091,
           0.2,
           0.05,
           escapement_unmark_cohort = 29.76119)

instUnclip(100,
           10,
           750,
           750,
           36.363636,
           20.909091,
           0.2,
           0.05,
           terminal_unmark_cohort = 60)
