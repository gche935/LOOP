#' Trust dataset
#'
#'
#' Simulated dataset for demonstrating LOOP functions.
#'
#' @format A data frame with 1000 rows and 19 variables:
#' \describe{
#'   \item{ID}{Numeric, case ID number.}
#'   \item{trust.1}{Numeric, Average of 6 Trust items at Time 1: from 1 (strongly disagree) to 5 (strongly agree).}
#'   \item{trust.2}{Numeric, Average of 6 Trust items at Time 2: from 1 (strongly disagree) to 5 (strongly agree).}
#'   \item{trust.3}{Numeric, Average of 6 Trust items at Time 3: from 1 (strongly disagree) to 5 (strongly agree).}
#'   \item{trust.4}{Numeric, Average of 6 Trust items at Time 4: from 1 (strongly disagree) to 5 (strongly agree).}
#'   \item{trust.5}{Numeric, Average of 6 Trust items at Time 5: from 1 (strongly disagree) to 5 (strongly agree).}
#'   \item{trust.6}{Numeric, Average of 6 Trust items at Time 6: from 1 (strongly disagree) to 5 (strongly agree).}
#'   \item{lonely.1}{Numeric, Average of 6 UCLS Loneliness Scale at Time 1: from 1 (never) to 4 (always).}
#'   \item{lonely.2}{Numeric, Average of 6 UCLS Loneliness Scale at Time 2: from 1 (never) to 4 (always).}
#'   \item{lonely.3}{Numeric, Average of 6 UCLS Loneliness Scale at Time 3: from 1 (never) to 4 (always).}
#'   \item{lonely.4}{Numeric, Average of 6 UCLS Loneliness Scale at Time 4: from 1 (never) to 4 (always).}
#'   \item{lonely.5}{Numeric, Average of 6 UCLS Loneliness Scale at Time 5: from 1 (never) to 4 (always).}
#'   \item{lonely.6}{Numeric, Average of 6 UCLS Loneliness Scale at Time 6: from 1 (never) to 4 (always).}
#'   \item{lifesat.1}{Numeric, Average of 5 Satisfaction with Life Scale at Time 1: from 1 (strongly disagree) to 7 (strongly agree).}
#'   \item{lifesat.2}{Numeric, Average of 5 Satisfaction with Life Scale at Time 2: from 1 (strongly disagree) to 7 (strongly agree).}
#'   \item{lifesat.3}{Numeric, Average of 5 Satisfaction with Life Scale at Time 3: from 1 (strongly disagree) to 7 (strongly agree).}
#'   \item{lifesat.4}{Numeric, Average of 5 Satisfaction with Life Scale at Time 4: from 1 (strongly disagree) to 7 (strongly agree).}
#'   \item{lifesat.5}{Numeric, Average of 5 Satisfaction with Life Scale at Time 5: from 1 (strongly disagree) to 7 (strongly agree).}
#'   \item{lifesat.6}{Numeric, Average of 5 Satisfaction with Life Scale at Time 6: from 1 (strongly disagree) to 7 (strongly agree).}
#' }
#' @source Simulated dataset (1000 observations) based on Chen, Y., Fang, Y., Yang, Y., & Dong, Y. (2026) Longitudinal associations among trust, loneliness, and life satisfaction among university students: A between- and within-person analysis.Journal of Personality, 94: 151-162. 
"Trust"


#' Smoking dataset
#'
#'
#' Simulated dataset for demonstrating GCLM() and STARTS() in longitudinal studies.
#'
#' @format A data frame with 500 rows and 13 variables:
#' \describe{
#'   \item{id}{Numeric, case id number.}
#'   \item{expose.1}{Numeric, Perceived exposure to smoking in movies at Time 1.}
#'   \item{expose.2}{Numeric, Perceived exposure to smoking in movies at Time 2.}
#'   \item{expose.3}{Numeric, Perceived exposure to smoking in movies at Time 3.}
#'   \item{expose.4}{Numeric, Perceived exposure to smoking in movies at Time 4.}
#'   \item{expose.5}{Numeric, Perceived exposure to smoking in movies at Time 5.}
#'   \item{expose.6}{Numeric, Perceived exposure to smoking in movies at Time 6.}
#'   \item{intens.1}{Numeric, Smoking intensity at Time 1.}
#'   \item{intens.2}{Numeric, Smoking intensity at Time 2.}
#'   \item{intens.3}{Numeric, Smoking intensity at Time 3.}
#'   \item{intens.4}{Numeric, Smoking intensity at Time 4.}
#'   \item{intens.5}{Numeric, Smoking intensity at Time 5.}
#'   \item{intens.6}{Numeric, Smoking intensity at Time 6.}
#' }
#' @source Dataset with 500 observations across six time points simulated with lavaan::simulateData() based on results of the Minnesota Adolescent Community Cohort (MACC) Study, 2000-2013, reported in Usami, S., Murayama, K., & Hamaker, E. L. (2019). A unified framework of longitudinal models to examine reciprocal relations. Psychological Methods, 24(5), 637-657.
"Smoking"


