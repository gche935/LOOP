#' Trust dataset
#'
#'
#' Simulated dataset for demonstrating LOOP functions.
#'
#' @format A data frame with 1000 rows and 19 variables:
#' \describe{
#'   \item{ID}{Numeric, case ID number.}
#'   \item{trust_T1}{Numeric, Average of 6 Trust items at Time 1: from 1 (strongly disagree) to 5 (strongly agree).}
#'   \item{trust_T2}{Numeric, Average of 6 Trust items at Time 2: from 1 (strongly disagree) to 5 (strongly agree).}
#'   \item{trust_T3}{Numeric, Average of 6 Trust items at Time 3: from 1 (strongly disagree) to 5 (strongly agree).}
#'   \item{trust_T4}{Numeric, Average of 6 Trust items at Time 4: from 1 (strongly disagree) to 5 (strongly agree).}
#'   \item{trust_T5}{Numeric, Average of 6 Trust items at Time 5: from 1 (strongly disagree) to 5 (strongly agree).}
#'   \item{trust_T6}{Numeric, Average of 6 Trust items at Time 6: from 1 (strongly disagree) to 5 (strongly agree).}
#'   \item{lonely_T1}{Numeric, Average of 6 UCLS Loneliness Scale at Time 1: from 1 (never) to 4 (always).}
#'   \item{lonely_T2}{Numeric, Average of 6 UCLS Loneliness Scale at Time 2: from 1 (never) to 4 (always).}
#'   \item{lonely_T3}{Numeric, Average of 6 UCLS Loneliness Scale at Time 3: from 1 (never) to 4 (always).}
#'   \item{lonely_T4}{Numeric, Average of 6 UCLS Loneliness Scale at Time 4: from 1 (never) to 4 (always).}
#'   \item{lonely_T5}{Numeric, Average of 6 UCLS Loneliness Scale at Time 5: from 1 (never) to 4 (always).}
#'   \item{lonely_T6}{Numeric, Average of 6 UCLS Loneliness Scale at Time 6: from 1 (never) to 4 (always).}
#'   \item{lifesat_T1}{Numeric, Average of 5 Satisfaction with Life Scale at Time 1: from 1 (strongly disagree) to 7 (strongly agree).}
#'   \item{lifesat_T2}{Numeric, Average of 5 Satisfaction with Life Scale at Time 2: from 1 (strongly disagree) to 7 (strongly agree).}
#'   \item{lifesat_T3}{Numeric, Average of 5 Satisfaction with Life Scale at Time 3: from 1 (strongly disagree) to 7 (strongly agree).}
#'   \item{lifesat_T4}{Numeric, Average of 5 Satisfaction with Life Scale at Time 4: from 1 (strongly disagree) to 7 (strongly agree).}
#'   \item{lifesat_T5}{Numeric, Average of 5 Satisfaction with Life Scale at Time 5: from 1 (strongly disagree) to 7 (strongly agree).}
#'   \item{lifesat_T6}{Numeric, Average of 5 Satisfaction with Life Scale at Time 6: from 1 (strongly disagree) to 7 (strongly agree).}
#' }
#' @source Simulated dataset (500 observations) based on Chen, Y., Fang, Y., Yang, Y., & Dong, Y. (2026) Longitudinal associations among trust, loneliness, and life satisfaction among university students: A between- and within-person analysis.Journal of Personality, 94: 151-162. 
"Trust"


#' Smoking dataset
#'
#'
#' Simulated dataset for demonstrating GCLM() and STARTS() in longitudinal studies.
#'
#' @format A data frame with 500 rows and 13 variables:
#' \describe{
#'   \item{id}{Numeric, case id number.}
#'   \item{expose_T1}{Numeric, Perceived exposure to smoking in movies at Time 1.}
#'   \item{expose_T2}{Numeric, Perceived exposure to smoking in movies at Time 2.}
#'   \item{expose_T3}{Numeric, Perceived exposure to smoking in movies at Time 3.}
#'   \item{expose_T4}{Numeric, Perceived exposure to smoking in movies at Time 4.}
#'   \item{expose_T5}{Numeric, Perceived exposure to smoking in movies at Time 5.}
#'   \item{expose_T6}{Numeric, Perceived exposure to smoking in movies at Time 6.}
#'   \item{intens_T1}{Numeric, Smoking intensity at Time 1.}
#'   \item{intens_T2}{Numeric, Smoking intensity at Time 2.}
#'   \item{intens_T3}{Numeric, Smoking intensity at Time 3.}
#'   \item{intens_T4}{Numeric, Smoking intensity at Time 4.}
#'   \item{intens_T5}{Numeric, Smoking intensity at Time 5.}
#'   \item{intens_T6}{Numeric, Smoking intensity at Time 6.}
#' }
#' @source Dataset with 500 observations across six time points simulated with lavaan::simulateData() based on results of the Minnesota Adolescent Community Cohort (MACC) Study, 2000-2013, reported in Usami, S., Murayama, K., & Hamaker, E. L. (2019). A unified framework of longitudinal models to examine reciprocal relations. Psychological Methods, 24(5), 637-657.
"Smoking"


#' Responsive dataset
#'
#'
#' Simulated daily diary dataset for demonstrating ML function.
#'
#' @format A data frame with 1620 rows and 13 variables:
#' \describe{
#'   \item{id}{Numeric, case id number.}
#'   \item{gender}{Numeric, Employee gender: 0 = female and 1 = male.}
#'   \item{age}{Numeric, Employee age.}
#'   \item{edu}{Numeric, Employee education level.}
#'   \item{tenure}{Numeric, Employee organizational tenure.}
#'   \item{proactive}{Numeric, Average of 4 proactive personality items: from 1 (strongly disagree) to 7 (strongly agree).}
#'   \item{time}{Numeric, Day of data collection.}
#'   \item{act_con}{Numeric, Average of 3 active-constructive items: from 1 (strongly disagree) to 7 (strongly agree).}
#'   \item{pas_con}{Numeric, Average of 3 passive-constructive items: from 1 (strongly disagree) to 7 (strongly agree).}
#'   \item{act_des}{Numeric, Average of 3 active-destructive responsiveness items: from 1 (strongly disagree) to 7 (strongly agree).}
#'   \item{pas_des}{Numeric, Average of 3 passive-destructive responsiveness items: from 1 (strongly disagree) to 7 (strongly agree).}
#'   \item{engage}{Numeric, Average of 9 work engagement items: from 1 (strongly disagree) to 7 (strongly agree).}
#'   \item{satisfaction}{Numeric, Average of 6 job satisfaction items: from 1 (strongly disagree) to 7 (strongly agree).}
#' }
#' @source Simulated dataset (180 observations) based on Liu, X., Dong, J., Yu, Y., Zheng, J., Zheng, X., & Lee, B. Y. (2026). From friction to fulfillment: Examining when and how spousal active-destructive responsiveness to employees' sharing of positive events benefits the work domain. Journal of Organizational Behavior, 47,  66-80. 
"Responsive"

