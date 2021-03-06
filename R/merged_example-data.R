#' Simulated BP+PTSD+Polygenic Risk Score Data with covariates
#'
#' Simulated data for use in testing functions in polygenicPTSD package
#' Includes outcome variables SBP and DBP, along with: Polygenic Risk Scores,
#' PTSD variables, 3 Principle components, sex, education, cohort,
#' antihistamine use, bmi, and trauma count.
#'
#' @docType data
#'
#' @usage data(merged_example)
#'
#' @format A dataframe with 2264 rows and 30 columns
#'  \describe{
#'   \item{FID}{Family ID, for matching purposes}
#'   \item{IID}{Individual ID}
#'   \item{...prs...}{Polygenic Risk Score (PRS) variables per outcome (SBP, DBP, and PRS) and cutoff (01, 02)}
#'   \item{SBP_meas}{Simulated SBP measurements (OUTCOME)}
#'   \item{DBP_meas}{Simulated DBP measurement outcome}
#'   \item{ptsd...}{Four PTSD measurements}
#'   \item{P1, P2, P3}{First 3 Principle Components}
#'   \item{sex}{Sex (0 = Female, 1 = Male)}
#'   \item{age}{Numeric, years}
#'   \item{educ}{...}
#'   ...
#' }
#'
#' @keywords datasets
#'
#'
#' @examples
#' data(merged_example)
#' classed <- htn_outcome_classify(merged_example)
#'
"merged_example"
