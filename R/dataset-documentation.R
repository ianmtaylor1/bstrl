#' Simulated Noisy Records
#'
#' A dataset containing several files of noisy simulated records. Records are
#' simulated using GeCo (Tran, Vatsalan, and Cristen (2013)) and organized into
#' files of 200 records each. The columns in each file consist of two ID columns
#' for validating links:
#'
#' \itemize{
#'   \item rec.id. Contains the entity number and duplicate number of each record. This is unique to a record.
#'   \item entity. Contains the entity number of which this record is a copy. Is identical for all records which are noisy duplicates of the same original.
#' }
#'
#' The columns also consist of fields used to perform linkage, into which 3 errors
#' have been randomly inserted:
#'
#' \itemize{
#'   \item given.name, surname. Text fields with potential typographical errors.
#'   \item age, occup, extra1, ..., extra10. Categorical fields with potential swapped category errors.
#' }
#'
#' Linkage may be performed on either the full dataset or on only a subset of the fields.
#'
#' @format A list of 7 data.frames. Each data.frame has 200 rows and 16 columns.
#' @seealso geco_small
#' @usage data(geco_30over_3err)
"geco_30over_3err"

#' Simulated Noisy Records (smaller set)
#'
#' A dataset containing several files of noisy simulated records. Records are
#' simulated using GeCo (Tran, Vatsalan, and Cristen (2013)) and organized into
#' files of 10 records each. These files are subsets of the larger dataset.
#' The columns in each file consist of two ID columns for validating links:
#'
#' \itemize{
#'   \item rec.id. Contains the entity number and duplicate number of each record. This is unique to a record.
#'   \item entity. Contains the entity number of which this record is a copy. Is identical for all records which are noisy duplicates of the same original.
#' }
#'
#' The columns also consist of fields used to perform linkage, into which 3 errors
#' have been randomly inserted:
#'
#' \itemize{
#'   \item given.name, surname. Text fields with potential typographical errors.
#'   \item age, occup, extra1, ..., extra10. Categorical fields with potential swapped category errors.
#' }
#'
#' Linkage may be performed on either the full dataset or on only a subset of the fields.
#'
#' @format A list of 7 data.frames. Each data.frame has 10 rows and 16 columns.
#' @seealso geco_30over_3err
#' @usage data(geco_small)
"geco_small"
