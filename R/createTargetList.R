#' Create target list
#'
#' Loads a .xlsx or .csv file that contains all the information needed to detect targeted compounds.
#'
#' @param input_directory_targets Full path to the file containing the database with targeted compounds.
#' @param pos_pattern Pattern of positive ions to be recognized in  ion_column
#' @param neg_pattern Pattern of negative ions to be recognized in  ion_column
#' @param ion_column Name of the column in which the type op ionisation is found
#' @param polarity Indicates the polarity, can be either "negative" or "positive"
#' @param columns_of_interest Column names of the columns with information that need to be included, can be as many as needed, but needs to include "ID","Name","m/z-value" and "RT (min)"
#'
#' @return A dataframe containing the targeted compounds
#' @export
#'
#' @examples

createTargetList <- function(input_directory_targets,pos_pattern,neg_pattern,polarity,ion_column,columns_of_interest){

  if (str_ends(input_directory_targets,"csv") == TRUE) {
    masslist <-
      read.csv(input_directory_targets,
               header = TRUE,
               sep = ",",check.names = TRUE)
  } else if (str_ends(input_directory_targets,"xlsx") == TRUE) {
    masslist <- read_xlsx(input_directory_targets)
  }

  masslist <- data.frame(masslist)

  masslist_negative <-
    masslist[grep(neg_pattern, masslist[,ion_column], fixed = T), ]
  masslist_positive <-
    masslist[grep(pos_pattern, masslist[,ion_column], fixed = T), ]

  #Remove unnecessary columns and rename

  masslist_positive <-
    masslist_positive[, columns_of_interest]
  masslist_negative <-
    masslist_negative[, columns_of_interest]
  colnames(masslist_positive) <- c("ID", "NAME", "m/z", "tr")
  colnames(masslist_negative) <- c("ID", "NAME", "m/z", "tr")

  #Set RT in seconds & make numeric

  masslist_positive$tr <- as.numeric(masslist_positive$tr) * 60
  masslist_positive$`m/z` <- as.numeric(masslist_positive$'m/z')

  masslist_negative$tr <- as.numeric(masslist_negative$tr) * 60
  masslist_negative$`m/z` <- as.numeric(masslist_negative$'m/z')

  #set ID as character

  masslist_positive$ID <- as.character(masslist_positive$ID)
  masslist_negative$ID <- as.character(masslist_negative$ID)

  if (polarity == "positive") {
    dbData <- masslist_positive
  }
  if (polarity == "negative") {
    dbData <- masslist_negative
  }

  dbData <- na.omit(dbData)

  return(dbData)

   }
