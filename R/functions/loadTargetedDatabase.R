
#' Load targeted database
#'
#' Loads a .xlsx or .csv database files that contains all the information needed to detect targeted compounds.
#'
#' @param path_to_file Full path to the file containing the database with targeted compounds. Polarity op ions needs to be indicated in the format: "[M+ION]+" or "[M-ION]-"  
#' @param polarity Indicates the polarity, can be either "negative" or "positive"
#' @param columnNames Column names of the columns with information that need to be included, can be as many as needed, but needs to include "ID","Name","m/z-value" and "RT (min)"
#'
#' @return A dataframe containing the targeted compounds
#' @export
#'
#' @examples

loadTargetedDatabase <- function(path_to_file,polarity,columnNames){

  if(grepl(".xlsx",path_to_file) == T){
    masslist <- read_excel(path_to_file)  
  }
  else if(grepl(".csv",path_to_file) == T ){
    masslist <- read.csv(path_to_file)
  }
  else {
    stop("The database input format should be .csv, .txt or .xlsx")
  }
  
  masslist_negative <- masslist[grep("]-",masslist$`Ion adduct`,fixed = T),]
  masslist_positive <- masslist[grep("]+",masslist$`Ion adduct`,fixed = T),]
  
  #Remove unnecessary columns and rename
  
  masslist_positive <- masslist_positive[,c("ID","Name","m/z-value","RT (min)")]
  masslist_negative <- masslist_negative[,c("ID","Name","m/z-value","RT (min)")]
  colnames(masslist_positive) <- c("ID","NAME","m/z","tr")
  colnames(masslist_negative) <- c("ID","NAME","m/z","tr")
  
  #Set RT in seconds & make numeric
  
  masslist_positive$tr <- as.numeric(masslist_positive$tr) *60
  masslist_positive$`m/z` <- as.numeric(masslist_positive$'m/z')
  
  masslist_negative$tr <- as.numeric(masslist_negative$tr) *60
  masslist_negative$`m/z` <- as.numeric(masslist_negative$'m/z')
  
  #set ID as character
  
  masslist_positive$ID <- as.character(masslist_positive$ID)
  masslist_negative$ID <- as.character(masslist_negative$ID)
  
  if(polarity == "positive"){
    dbData <- masslist_positive
  }
  if (polarity == "negative"){
    dbData <- masslist_negative
  }
  
  return(dbData)
  
}
