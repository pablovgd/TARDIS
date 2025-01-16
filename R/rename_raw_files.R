#' Rename Thermo .raw files
#'Script to read XCalibur sequence list and rename files according to run type
#'(sample, QC, blank, etc...)
#' @param data_path_raw_files path to raw files. Important: all runs from the
#'     sequence list need to have a corresponding .raw file in your input folder.
#' @param data_path_list Path to exported sequence list in .csv format,
#'     all columns may be exported
#'
#' @importFrom stringr str_replace_all
#'
#' @export
#'
#' @author Pablo Vangeenderhuysen
renameRawFiles <- function(data_path_raw_files,data_path_list){
file_names_original <- list.files(data_path_raw_files)
file_names <- gsub(".raw","",list.files(data_path_raw_files))
#gather original file names and remove .raw
sequence <- read.csv(data_path_list,skip = 1,check.names = FALSE) #read sequence list
sequence[,"save_name"] = NA
for(k in 1:length(file_names)){
  if(file_names[k] == sequence$`File Name`[k]){
    sequence$save_name[k] <- paste0(file_names[k],"_",
                                    sequence$`Sample Name`[k],".raw")
  }
}
sequence$save_name <- str_replace_all(sequence$save_name, "/", "")
new_names <- sequence$save_name
file.rename(paste0(data_path_raw_files, file_names_original),
            paste0(data_path_raw_files, new_names))
}
