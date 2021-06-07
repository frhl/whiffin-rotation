#' @title read SaintExpress workbooks 
#' @description reads in a workbook using the readxl tidyverse library.
#' @param workbook string. Path to workbook.
#' @param sheets string(s). What sheets should be read? Default is reading all.
#' @param rename_sheets strings(s). Should sheetes be renamed?
#' @param readxl_function what function should be used to read in objects. Takes
#' the first two arguments as inputs.
#' @export

read_workbook <- function(workbook, sheets = NULL, rename_sheets = NULL, readxl_function = read_xlsx){
  
  if (is.null(sheets)) sheets = excel_sheets(workbook)
  
  # get sheet mapping
  sheet_mapping <- ifelse(is.null(rename_sheets), list(as.list(sheets)), list(as.list(rename_sheets)))[[1]]
  names(sheet_mapping) <- sheets
  
  # iterate over every sheet
  result <- lapply(sheets, function(sheet) {
    in_data = readxl_function(workbook, sheet)
    in_data$Sheet <- sheet_mapping[[sheet]]
    return(in_data)
  })
  
  names(result) <- unlist(ifelse(is.null(rename_sheets), list(sheets), list(rename_sheets)))
  return(result)
}


