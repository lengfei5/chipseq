##################################################
## Project: General purpose 
## Script purpose: convert excel sheet (e.g. .csv, .xlsx) to a text file 
## Usage example : Rscript convertExcel2txt.R example.xlxs 
## Author: Jingkui Wang (jingkui.wang@imp.ac.at)
## Date of creation: Mon Dec  4 10:26:20 2017
##################################################
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0){
  stop("At least one argument must be supplied (input_file).xlsx", call.=FALSE)
}else{
  library(openxlsx)
  library(tools)
  for(n in 1:length(args))
  {
    filename = file_path_sans_ext(args[n])
    extension = file_ext(args[n])
    if(extension == "xlsx"){
      sheet = read.xlsx(args[n], sheet = 1, colNames = TRUE,  rowNames = FALSE)
    }else{
      if(extension == "csv"){
        sheet = read.csv(args[n], header = TRUE)
      }else{
        stop("Input file much .csv or .xlsx", call. = FALSE)
      }
    }
    
    newff = data.frame(sheet[ ,1], apply(sheet[, -1], 1, function(x) paste0(as.character(x), collapse = "_")))
    colnames(newff) = c("sampleID", "fileName")
    write.table(newff, file=paste0(filename, ".txt"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
    #print(sheet[1,])
    
    #if(extension == "txt"){
    #  sheet = read.xlsx(args[n], sheet = 1, colNames = TRUE,  rowNames = FALSE)
    #}
    
    #if (length(args)==1){
    # default output file
    #args[2] = "out.txt"
    #}
    
  }
}