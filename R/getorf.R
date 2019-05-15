#' Find and extract open reading frames (ORFs)  from nucleotide sequences
#'
#' This function finds and outputs the sequences of open reading frames (ORFs) in one or more nucleotide sequences
#'
#' @param getorf.path  Path of getorf excutable file path,path.getorf excutable file will provide along with the package,if path is null then this function will take default path from syste path.
#' @param input.file  inputfile should have DNA query sequence filename.
#' @param output.file name of output filename.
#' @param error If \code{TRUE}, throw an error if programes is not installed.
#' @param verbose If \code{TRUE}, print progress reports.
#'
#' @export
#' @return return output filename
#'
#'  @examples
#'
#' get_orf( getorf.path=NULL, input.file="filename", output.file = NULL)
#'
get_orf <- function(getorf.path=NULL, input.file = NULL, output.file = NULL, error = TRUE, verbose = FALSE) {
  getorf.path <- external_programe_existance(programe = "getorf", path = getorf.path)
  check_file_existance(input.file = input.file)
  if (!is.null(output.file)){
    out_file <- output.file
  } else {
    output_filename <- modify_filename(file = input.file, programe_name = "orf.fasta", concat_string = "_")
    out_file <- output_filename
  }
  orfcommand <- c(getorf.path, input.file, outputfile = out_file)
  system2(orfcommand)
  return(out_file)
}



