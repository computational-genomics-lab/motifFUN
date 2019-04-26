#' Prediction of transmembrane helices in proteins
#'
#' This function Predicts the presence of transmembrane helices and extract and returns protiens with no transmembrane helices
#' @param tmhmm.path Localpath of folder containing tmhmm binary excutable file or the excutable file itself, if path not specified, then tmhmm must be in program search path
#' @param file Provide input fasta file name
#'
#' @export
#' @return Returns dataframe which contain summary of tmhmm result file
#'
#' @examples
#'  \dontrun{
#'
#' tmhmm_data <- get_tmhmm(tmhmm.path = NULL, file = NULL)
#' }
#'
get_tmhmm <- function(tmhmm.path = NULL, file = NULL) {
  tmhmm.path <- external_programe_existance(programe = "tmhmm", path = tmhmm.path)
  check_file_existance(input.file = file)
  output_filename <- modify_filename(file = input_file, programe_name = "tmhmm")
  tmhmm_inputfiles <- prepare_inputdata(file_input = input_file, dir_name = "tmhmmdir")
  tmhmm_run <- run_tmhmm(tmhmm.path = tmhmm.path, input_file.list = file, output = output_filename)
  tmhmm_summary_data <- tmhmm_summary(tmhmm_outfile = output_filename)
  return(tmhmm_summary_data)
}






