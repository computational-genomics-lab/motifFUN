#' Predicts the secretory proteins
#'
#' This function predicts the presence and location of signal peptide cleavage sites in amino acid sequences from different organisms(eukaryotes,gram+,gram-).
#'
#' @param signalp.path Localpath of folder containing signalp binary excutable file or the the the excutable file itself, if path not specified, then signalp must be in program search path
#' @param signalp.version please provide signalp version to use
#' @param input_file Please provide fasta input file
#' @param org.type Provide organism type [euk,gram+,gram-]
#'
#' @export
#' @return Provide  dataframe its contain summary of the signalp results file
#'
#'  @examples
#'
#'  signalp <- get_signalp(signalp.path = NULL, signalp.version = NULL, input_file = NULL, org.type = NULL)
#'
get_signalp <- function(signalp.path = NULL, signalp.version = NULL, input_file = NULL, org.type = NULL){
  signalp.path <- external_programe_existance(programe = "signalp", path = signalp.path)
  check_file_existance(input.file = input_file)
  if (missing(signalp.version)) {
    stop('Missing argument: version.', call. = FALSE)
  }
  if (is.null(org.type)) {
    stop(paste0("Please provide organism type"), call. = FALSE)
  }
  program_version <- check_signalp_version(signalp.path)
  check_signalp_version_match(program_version, signalp.version)
  output_filename <- modify_filename(file = input_file, programe_name = "signalp3")
  inputfile_list <- prepare_inputdata_signalp(signalp.version, input_file)
  signalp <- run_signalp(programe_path = signalp.path, version = signalp.version, input_file.list = inputfile_list, org.type = org.type, output = output_filename)
  if (signalp.version == 3){
    signalp_summary_data <- summary_signalp3(signalp3_outname = output_filename, original_input = input_file)
  } else if (signalp.version == 5){
    signalp_summary_data <- siganlp5_summary(signalp5out = signalp)
  }
  return(signalp_summary_data)

}





