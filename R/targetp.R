#' Predicts sub cellular location of eukaryotic proteins
#'
#' This function predicts the sub cellular location of eukaryotic proteins. The location assignment is based on the predicted presence of any of the N-terminal presequences:chloroplast transit peptide (cTP), mitochondrial targeting peptide (mTP) or secretory pathway signal peptide (SP).
#'
#' @param targetp.path Localpath of folder containing target binary excutable file or the excutable file itself, if path not specified, then targetp must be in program search path.
#' @param organismgroup Please provide oraganismgroup for targetp: -P|-N-use plant/non-plant networks mandatory, -c-include cleavage site predictions, -h-print this note and exit,  -v-print ve           rsion info and exit, -p-chloroplast prediction cutoff, default, -s-secretory pathway prediction cutoff, default, -t-mitochondrial prediction cutoff, default, -o-other location pred           iction cutoff, default.
#' @param file Please provide fasta input file
#' @export
#' @return Returns dataframe which contains summary of targetp result file
#'
#'  @examples
#'
#' get_targetp(targetp.path = NULL, organismgroup = NULL, file = input_file)
#'
get_targetp <- function(targetp.path = NULL, organismgroup = NULL, file = NULL) {
  targetp.path <- external_programe_existance(programe = "targetp", path = targetp.path)
  if (is.null(organismgroup)){
    stop(paste0("Please provide organismgroup name"))
  }
  check_file_existance(input.file = file)
  output_filename <- modify_filename(file = file, programe_name = "targetp")
  targetp_inputfiles <- prepare_inputdata(file_input = file, dir_name = "targetpdir")
  targetp_run <- run_targetp(programe_path = targetp.path, input_file.list = targetp_inputfiles, org.group = organismgroup, output = output_filename)
  targetp_summary_data <- targetp_summary(targetp_outfile = output_filename, organism = organismgroup)
  return(targetp_summary_data)
}



