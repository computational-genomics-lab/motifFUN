#' Returns non-redundant sequences from hmm.search or regex.search and generates a motif table
#'
#' This function summarize the results from pattern.search or hmm.search along with functional domains.
#' @param hmm.result A list of SeqFastadna objects obtained from regex.search or hmm.search
#' @param reg.pat A character string indicating the PATT_REG pattern for the motif. The specification of the PATT_REG pattern in must be in regex format.
#' @param signalp_version signalp version
#' @param input_file file name used in the functional domains(signalp, targtp, tmhmm)
#' @export
#' @return returns summary table
#'
#' @examples
#' \dontrun{
#'
#' summary_motifs(hmm.result, reg.pat=NULL, signalp_version = NULL, input_file = NULL)
#'}
#'
summary_motifs <- function (hmm.result, reg.pat=NULL, signalp_version = NULL, input_file = NULL){
  output_tmhmm <- modify_filename(file = input_file, programe_name = "tmhmm")
  tmhmm_summary_data <- tmhmm_summary(tmhmm_outfile = output_tmhmm)
  output_targetp  <- modify_filename(file = input_file, programe_name = "targetp")
  targetp_summary_data <- targetp_summary(targetp_outfile = output_targetp, organism = "-P")
  if (signalp_version == 3){
    output_signalp3 <- modify_filename(file = input_file, programe_name = "signalp3")
    signalp_summary_data <- summary_signalp3(signalp3_outname = output_signalp3, original_input = input_file)
  } else {
    output <- modify_filename(file = input_file, programe_name = "summary.signalp5", concat_string = "_")
    siganlp5_out_file <- readLines(output)
    signalp_summary_data <- utils::read.table(text = siganlp5_out_file, sep = "\t")
    signalp_summary_data$signal_Peptide <- signalp_summary_data$V2
  }
  motif.out <- list()
  if (length(hmm.result) == 4) {
    hmm.result <- hmm.result[names(hmm.result) %in% c("REGEX","HMM","HMM_Table")]
  }
  if (length(hmm.result) == 3) {
    summ.dataf <- c(hmm.result[[1]],hmm.result[[2]])
    consensus.seq <- summ.dataf[!duplicated(seqinr::getName(summ.dataf))]
    motif.out$consensus.sequences <- consensus.seq
  } else {
    consensus.seq <- hmm.result
  }
  sequences <- lapply(consensus.seq, function (x) paste(unlist(x),collapse = ""))
  if (is.null(reg.pat)){
    stop("No custom REGEX pattern found.\n The 'custom' option requires a mandatory REGEX pattern")
  } else {
    custom.num <- as.numeric(gsub(pattern = "", replacement = "", unlist(lapply(lapply(sequences, function (x) unlist(gregexpr(x, pattern=reg.pat, perl = T,ignore.case = T))), function (x) length(x)))))
    custom.motif <- unlist(lapply(lapply(sequences, function (x) unlist(gregexpr(x, pattern=reg.pat, perl = T,ignore.case = T))), function (x) paste(x,collapse = ",")))
    motifs <- data.frame(seqinr::getName(consensus.seq),custom.num,custom.motif, stringsAsFactors = F)
    rownames(motifs) <- NULL
  }
  motifs[,2][is.na(motifs[,3])] <- 0
  colnames(motifs) <- c("Sequence ID","Motif number","Motifposition")
  if (nrow(tmhmm_summary_data) != nrow(motifs)){
    motifs <- motifs[!motifs$`Motifposition` == -1, ]
  }
  motifs <- data.frame(motifs, length = unlist(lapply(motifs[,1], function (x) length(unlist(consensus.seq[seqinr::getName(consensus.seq) %in% x])))))
  motif.out$motif.table <- motifs
  functional_data <- data.frame(matrix(nrow = nrow(motif.out$motif.table), ncol = 4))
  colnames(functional_data) <- c("targetp_localization", "NUMBER_OF_PREDICTED_TMHS", "signalpeptide", "signalp_HMM_probability")
  functional_data$targetp_localization <- targetp_summary_data$LOCALIZATION
  functional_data$NUMBER_OF_PREDICTED_TMHS <- tmhmm_summary_data$NUMBER_OF_PREDICTED_TMHS
  functional_data$signalpeptide <- signalp_summary_data$signal_Peptide
  functional_data$signalp_HMM_probability <- signalp_summary_data$HMM_probability
  motif.out <- cbind( motif.out$motif.table, functional_data)
  return(motif.out)
}
