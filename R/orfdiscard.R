#' Descard ORF file
#'
#'  This function descard the ORF protien sequences based on its length. The length cutoff is managed by lowerlimit and upperlimit.
#'
#' @param orf_file ORF file created from getorf function.
#' @param upper.limit enter sequnece length it will descard the sequences within specified range.
#' @param lower.limit enter sequence length.
#' @param output_file output filename defined by user.
#'
#' @export
#' @return This function returns output file name
#'
#'  @examples
#'
#' orf_discard(orf_file= "filename", upper.limit= 1800, lower.limit= 100,  output_file = NULL)
#'
#'
orf_discard <- function(orf_file = NULL, upper.limit = NULL, lower.limit = NULL, output_file = NULL) {
  check_file_existance(input.file = orf_file)
  if (is.null(upper.limit)){
    stop(paste0("Please enter upper limit"), call. = FALSE)
  }
  if (is.null(lower.limit)){
    stop(paste0("Please enter lower limit"), call. = FALSE)
  }
  sequence <- list()
  ID <- list()
  word_count <- list()
  extract_seq <- list()
  read_seqence <- seqinr::read.fasta(orf_file)
  cat("discarding...\n")
  for (i in 1:length(read_seqence)) {
    extract_seq[i] <- seqinr::getSequence(read_seqence[i], as.string = TRUE)
    word_count[i] <- sapply(gregexpr("\\w",extract_seq[i]), length)
  }
  for (j in 1:length(word_count)) {
    if ((word_count[j] > lower.limit) & (word_count[j] < upper.limit)){
      sequence[j] <- read_seqence[j]
      ID[j] <- seqinr::getName(read_seqence[j])
    }
  }
  ID_null <- which(!sapply(ID, is.null))
  ID <- ID[ID_null]
  names(ID) <- ID_null
  Seq_null <- which(!sapply(sequence, is.null))
  sequence <- sequence[Seq_null]
  names(sequence) <- Seq_null
  if (!is.null(output_file)){
    out_file <- output_file
  } else {
    output_filename <- modify_filename(file = orf_file, programe_name = "descard.fasta", concat_string = "_")
    out_file <- output_filename
  }
  seqinr::write.fasta(sequence,ID, file.out = out_file)
  cat("completed.\n")
  filepath <- paste0(getwd(), "/", out_file)
  return(filepath)
}




