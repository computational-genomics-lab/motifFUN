#' Searching for motifs using regular expressions
#'
#' This function searching the motif of interest by user defined regular expression against input fasta sequences.
#'
#' @param fasta.file Provide protien sequence fasta file.
#' @param reg.pat Provide R regular expression pattern for searching motifs.
#'
#' @export
#'
#' @return A list contain motif sequences and motif dataframe
#'
#'  @examples
#'
#'  pattern.search(fasta.file = NULL, reg.pat = NULL)
#'
pattern.search <- function(fasta.file = NULL, reg.pat = NULL) {
  ORF <- seqinr::read.fasta(fasta.file)
  seq <- lapply(ORF, function (x) paste(unlist(x),collapse = ""))
  if (missing(fasta.file)) {
  stop('Missing argument: fasta.file', call. = FALSE)
  }
  if (is.null(reg.pat)) {
    stop('plese provide regular expression', call. = FALSE)
  }
  if (missing(fasta.file)) {
    stop('please provide input file', call. = FALSE)
  }
  if (is.null(reg.pat)) {
    stop('plese provide regular expression', call. = FALSE)
  }
  regex <- list()
  j=0
  for (i in 1:length(seq)){
    regex[[i]] <- unlist(gregexpr(seq[[i]], pattern = reg.pat, perl = T ,ignore.case = T))
    if (length(regex[[i]])>0) {
      var <- length(regex[[i]])
      sort(var)
     for (i in var) {
        if(i>j){
         j=i
        }
      }
    }
  }
  col_names <- list()
  for (i in 1:j) {
    col_names[i]<- paste("Pos")
  }
  col_no <- j+1
  regex1 <-data.frame(matrix( nrow = length(ORF), ncol = (j+1)))
  colnames(regex1) <- c(col_names, "Motif_number")
  for (i in 1:length(ORF)) {
    if (length(unlist(regex[i])) == 1){
      regex1[i,1] <- unlist(regex[i])
    } else {
      regex1[i,1] <- unlist(regex[i])[1]
      regex1[i,2] <- unlist(regex[i])[2]
    }
    regex1[i,col_no]<- length(unlist(regex[i]))
  }
  regex1$SEQUENCE_ID <- names(seq)
  regex1 <- regex1[!regex1$Pos < 0, ]
  motif_data <- data.frame(regex1)
  regex1 <-ORF[seqinr::getName(ORF) %in% regex1$SEQUENCE_ID]
  if (length(regex1) == 0) {
    stop(paste0("No motif sequences are found."))
  }
  output_filename <- modify_filename(file = fasta.file, programe_name = "REGEX.fasta", concat_string = "_")
  seqinr::write.fasta(sequences =regex1, names = names(regex1), file.out = output_filename)
  print(length(regex1))
  return(regex1)
}


