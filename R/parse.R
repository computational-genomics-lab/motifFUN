#' Parses single fasta file into multiple files
#'
#' This function Parses the large fasta file and store the sequences into multiple files according to user requirement
#'
#' @param path_in A comlpete path to the input fasta sequence file
#' @param path_out A Path where resulting files will be store.Path should be contain prefix name of input file on which will append .fa extenstion
#' @param num_seq Integer defining howmany sequences want to seperate in to one file
#' @param trim Logical, should the sequences be trimmed to 4000 amino acids to bypass the CBS server restrictions. Defaults to FALSE.
#' @param trunc Integer, truncate the sequences to this length.
#' @param id Logical, should the protein id's be returned. Defaults to FALSE.
#'
#' @export
#' @return if id = FALSE, A Character vector of the paths to the resulting .FASTA formatted files.
#'
#'
#'  @examples
#' parse_file(path_in= "path of input file", path_out= "path to store outputfiles", num_seq= 300,trim = FALSE, trunc = 0, id = TRUE)


parse_file <- function(path_in =  NULL, path_out = NULL, num_seq = NULL, trim = FALSE, trunc = 0, id = TRUE){
  if (!file.exists(path_in)) {
    stop("cannot find file in the specified path_in")
  }
  if (missing(num_seq)) {
    num_seq <- 20000
  }
  if (length(num_seq) > 1) {
    num_seq <- 20000
    warning("num_seq should be of length 1, setting to default: num_seq = 20000")
  }
  if (!is.numeric(num_seq)) {
    num_seq <- as.numeric(num_seq)
    warning("num_seq is not numeric, converting using 'as.numeric'")
  }
  if (is.na(num_seq)) {
    num_seq <- 20000
    warning("num_seq was set to NA, setting to default: num_seq = 20000")
  }
  if (is.numeric(num_seq)) {
    num_seq <- floor(num_seq)
  }
  if (!missing(trunc)) {
    if (length(trunc) > 1) {
      stop("trunc should be of length 1.")
    }
    if (!is.numeric(trunc)) {
      stop("trunc is not numeric.")
    }
    if (is.na(trunc)) {
      stop("trunc was set to NA.")
    }
    if (is.numeric(trunc)) {
      trunc <- floor(trunc)
    }
    if (trunc < 0) {
      stop("trunc was set to a negative number.")
    }
    if (trunc == 0) {
      trunc <- 1000000L
    }
  }
  if (missing(trim)) {
    trim <- FALSE
  }

  if (length(trim) > 1) {
    trim <- FALSE
  }
  if (!is.logical(trim)) {
    trim <- as.logical(trim)
    warning("trim is not logical, converting using 'as.logical'")
  }
  if (is.na(trim)) {
    trim <- FALSE
    warning("trim was set to NA, setting to default: trim = FALSE")
  }
  if (missing(id)) {
    id <- FALSE
  }
  if (length(id) > 1) {
    id <- FALSE
    warning("id should be of length 1, setting to default: id = FALSE")
  }
  if (!is.logical(id)) {
    id <- as.logical(id)
    warning("id is not logical, converting using 'as.logical'")
  }
  if (is.na(id)) {
    id <- FALSE
    warning("id was set to NA, setting to default: id = FALSE")
  }
  temp_file <- seqinr::read.fasta(file = path_in, seqtype = "AA")
  name <- names(temp_file)
  if (!missing(trunc)) {
    temp_file <- lapply(temp_file, function(x) {
      len <- length(x)
      if (len > trunc) {
        out <- x[1:trunc]

      }
      else {
        out <- x
      }
      return(out)
    })
  }
  if (trim == TRUE && missing(trunc)) {
    temp_file <- lapply(temp_file, function(x) {
      len <- length(x)
      if (len > 4000) {
        out <- x[1:4000]
      }
      else {
        out <- x
      }
      return(out)
    })
  }
  len <- length(temp_file)
  print(len)
  splt <- num_seq
  print(splt)
  pam <- ((seq(len) - 1)%/%splt) + 1
  print(pam)
  m_split <- split(temp_file, pam)
  #print(m_split)
  file_list <- vector("character", length(m_split))
  for (i in 1:length(m_split)) {
    seqinr::write.fasta(sequences = m_split[[i]], names = names(m_split[[i]]),
                        file.out = paste(path_out, i, ".fa", sep = ""))
    file_list[i] <- paste(path_out, i, ".fa", sep = "")
  }
  if(id){
    return(list(id = name, file_list = file_list))
  } else {
    return(file_list)
  }

}




