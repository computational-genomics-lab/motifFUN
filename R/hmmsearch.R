#' Searching for motifs using HMM searches
#'
#' This function uses MAFFT and HMMER to search for sequences with protien motifs using hidden markov models.
#'
#' @param original.seq The absolute path for the original six-frame translation FASTA file
#' @param regex.seq A list[1] of objects returns from pattern.search. The HMM profile will be constructed using these sequences
#' @param alignment.file The absolute path for an alignment file of the sequences to build the hmmer profile from.It's recommended that the alignment file cointains the same sequences than the `regex.seq` files. If the user provides the absolute path, *motiFUN* won't use MAFFT to align the sequences and will use the alignment for the HMMER searches. If no alignment file is provided, *motifFUN* will use MAFFT to align the sequences from `regex.seq` and run HMMER.
#' @param save.alignment Save the alignment in the returning object. The MAFFT alignment will be saved as the first element of the returned object
#' @param mafft.path Local path of folder containing the MAFFT binary executable file or the executable file itself. If not specified, then MAFFT must be in the program search path.
#' @param num.threads Number of threads to be used by MAFFT.
#' @param hmm.tresh  Set the bit score cutoff for the per-sequence ranked hit list to a real number. (Default = 0)
#' @param hmm.path Local path of  folder containing the HMMER binaries.  If not specified, then HMMER executables must be in the program search path.
#' @param seed The seed to used with HMMER commands. Set this to get the same output each time.
#'
#' @export
#' @return Returns total sequence.
#'
#'  @examples
#'
#' hmm.search(original.seq = fasta.file, regex.seq = PATT_REG, mafft.path = "complete path", hmm.path = "complete path")
#'
hmm.search <-  function(original.seq, regex.seq, alignment.file = NULL, save.alignment = FALSE, mafft.path = NULL, num.threads = 2,  hmm.tresh = 0, hmm.path = NULL, seed = 12345){
  set.seed(seed)
  sequences <- regex.seq
  if (unique(unlist(lapply(sequences, class))) != "SeqFastadna") {
    stop("The object is not a list of sequences read by seqinr.")
  }
  if (length(regex.seq) < 4) {
    stop("Not enough sequences for HMM step. At least 4 sequences are required.")
  }
  if (is.null(alignment.file) == F & save.alignment == T ) {
    stop("Alignment is provided. No new alignment to save")
  }
  if (class(hmm.tresh) != "numeric") {
    stop("hmm.thres requires a real number to set the bit-score threshold.")
  }
  #mafft alignment
  mafft.out.name <- get_mafft(alignment.file = alignment.file, sequences = sequences, mafft.path = mafft.path, num.threads = num.threads)
  cat ("MAFFT Alignment file is ", mafft.out.name, "\n")
  if (file.exists(mafft.out.name) == F) {
     stop("No MAFFT alignment found")
  }
  #HMM for windows
  if (Sys.info()[['sysname']] %in% "Windows") {
    #HMM build
    hmmbuild <- hmm_build(hmm_path = hmm.path, seed = seed)
    #HMM Press
    hmmpress <- hmm_press(hmm.path)
    #HMM search
    hmmsearch <- hmmsearch(hmmsearch_path = hmm.path, seed = seed,hmm.tresh = hmm.tresh, original.seq = original.seq)
  } else{
    #HMM build
    hmmbuild <- hmm_build(hmm_path = hmm.path, seed = seed, mafft_out = mafft.out.name)
    #HMM Press
    hmm_press(path_hmm = hmm.path)
    #HMM search
    hmmsearch <- hmm_search (hmmsearch_path = hmm.path, seed = seed, hmm.tresh = hmm.tresh,  original.seq = original.seq)
  }
  #Reading HMM results and producing motif candidate
  total.seq <- reading_hmm.results(original_seq = original.seq, hmmsearch.out = hmmsearch, hmmbuild.out = hmmbuild, regex.seq = regex.seq)
  #Generating HMM sequence file
  HMM_out.file <- HMM_fasta(total_seq = total.seq, original.seq = original.seq)

  return(total.seq)
}


