#' For checking programe existance and set default path
#'
#' This function check the program and file paths.
#'
#' @param programe parameter for programe name.
#' @param path for path of the programe.
#' @param error If \code{TRUE}, throw an error if programs is not installed.
#' @param verbose If \code{TRUE}, print progress reports.
#'
#' @export
#' @return provide path of the programe
#'
#'
#'  @examples
#'
#' external_programe_existance (programe = NULL, path = NULL, error = TRUE, verbose = FALSE)
#'

### checking external programes existance

external_programe_existance <- function(programe = NULL, path = NULL, error = TRUE, verbose = FALSE){
  if (is.null(path)) {
    path <- unname(Sys.which(programe))
  } else {
    path <- file.path(path)
  }

### Check if TMHMM is installed #########

  is_installed <- file.exists(path) == T
  if (! is_installed && error) {
    if (is.null(path)) {
      stop(paste0( programe, "path not found in your computer's search path.",
                "'\n Please check that", programe , "is installed and in the search path or specify the path to the", programe, "installation using the", path, "option."), call. = FALSE)
   } else {
      stop(paste0(programe, "path not found in the specified path: '", path,
                "'\n Please check your", programe, "installation."), call. = FALSE)
    }
  }
  return(path)
}


### check input file existance

#' For checking file existance
#'
#' This function check the file paths.
#'
#' @param input.file input sequence file
#' @param error If \code{TRUE}, throw an error if programs is not installed.
#' @param verbose If \code{TRUE}, print progress reports.
#'
#' @export
#' @return provide logical
#'
#'
#'  @examples
#'  check_file_existance(input.file = NULL, error = TRUE, verbose = FALSE)


### check input file existance
check_file_existance <- function(input.file = NULL, error = TRUE, verbose = FALSE) {
  if (is.null(input.file)){
    stop(paste0("Please provide input file"), call. = FALSE)
  }

  if (missing(input.file)){
    stop('input.file is missing', call. = FALSE)
  }
  read_file <- readLines(input.file)
  check_file <- grepl(">", read_file[1])

  if (check_file == FALSE){
    stop(paste0("Please provide FASTA formated file"))
  }
  return(check_file)
}


#' For checking signalp version
#'
#' This function check check signalp version
#'.
#' @param path for path of the programe.
#'
#' @export
#' @return provide version of the signalp programe
#'
#'
#'  @examples
#'
#'check_signalp_version(path = NULL)
#'

### check signalp version
check_signalp_version <- function(path = NULL){
  system2(path, "-version", stdout = "version.txt")
  sinalp_version_text <- readLines("version.txt")
  path <- grepl("5", sinalp_version_text)
  if (path == TRUE) {
    version = 5
  }
  path <- grepl("3", sinalp_version_text)
  if (path == TRUE) {
    version = 3
  }

  file.remove("version.txt")
  return(version)
}


#' For check signalp program and version match
#'
#' This function check signalp program and version match
#'.
#' @param program_version for version of the programe.
#' @param input_version Provided input version by user
#' @export
#' @return return logical condition
#'
#'
#'  @examples
#'
#'check_signalp_version_match(program_version=NULL, input_version= NULL)
#'

### check signalp program and version match
check_signalp_version_match <- function(program_version=NULL, input_version= NULL){
  if (program_version != input_version){
    stop('SiganlP program path and input signalp version path is not matching', call. = FALSE)
  }
  return(TRUE)
}


#' modifying filenames and set default names
#'
#' This function modifying filenames and set default names
#'.
#' @param file for input filename.
#' @param programe_name Programe name
#' @param concat_string used for concatination of string in to filename
#' @export
#' @return return output filename
#'
#'  @examples
#'
#'modify_filename(file = NULL, programe_name = NULL, concat_string = ".")
#'


## modifying filenames and set default names
modify_filename <- function(file = NULL, programe_name = NULL, concat_string = ".") {
  filename <- basename(file)
  file_without_ext <- tools::file_path_sans_ext(filename)
  default_fileout_name <- paste0(file_without_ext, concat_string, programe_name)
  return(default_fileout_name)
}



#' Preparation of input data to signalp
#'
#' This function is used for preparation of input data to signalp
#'.
#' @param file for input filename.
#' @param version version of signalp programe.
#' @export
#' @return return input filename
#'
#'
#'  @examples
#'prepare_inputdata_signalp(version = NULL, file = NULL)
#'
#'

### preparation of input data to signalp
prepare_inputdata_signalp <- function(version = NULL, file = NULL){
  if (version == 3) {
    sequence <- readLines(file)
    Id <- grep(">", sequence, value = TRUE)
    sequence_count_cutoff <- 300
    if (length(Id) > sequence_count_cutoff){
      destdir <- dir.create(file.path("signalpdir"), recursive = TRUE)
      file_parse <- parse_file(path_in = file, path_out = "signalpdir/parse.fasta",num_seq = sequence_count_cutoff, trim = FALSE,trunc = 0, id = TRUE)
      file_list <- file_parse$file_list
      #path <- file.path(getwd())
      #get_path<- file.path(getwd(),"signalpdir")
      #setwd(get_path)
      #temp = paste(list.files(pattern="*.fa", full.names = TRUE))
      #numbers = as.numeric(regmatches(temp, regexpr("[0-9]+", temp)))
      #temp <- temp[order(numbers)]
    } else {
      file_list <- c(file)
    }
  }
  if (version == 5){
    file_list <-c(file)
  }
  return(file_list)
}

#' Running signalp
#'
#' This function is used for Running signalp
#'
#' @param programe_path Full path of programe.
#' @param version version of signalp programe.
#' @param input_file.list input files
#' @param org.type For type of organism
#' @param output output file name
#' @export
#' @return return output filename
#'
#'
#'  @examples
#'run_signalp(programe_path =NULL, version = NULL, input_file.list = NULL, org.type = NULL, output = NULL)
#'

### Running signalp
run_signalp <- function(programe_path =NULL, version = NULL, input_file.list = NULL, org.type = NULL, output = NULL){
  if (version == 5) {
    signalp.command  <- c(programe_path, "-fasta", input_file.list)
    system2(signalp.command, org.type)
    ###need to add something
    input_file <- input_file.list
    output <- modify_filename(file = input_file, programe_name = "summary.signalp5", concat_string = "_")
  } else if(version == 3){
    file_count <- length(input_file.list)
    if (file_count > 1){
      for (i in 1:file_count){
        signalp.command  <- c(programe_path ,"-t", org.type, "-f","summary","-trunc", "70", input_file.list[i])
        utils::write.table(system2(signalp.command, stdout = TRUE), file = output, row.names = FALSE, append = TRUE, col.names = FALSE, sep = ", ")
      }
    } else {
      signalp.command  <- c(programe_path ,"-t", org.type, "-f","summary", "-trunc", "70", input_file)
      system2(signalp.command, stdout = output)
    }
  }
  return(output)
}


#' Dataframe for signalp
#'
#' This function is used for producing Dataframe for signalp
#'
#' @param signalp3_outname signalp version 3 output filename.
#' @param original_input orginal input file name.
#' @export
#' @return return Datframe
#'
#'
#'  @examples
#'summary_signalp3(signalp3_outname = NULL, original_input = NULL)
#'
#'

### Dataframe for signalp
####signalp3 summary ###
summary_signalp3 <- function(signalp3_outname = NULL, original_input = NULL){
  out_file <- readLines(signalp3_outname)
  source_file <- seqinr::read.fasta(original_input)
  ID <- seqinr::getName(source_file)
  geneId <- data.frame(ID)
  Dscore_extract <-grep("  D", out_file ,value = TRUE)
  Dscore_extract <- unlist(strsplit(Dscore_extract, "\""))
  Dscore_extract <- utils::read.table(text = Dscore_extract)
  probabilty <- grep("Signal peptide probability", out_file, value = TRUE)
  prediction <- grep("Prediction:", out_file, value = TRUE)
  prediction <- gsub("Prediction:", "",prediction)
  Max_cleavage_site_probability <- grep("Max cleavage site probability:", out_file, value = TRUE)
  Max_cleavage_site_probability <- gsub("Max cleavage site probability:", "", Max_cleavage_site_probability)
  probabilty <- unlist(strsplit(probabilty, "\""))
  probabilty <- utils::read.table(text = probabilty)
  Dscore_extract$V6 <- probabilty$V4
  Dscore_extract$V7 <- prediction
  Dscore_extract$V8 <- Max_cleavage_site_probability
  colnames(Dscore_extract) <- c("Measure","Position","Value","Cutoff","signal_Peptide", "HMM_probability","HMM_prediction", "cleavage_site_probability")
  summary_signal_out <- cbind(geneId, Dscore_extract)
  knitr::kable(summary_signal_out)
  return(summary_signal_out)
}



#' Dataframe for signalp
#'
#' This function is used for producing Dataframe for signalp
#'
#' @param signalp5out signalp version 5 output filename.
#' @export
#' @return return Datframe
#'
#'
#'  @examples
#'siganlp5_summary(signalp5out = NULL)
#'
#### signalp5 summry
siganlp5_summary <- function(signalp5out = NULL){
  siganlp5_out_file <- readLines(signalp5out)
  sigpeptide_extract <- grep("SP", siganlp5_out_file, value = TRUE)
  summary_data <- utils::read.table(text = sigpeptide_extract)
  colnames(summary_data) <- c("ID","Prediction","SP(Sec/SPI)", "OTHER",	"CSPosition","v1", "v2","v3", "v4", "v5")
  #summary_data <- subset(summary_data, summary_data$v5 > 0.9)
  summary_data$CSPosition <- paste(summary_data$v1, summary_data$v2, summary_data$v3, summary_data$v4, summary_data$v5)
  drops <- c("v1", "v2","v3", "v4", "v5")
  summary_signal_out <- summary_data[ , !(names(summary_data) %in% drops)]
  original_file <- seqinr::read.fasta(input_file)
  complete_seq <- original_file[seqinr::getName(original_file) %in% summary_signal_out$ID]
  #write.fasta(sequences = complete_seq, names = names(complete_seq), file.out = signal_pred)
  return(summary_signal_out)
}

#### targetp section ########

#'input data preparation
#'
#' This function is used for preparation of input data
#'
#' @param file_input input file name
#' @param dir_name For storing parsing files
#' @export
#' @return file names
#'
#'
#'  @examples
#'prepare_inputdata(file_input = NULL, dir_name = NULL)
#'
##input data preparation ##
prepare_inputdata <- function(file_input = NULL, dir_name = NULL){
  sequence <- readLines(file_input)
  Id <- grep(">", sequence, value = TRUE)
  dirname_with_filename <- paste0(dir_name, "/", "parse.fasta")
  sequence_count_cutoff <- 300
  if (length(Id) > sequence_count_cutoff){
    destdir <- dir.create(file.path(dir_name), recursive = TRUE)
    file_parse <- parse_file(path_in = file_input, path_out = dirname_with_filename, num_seq = sequence_count_cutoff, trim = FALSE,trunc = 0, id = TRUE)
    file_list <- file_parse$file_list
  } else{
    file_list <- list(file_input)
  }
  return(file_list)
}


#' FOR running targetp
#'
#' This function is used forrun targetp
#'
#' @param programe_path Programe complete path
#' @param input_file.list Input filenames
#' @param org.group Organism grop name
#' @param output output filenam
#' @export
#' @return Resulted output file name
#'
#'
#'  @examples
#'run_targetp(programe_path =NULL, input_file.list = NULL, org.group = NULL, output = NULL)
#'
##### run targetp

run_targetp <- function(programe_path =NULL, input_file.list = NULL, org.group = NULL, output = NULL) {
  file_count <- length(input_file.list)
  if (file_count > 1) {
    for (i in 1:file_count) {
      target.command <- c(programe_path, org.group, input_file.list[i])
      utils::write.table(system2(target.command, stdout = TRUE), file = output, row.names = FALSE, append = TRUE,col.names = FALSE, sep = ",")
    }
  } else {
    target.command <- c(programe_path, org.group, input_file.list)
    system2(target.command, stdout = output)
  }
  return(output)
}

#' Dataframe targetp
#'
#' This function is used for producing Dataframe targetpp
#'
#' @param targetp_outfile Resulted targetp output filename
#' @param organism organism name (Plant or nonplant)
#' @param error If \code{TRUE}, throw an error if programs is not installed.
#' @param verbose If \code{TRUE}, print progress reports.
#' @export
#' @return Resulted output file name
#'
#'
#'  @examples
#' targetp_summary(targetp_outfile = NULL)
#'
#'
### dataframe targetp ##

targetp_summary <- function(targetp_outfile = NULL, organism = NULL, error = TRUE, verbose = FALSE) {
  targetp_file <- readLines(targetp_outfile)
  extra <- data.frame(beg=which(grepl("### targetp v1.1 prediction results ##################################", targetp_file)), end=which(grepl("RC",targetp_file)))
  extra_text <- vector("list", 1)
  for (i in 1:nrow(extra)) {
    extra_text[[i]] <- extra$beg[i]:extra$end[i]
  }
  extra_text <- unlist(extra_text)
  start <- list()
  for (num in 1:length(extra_text)) {
    start[num] <- targetp_file[extra_text[num]]
  }
  out_file_without_extra <- targetp_file[which(!targetp_file %in% start)]
  out_file1 <- grep("---", out_file_without_extra, value = TRUE)
  out_file2 <- grep("cutoff", targetp_file, value = TRUE)
  out_file1 <- gsub("\\s",",",out_file_without_extra)
  out_file2<-strsplit(out_file1,"[, ]+")

  targetp_summary <- data.frame(matrix( nrow = length(out_file2), ncol = 2))
  colnames(targetp_summary) <- c("NAME","LOCALIZATION")
  if (organism == "-P") {
    for (i in 1:length(out_file2)) {
      targetp_summary[i,2] <- out_file2[[i]][7]
      targetp_summary[i,1] <- out_file2[[i]][1]
    }
  } else {
    for (i in 1:length(out_file2)) {
      targetp_summary[i,2] <- out_file2[[i]][6]
      targetp_summary[i,1] <- out_file2[[i]][1]
    }
  }
  targetp_summary <- targetp_summary[!is.na(targetp_summary$LOCALIZATION), ]
  return(targetp_summary)
}



#' FOR running tmhmm
#'
#' This function is used forrun tmhmm
#'
#' @param tmhmm.path Programe complete path
#' @param input_file.list Input filenames
#' @param output output filename
#'
#' @export
#' @return Resulted output file name
#'
#'
#'  @examples
#'run_tmhmm(tmhmm.path = NULL, input_file.list = NULL, output = NULL)
#'
########## tmhmm #####

run_tmhmm <- function(tmhmm.path = NULL, input_file.list = NULL, output = NULL) {
  file_count <- length(input_file.list)
  if (file_count > 1) {
    for (i in 1:file_count) {
      tmhmm.command<- c(tmhmm.path,input_file.list[i])
      utils::write.table(system2(tmhmm.command, stdout = TRUE), file = output, row.names = FALSE, append = TRUE, col.names = FALSE, sep = ", ")
    }
  } else {
    tmhmm.command <- c(tmhmm.path, input_file.list)
    system2(tmhmm.command, stdout = output)
  }
  return(output)
}


#' Dataframe tmhmm
#'
#' This function is used for producing Dataframe tmhmm
#'
#' @param tmhmm_outfile Resulted tmhmm output filename
#' @export
#' @return Resulted output file name
#'
#'
#'  @examples
#' tmhmm_summary(tmhmm_outfile = NULL)
#'
#### tmhmm dataframe
tmhmm_summary <- function(tmhmm_outfile = NULL) {
  tmhmm_file <- readLines(tmhmm_outfile)
  extract <- grep("predicted TMHs",tmhmm_file,value = TRUE)
  extra_extract <- gsub("Number of predicted TMHs","",extract)
  hash_extract <- gsub("#","",extra_extract)
  out_file <-unlist(strsplit(hash_extract,":"))
  tmhmm_summary <-data.frame(matrix(out_file, ncol = 2, byrow = TRUE))
  #tmhmm_summary <-data.frame(matrix(nrow = nrow(out_file), ncol = 2))
  colnames(tmhmm_summary)<- c("ID","NUMBER_OF_PREDICTED_TMHS")
  return(tmhmm_summary)
}


#' MAFFT Alignment
#'
#' This function is perfoms multiple sequence alignment when alignment file is not provided by user
#'
#' @param alignment.file MAFFT alignment resulted file
#' @param sequences PATT_REG sequence reading and stored in sequences
#' @param mafft.path complete path of MAFFT
#' @param num.threads Number of threads to be used by MAFFT.
#' @export
#' @return Resulted output file name
#'
#'
#'  @examples
#' get_mafft(alignment.file = NULL, sequences = NULL, mafft.path = NULL, num.threads = NULL)
#'
## get MAFFT alignment file
get_mafft <- function(alignment.file = NULL, sequences = NULL, mafft.path = NULL, num.threads = NULL){
  mafft.out.name <- c("MAFFT.fasta")
  if (is.null(alignment.file)== T) {
    file.name <- c("REGEX.fasta")
    seqinr::write.fasta(sequences = seqinr::getSequence(sequences), names=seqinr::getName(sequences), file.out = file.name)
    mafft.path <- external_programe_existance(programe = "mafft", path = mafft.path)
    mafft.out.name <- run_mafft(mafft_path = mafft.path, num_threads = num.threads, file_name = file.name, mafft.out_file = mafft.out.name)
  } else {
    mafft.out.name <- alignment.file
  }
  return(mafft.out.name)
}

## mafft alignment ###

#' Running MAFFT Alignment
#'
#' This function is perfoms multiple sequence alignment
#'
#' @param mafft_path complete path of MAFFT
#' @param num_threads Number of threads to be used by MAFFT.
#' @param file_name PATT_REG sequence
#' @param mafft.out_file MAFFT ouputput filename
#' @export
#' @return Resulted output file name
#'
#'
#'  @examples
#' run_mafft(mafft_path = NULL, num_threads = NULL, file_name = NULL, mafft.out_file = NULL)
#'
run_mafft <- function(mafft_path = NULL, num_threads = NULL, file_name = NULL, mafft.out_file = NULL){
  file_read <- read.fasta(file_name)
  file_length <- length(file_read)
  if (file_length < 1000){
    mafft.command <- c(mafft_path, "--legacygappenalty", "--genafpair", "--maxiterate", "1000", "--thread",num_threads, "--quiet", file_name)
  } else {
    mafft.command <- c(mafft_path, "--parttree", "--thread", num_threads, file_name)
  }
  system2(mafft.command, stdout = mafft.out_file)
  return(mafft.out_file)
}


## HMMER

#' HMM build
#'
#' This function is build model by using MAFFT alignment sequence file
#'
#' @param hmm_path  complete path of HMMER
#' @param seed The seed to used with HMMER commands. Set this to get the same output each time
#' @param mafft_out MAFFT alignment output file name
#'
#' @export
#' @return Resulted output file name
#'
#'  @examples
#' hmm_build(hmm_path = NULL, seed = NULL,  mafft_out = NULL)
#'
hmm_build <- function(hmm_path = NULL, seed = NULL,  mafft_out = NULL){
  hmmbuild.out <- c("hmmbuild.hmm")
  hmm_path <- external_programe_existance(programe = "hmmbuild", path = hmm_path)
  hmmbuild_command <- c(hmm_path, "--amino", "--seed", seed, hmmbuild.out,  mafft_out)
  system2(hmmbuild_command, stdout = F)
  return(hmmbuild.out)
}

### hmmpress

#' Prepare a profile database for hmmscan
#'
#' This function is Constructs binary compressed datafiles for hmmscan, starting from a profile database hmmfile in standard HMMER3 format.
#'
#' @param path_hmm  complete path of HMMER
#'
#' @export
#' @return Resulted output file name
#'
#'  @examples
#' hmm_press <- function(path_hmm = NULL)
#'

hmm_press <- function(path_hmm = NULL){
  hmmbuild.out <- c("hmmbuild.hmm")
  path_hmm <- external_programe_existance(programe = "hmmpress", path = path_hmm)
  hmmpress_command <- c(path_hmm, "-f", hmmbuild.out)
  system2(hmmpress_command)
  return(TRUE)
}

## hmm search

#' search profile(s) against a sequence database
#'
#' hmmsearch is used to search one or more profiles against a sequence database.
#'
#' @param hmmsearch_path  complete path of HMMER
#' @param seed The seed to used with HMMER commands. Set this to get the same output each time.
#' @param hmm.tresh Set the bit score cutoff for the per-sequence ranked hit list to a real number. (Default = 0)
#' @param original.seq The absolute path for the original six-frame translation FASTA file
#'
#' @export
#' @return Resulted output file name
#'
#'  @examples
#' hmm_search(hmmsearch_path = NULL, seed = NULL,hmm.tresh = NULL, original.seq = NULL)
#'
hmm_search <- function(hmmsearch_path = NULL, seed = NULL,hmm.tresh = NULL, original.seq = NULL){
  hmmbuild.out <- c("hmmbuild.hmm")
  hmmsearch.out <- c("hmmsearch.txt")
  hmmsearch_path <- external_programe_existance(programe = "hmmsearch", path = hmmsearch_path)
  hmmsearch_command <- c(hmmsearch_path,
                         "--seed", seed,
                         "-T", hmm.tresh,
                         "--tblout", hmmsearch.out,
                         hmmbuild.out,
                         original.seq)
  system2(hmmsearch_command)
  return(hmmsearch.out)
}

### HMMER for windows

### Converts FASTA to STOCKHOLM for windows

#' Converts FASTA to STOCKHOLM
#'
#' The HMMER binary version for windows does not recognize FASTA files to build the
#' hmm profile. This function creates a STOCKHOLM file readable by HMMER 3.0
#'
#' @param fasta.file FASTA object to be converted
#'
#' @return TRUE/FALSE
#'
#' @keywords general
#'
fasta_to_stockholm <- function(fasta.file){
  stock.name <- gsub(fasta.file, pattern = ".fasta", replacement = ".stockholm")
  seq <- seqinr::read.fasta(fasta.file)
  seq.seq <- lapply(seqinr::getSequence(seq), function (x) (paste0(x,collapse="")))
  seq.names <- seqinr::getName(seq)
  seq.final <- list()
  for (i in 1:length(seq)){
    seq.final[[i]] <- paste(seq.names[[i]],seq.seq[[i]], sep = " ")
  }
  seq.final <- unlist(seq.final)
  writeLines(c("# STOCKHOLM 1.0", seq.final,"//"), con = stock.name, sep = "\n")
}

## HMM build

#' HMM build
#'
#' This function is build model by using MAFFT alignment sequence file
#'
#' @param hmm_path  complete path of HMMER
#' @param seed The seed to used with HMMER commands. Set this to get the same output each time.
#' @param mafft.out.name mafft alignment filename.
#'
#' @return Resulted output file name
#'
#' @keywords general
#'
hmmbuild <- function(hmm_path = NULL, seed = NULL, mafft.out.name = NULL){
  hmmbuild.out <- c("hmmbuild.hmm")
  fasta_to_stockholm(fasta.file = mafft.out.name)
  stock.name <- gsub(mafft.out.name, pattern = ".fasta", replacement = ".stockholm")
  hmmbuild_path <- external_programe_existance(programe = "hmmbuild.exe", path = hmmbuild_path)
  hmmbuild_command <- c(hmmbuild_path,
                        "--amino",
                        "--seed", seed,
                        hmmbuild.out,
                        stock.name)
  system2(hmmbuild_command, stdout = F)
  return(hmmbuild.out)
}

## HMM Press

#' Prepare a profile database for hmmscan
#'
#' This function is Constructs binary compressed datafiles for hmmscan, starting from a profile database hmmfile in standard HMMER3 format.
#'
#' @param hmmpress_path  complete path of HMMER
#'
#' @return Resulted output file name
#'
#' @keywords general
#'
hmmpress <- function(hmmpress_path = NULL) {
  hmmbuild.out <- c("hmmbuild.hmm")
  hmmpress_path <- external_programe_existance(programe = "hmmpress.exe", path = hmmpress_path)
  hmmpress_command <- c(hmmpress_path,
                        "-f",
                        hmmbuild.out)

  system2(hmmpress_command)
  return(TRUE)
}

##  HMM search

#' search profile(s) against a sequence database
#'
#' hmmsearch is used to search one or more profiles against a sequence database.
#'
#' @param hmmsearch_path  complete path of HMMER
#' @param seed The seed to used with HMMER commands. Set this to get the same output each time.
#' @param hmm.tresh Set the bit score cutoff for the per-sequence ranked hit list to a real number. (Default = 0)
#' @param original.seq The absolute path for the original six-frame translation FASTA file
#'
#' @return Resulted output file name
#'
#' @keywords general
#'

hmmsearch <- function(hmmsearch_path = NULL, seed = NULL,hmm.tresh = NULL, original.seq = NULL){
  hmmbuild.out <- c("hmmbuild.hmm")
  hmmsearch.out <- c("hmmsearch.txt")
  hmmsearch_path <- external_programe_existance(programe = "hmmsearch.exe", path = hmmsearch_path)
  hmmsearch_command <- c(hmmsearch_path,
                         "--seed", seed,
                         "-T", hmm.tresh,
                         "--tblout", hmmsearch.out,
                         hmmbuild.out,
                         original.seq)

  system2(hmmsearch_command)
  return(hmmsearch.out)
}


### reading hmm results

#' Reading HMM resulted file
#'
#' hmmsearch is used to search one or more profiles against a sequence database.
#'
#' @param original_seq The absolute path for the original six-frame translation FASTA file
#' @param hmmsearch.out hmmsearch output file
#' @param hmmbuild.out hmmbuild output file name
#' @param regex.seq pattern search function returns the sequence stored in PATT_REG variable.
#' @param save.alignment alignment file
#'
#' @export
#' @return returns total sequence (HMM, PATT_REG)
#'
#'  @examples
#' reading_hmm.results(original_seq = NULL, hmmsearch.out = NULL, hmmbuild.out = NULL, regex.seq = NULL, save.alignment = FALSE)
#'

reading_hmm.results <- function(original_seq = NULL, hmmsearch.out = NULL, hmmbuild.out = NULL, regex.seq = NULL, save.alignment = FALSE){
  hmm.hits <- utils::read.delim(hmmsearch.out, comment.char = "#", header = F, sep = "", blank.lines.skip = T, stringsAsFactors = F)[,1]
  if (Sys.info()[['sysname']] %in% "Windows"){
    hmm.table <- utils::read.table(hmmbuild.out, blank.lines.skip = T, skip = 14, sep = "", fill = T, stringsAsFactors = F)
  }
  hmm.table <- utils::read.table(hmmbuild.out, blank.lines.skip = T, skip = 16, sep = "", fill = T, stringsAsFactors = F)
  total.seq <- seqinr::read.fasta(original_seq)
  hmm.seq <- total.seq[seqinr::getName(total.seq) %in% hmm.hits]
  if (save.alignment == T){
    alin.seq <- seqinr::read.fasta(mafft.out.name)
    total.seq <- list(alin.seq, regex.seq, hmm.seq, hmm.table)
    names(total.seq) <- c("Alignment","REGEX","HMM","HMM_Table")
    cat(paste0("\nHMM search done. \n---\n\nTotal of sequences found in PATT_REG: ", length(total.seq[[2]]), "\n"))
    cat(paste0("Total of sequences found in HMM: ",  length(total.seq[[3]]), "\n"))
    cat(paste0("Total of redundant hits: ", sum(duplicated(unlist(lapply(total.seq[c(2,3)], function (x) seqinr::getName(x))))),"\n"))
    cat(paste0("Number of motif candidates: ",length(unique(unlist(lapply(total.seq[c(2,3)], function (x) seqinr::getName(x)))))),"\n")
    cat("---\n")
    cat("Alignment saved (as first list element).\n")
    cat("---\n")
  } else {
    total.seq <- list(regex.seq, hmm.seq, hmm.table)
    names(total.seq) <- c("REGEX","HMM","HMM_Table")
    cat(paste0("\nHMM search done. \n---\n\nTotal of sequences found in PATT_REG: ", length(total.seq[[1]]), "\n"))
    cat(paste0("Total of sequences found in HMM: ",  length(total.seq[[2]]), "\n"))
    cat(paste0("Total of redundant hits: ", sum(duplicated(unlist(lapply(total.seq[c(1,2)], function (x) seqinr::getName(x))))),"\n"))
    cat(paste0("Number of  motif candidates: ",length(unique(unlist(lapply(total.seq[c(1,2)], function (x) seqinr::getName(x)))))),"\n")
    cat("---\n")
  }
  cat("hmm.search finished!\n")
  cat("---\n")
  return(total.seq)
}

## Generating HMM fasta file

#' Generating HMM resulted fasta file
#'
#' This fuction is used for generating non redundant motif candidates fasta file.
#'
#' @param total_seq sequence containing (PATT_REG, HMM, HMM_table)
#' @param original.seq The absolute path for the original six-frame translation FASTA file.
#'
#' @export
#' @return returns ouput fasta file name
#'
#'  @examples
#' HMM_fasta(total_seq = NULL, original.seq = NULL)
#'
HMM_fasta <- function(total_seq = NULL, original.seq = NULL) {
  motif.out <- list()
  if (length(total_seq) == 4){
    hmm.result <- total_seq[names(total_seq) %in% c("REGEX","HMM","HMM_Table")]
  }
  if (length(total_seq) == 3) {
    summ.data <- c(total_seq[[1]],total_seq[[2]])
    consensus.seq <- summ.data[!duplicated(seqinr::getName(summ.data))]
    motif.out$consensus.sequences <- consensus.seq
  } else {
    consensus.seq <- hmm.result
  }
  #setwd(original.dir)
  output_filename <- modify_filename(file = original.seq, programe_name = "HMM.fasta", concat_string = "_")
  seqinr::write.fasta(sequences = motif.out$consensus.sequences, names = names(motif.out$consensus.sequences), file.out = output_filename)
  return(output_filename)
}
