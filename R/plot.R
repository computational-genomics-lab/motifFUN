#' Plots the relative frequencies of each position for hmmsearch table.
#'
#' This function plots the results from hmm.search as a pointplot with amino acids in the x axis and the relative frequency of each amino acid in the y axis
#' @param hmm_data The HMM profile table resulting from hmm.search
#' @keywords PATT_REG  plot
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' hmm.plot(hmm_data = NULL)
#' }

hmm.plot <- function (hmm_data = NULL) {
  hmm <- hmm_data
  colnames(hmm) <- hmm[1,]
  hmm <- hmm[-(1:5),]
  hmm[,1] <- suppressWarnings(as.numeric(hmm[,1]))
  hmm <- hmm[hmm[,1]%%1==0, ]
  hmm <- suppressWarnings(hmm[!is.na(as.numeric(hmm[,2])), ])
  rownames(hmm) <- hmm[,1]
  hmm <- hmm[,-1]
  hmm <- data.frame(sapply(hmm, function(x) as.numeric(as.character(x))), stringsAsFactors = F)
  hmm.sums <- apply(hmm,2,function (x) max(x)/x)
  hmm.sums <- apply(hmm.sums,2,function (x) x/sum(x))
  hmm.sums <- apply(hmm.sums, 2, function (x) x - mean(as.numeric(as.character(x))))
  hmm.sums[hmm.sums < 0] <- 0
  hmm.m <- reshape2::melt(t(hmm.sums))
  colnames(hmm.m) <- c("element","position","bits")
  hmm.m$bits <- as.numeric(as.character(hmm.m$bits))
  logo <- ggplot2::ggplot(hmm.m, aes(x = position, y = bits)) +
    ggplot2::geom_point(aes(color = element), size = 2) + ggplot2::geom_text(ggplot2::aes_string(label="element", size="bits"), position='stack') +
    viridis::scale_fill_viridis(discrete=TRUE) + ggplot2::theme_bw() + ggplot2::ylab("Relative Frequency (bits)") + ggplot2::guides(fill=FALSE)
  return(logo)
}
