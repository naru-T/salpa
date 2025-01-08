#' Perform Dynamic Time Warping (DTW) on two sequences
#'
#' This function takes two sequences, removes any rows with NA values, and performs
#' Dynamic Time Warping (DTW) to calculate the alignment distance between the cleaned sequences.
#' The alignment distance is then normalized by the length of the first sequence.
#'
#' @param seq1 A numeric vector representing the first sequence.
#' @param seq2 A numeric vector representing the second sequence.
#' @return A numeric value representing the normalized DTW alignment distance between the two sequences.
#' @examples
#' seq1 <- c(1, 2, 3, 4, 5)
#' seq2 <- c(2, 3, 4, 5, 6)
#' perform_dtw(seq1, seq2)
#' @importFrom dtw dtw
#' @export
perform_dtw <- function(seq1, seq2) {
    df <- data.frame(seq1, seq2)
    df <- na.omit(df)
    seq1_clean <- df$seq1
    seq2_clean <- df$seq2
    alignment <- dtw::dtw(seq1_clean, seq2_clean)$distance
    normalized_alignment <- alignment / length(seq1_clean)
    return(normalized_alignment)
}