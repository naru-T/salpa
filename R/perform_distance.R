#' Calculate Distance Metric Between Two Sequences
#'
#' This function takes two sequences and calculates a selected distance metric
#' commonly used in trajectory analysis. Each metric handles NA values appropriately:
#' 1. DTW - uses step patterns that can skip NA values
#' 2. Euclidean - computed only on complete pairs
#' 3. Manhattan - computed only on complete pairs
#' 4. Correlation - uses pairwise complete observations
#' 5. Area - interpolates NA values using linear interpolation
#' 6. Hausdorff - considers only non-NA values in both sequences
#'
#' @param seq1 A numeric vector representing the first sequence.
#' @param seq2 A numeric vector representing the second sequence.
#' @param method Character string specifying which distance metric to calculate.
#'        Options are "dtw", "euclidean", "manhattan", "correlation", "area", "hausdorff".
#'        Default is "dtw".
#' @return A list containing the selected distance metric and diagnostic information:
#'   \item{distance}{The calculated distance value}
#'   \item{completeness}{Proportion of complete (non-NA) pairs}
#'   \item{valid}{Logical indicating if the metric was successfully computed}
#' @examples
#' seq1 <- c(1, 2, NA, 4, 5)
#' seq2 <- c(2, NA, 4, 5, 6)
#' # Calculate DTW distance
#' result <- perform_distance(seq1, seq2)
#' # Calculate Euclidean distance
#' result <- perform_distance(seq1, seq2, method = "euclidean")
#' print(result$distance)
#' print(result$completeness)
#' print(result$valid)
#' @importFrom dtw dtw
#' @importFrom stats cor complete.cases na.omit approx
#' @export
perform_distance <- function(seq1, seq2, method = "dtw") {
    # Input validation
    if (length(seq1) != length(seq2)) {
        stop("Sequences must have equal length")
    }
    
    # Available methods
    available_methods <- c("dtw", "euclidean", "manhattan", "correlation", "area", "hausdorff")
    
    # Validate method parameter
    if (!method %in% available_methods) {
        stop("Invalid method: ", method, 
             "\nAvailable methods are: ", paste(available_methods, collapse = ", "))
    }
    
    # Create data frame and identify complete cases
    df <- data.frame(seq1, seq2)
    complete_cases <- complete.cases(df)
    n_complete <- sum(complete_cases)
    n_total <- length(seq1)
    completeness <- n_complete / n_total
    
    # Initialize result
    result <- list(
        distance = NA_real_,
        completeness = completeness,
        valid = FALSE
    )
    
    # Only proceed if we have sufficient complete pairs
    min_pairs <- if(method == "correlation") 2 else 1
    if (n_complete < min_pairs) {
        warning("Insufficient complete pairs for ", method, " calculation")
        return(result)
    }
    
    # Get complete pairs for point-wise metrics
    complete_seq1 <- df$seq1[complete_cases]
    complete_seq2 <- df$seq2[complete_cases]
    
    # Calculate selected metric
    tryCatch({
        result$distance <- switch(method,
            "dtw" = {
                alignment <- dtw::dtw(seq1, seq2, 
                                    step.pattern = dtw::asymmetric,
                                    keep.internals = TRUE)
                alignment$distance / n_complete
            },
            "euclidean" = sum(sqrt((complete_seq1 - complete_seq2)^2)) / n_complete,
            "manhattan" = sum(abs(complete_seq1 - complete_seq2)) / n_complete,
            "correlation" = 1 - abs(cor(complete_seq1, complete_seq2)),
            "area" = {
                idx <- seq_along(seq1)
                seq1_interp <- approx(idx[!is.na(seq1)], seq1[!is.na(seq1)], idx)$y
                seq2_interp <- approx(idx[!is.na(seq2)], seq2[!is.na(seq2)], idx)$y
                sum(abs(seq1_interp - seq2_interp)) / n_total
            },
            "hausdorff" = {
                dist_matrix <- outer(complete_seq1, complete_seq2, 
                                   function(x, y) abs(x - y))
                max(c(apply(dist_matrix, 1, min),
                     apply(dist_matrix, 2, min)))
            }
        )
        result$valid <- TRUE
    }, error = function(e) {
        warning(method, " calculation failed: ", e$message)
    })
    
    return(result)
}