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
#'        Default is "euclidean".
#' @return A list containing the selected distance metric and diagnostic information:
#'   \item{distance}{The calculated distance value}
#'   \item{completeness}{Proportion of complete (non-NA) pairs}
#'   \item{valid}{Logical indicating if the metric was successfully computed}
#' @examples
#' seq1 <- c(1, 2, NA, 4, 5)
#' seq2 <- c(2, NA, 4, 5, 6)
#' # Calculate Euclidean distance (default)
#' result <- perform_distance(seq1, seq2)
#' # Calculate DTW distance
#' result <- perform_distance(seq1, seq2, method = "dtw")
#' print(result$distance)
#' print(result$completeness)
#' print(result$valid)
#' @importFrom dtw dtw
#' @importFrom stats cor complete.cases na.omit approx
#' @export
perform_distance <- function(seq1, seq2, method = "euclidean") {
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
    
    # Calculate selected metric with fallback to R implementations
    tryCatch({
        result$distance <- switch(method,
            "dtw" = {
                alignment <- dtw::dtw(seq1, seq2, 
                                    step.pattern = dtw::asymmetric,
                                    keep.internals = TRUE)
                alignment$distance / n_complete
            },
            "euclidean" = {
                # Try C++ first, fallback to R
                tryCatch({
                    euclidean_distance_cpp(seq1, seq2)
                }, error = function(e) {
                    # R implementation: Sum of Euclidean distances (matches v0.0.1.4)
                    if (n_complete == 0) NA_real_ else sum(sqrt((complete_seq1 - complete_seq2)^2)) / n_complete
                })
            },
            "manhattan" = {
                # Try C++ first, fallback to R
                tryCatch({
                    manhattan_distance_cpp(seq1, seq2)
                }, error = function(e) {
                    # R implementation: MAE
                    if (n_complete == 0) NA_real_ else mean(abs(complete_seq1 - complete_seq2))
                })
            },
            "correlation" = 1 - abs(cor(complete_seq1, complete_seq2)),
            "area" = {
                # Try C++ first, fallback to R
                tryCatch({
                    area_distance_cpp(seq1, seq2)
                }, error = function(e) {
                    # R implementation: simple interpolation and area difference
                    interp_seq1 <- na.omit(approx(seq1, n = length(seq1))$y)
                    interp_seq2 <- na.omit(approx(seq2, n = length(seq2))$y)
                    if (length(interp_seq1) == 0 || length(interp_seq2) == 0) {
                        NA_real_
                    } else {
                        mean(abs(interp_seq1 - interp_seq2))
                    }
                })
            },
            "hausdorff" = {
                # Try C++ first, fallback to R
                tryCatch({
                    hausdorff_distance_cpp(seq1, seq2)
                }, error = function(e) {
                    # R implementation: basic Hausdorff distance
                    clean_seq1 <- complete_seq1
                    clean_seq2 <- complete_seq2
                    if (length(clean_seq1) == 0 || length(clean_seq2) == 0) {
                        NA_real_
                    } else {
                        max_dist <- 0
                        for (val1 in clean_seq1) {
                            min_dist <- min(abs(val1 - clean_seq2))
                            max_dist <- max(max_dist, min_dist)
                        }
                        for (val2 in clean_seq2) {
                            min_dist <- min(abs(val2 - clean_seq1))
                            max_dist <- max(max_dist, min_dist)
                        }
                        max_dist
                    }
                })
            }
        )
        result$valid <- TRUE
    }, error = function(e) {
        warning(method, " calculation failed: ", e$message)
    })
    
    return(result)
}