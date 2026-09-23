#include <Rcpp.h>
#include <algorithm>
#include <numeric>
#include <vector>
#include <cmath>

using namespace Rcpp;

// Fast Euclidean distance calculation (Sum of distances - matches v0.0.1.4)
// [[Rcpp::export]]
double euclidean_distance_cpp(NumericVector seq1, NumericVector seq2) {
    int n = seq1.size();
    if (n != seq2.size()) {
        stop("Sequences must have equal length");
    }
    
    double sum = 0.0;
    int count = 0;
    
    for (int i = 0; i < n; i++) {
        if (!NumericVector::is_na(seq1[i]) && !NumericVector::is_na(seq2[i])) {
            double diff = seq1[i] - seq2[i];
            sum += sqrt(diff * diff);  // Fixed: Sum of individual Euclidean distances (matches v0.0.1.4)
            count++;
        }
    }
    
    if (count == 0) {
        return NA_REAL;
    }
    
    return sum / count;  // Fixed: Average of Euclidean distances (matches v0.0.1.4 formula)
}

// Fast Manhattan distance calculation
// [[Rcpp::export]]
double manhattan_distance_cpp(NumericVector seq1, NumericVector seq2) {
    int n = seq1.size();
    if (n != seq2.size()) {
        stop("Sequences must have equal length");
    }
    
    double sum = 0.0;
    int count = 0;
    
    for (int i = 0; i < n; i++) {
        if (!NumericVector::is_na(seq1[i]) && !NumericVector::is_na(seq2[i])) {
            sum += std::abs(seq1[i] - seq2[i]);
            count++;
        }
    }
    
    if (count == 0) {
        return NA_REAL;
    }
    
    return sum / count;
}

// Optimized Hausdorff distance calculation with sorted vectors
// [[Rcpp::export]]
double hausdorff_distance_cpp(NumericVector seq1, NumericVector seq2) {
    // Remove NA values and sort for faster searching
    std::vector<double> clean_seq1, clean_seq2;
    
    for (int i = 0; i < seq1.size(); i++) {
        if (!NumericVector::is_na(seq1[i])) {
            clean_seq1.push_back(seq1[i]);
        }
    }
    
    for (int i = 0; i < seq2.size(); i++) {
        if (!NumericVector::is_na(seq2[i])) {
            clean_seq2.push_back(seq2[i]);
        }
    }
    
    if (clean_seq1.empty() || clean_seq2.empty()) {
        return NA_REAL;
    }
    
    // Sort for faster binary search-like operations
    std::sort(clean_seq1.begin(), clean_seq1.end());
    std::sort(clean_seq2.begin(), clean_seq2.end());
    
    double max_dist = 0.0;
    
    // Find max of min distances from seq1 to seq2
    for (double val : clean_seq1) {
        // Binary search for closest value
        auto it = std::lower_bound(clean_seq2.begin(), clean_seq2.end(), val);
        double min_dist = INFINITY;
        
        if (it != clean_seq2.end()) {
            min_dist = std::min(min_dist, std::abs(val - *it));
        }
        if (it != clean_seq2.begin()) {
            --it;
            min_dist = std::min(min_dist, std::abs(val - *it));
        }
        
        max_dist = std::max(max_dist, min_dist);
    }
    
    // Find max of min distances from seq2 to seq1
    for (double val : clean_seq2) {
        // Binary search for closest value
        auto it = std::lower_bound(clean_seq1.begin(), clean_seq1.end(), val);
        double min_dist = INFINITY;
        
        if (it != clean_seq1.end()) {
            min_dist = std::min(min_dist, std::abs(val - *it));
        }
        if (it != clean_seq1.begin()) {
            --it;
            min_dist = std::min(min_dist, std::abs(val - *it));
        }
        
        max_dist = std::max(max_dist, min_dist);
    }
    
    return max_dist;
}

// Optimized area distance calculation with efficient interpolation
// [[Rcpp::export]]
double area_distance_cpp(NumericVector seq1, NumericVector seq2) {
    int n = seq1.size();
    if (n != seq2.size()) {
        stop("Sequences must have equal length");
    }
    
    // Helper function to interpolate a sequence
    auto interpolate_sequence = [](const NumericVector& seq) -> std::vector<double> {
        int n = seq.size();
        std::vector<double> result(n);
        
        // First pass: copy non-NA values
        for (int i = 0; i < n; i++) {
            result[i] = NumericVector::is_na(seq[i]) ? NAN : seq[i];
        }
        
        // Second pass: forward fill
        for (int i = 1; i < n; i++) {
            if (std::isnan(result[i]) && !std::isnan(result[i-1])) {
                result[i] = result[i-1];
            }
        }
        
        // Third pass: backward fill
        for (int i = n-2; i >= 0; i--) {
            if (std::isnan(result[i]) && !std::isnan(result[i+1])) {
                result[i] = result[i+1];
            }
        }
        
        // If still NaN, set to 0
        for (int i = 0; i < n; i++) {
            if (std::isnan(result[i])) {
                result[i] = 0.0;
            }
        }
        
        return result;
    };
    
    std::vector<double> interp_seq1 = interpolate_sequence(seq1);
    std::vector<double> interp_seq2 = interpolate_sequence(seq2);
    
    // Calculate area difference using vectorized operation
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += std::abs(interp_seq1[i] - interp_seq2[i]);
    }
    
    return sum / n;
}

// WOA position update (Mirjalili 2016 correct formula, b=1 fixed)
// rand_idx is 0-indexed
// [[Rcpp::export]]
NumericMatrix woa_update_positions_cpp(
        NumericMatrix whales,
        NumericVector best_pos,
        IntegerVector rand_idx,
        NumericVector r1,
        NumericVector r2,
        NumericVector l_vec,
        NumericVector p_vec,
        double a,
        double b,
        double p_threshold,
        NumericVector lower_bounds,
        NumericVector upper_bounds) {

    int pop_size = whales.nrow();
    int n_dim    = whales.ncol();
    NumericMatrix new_whales(pop_size, n_dim);

    for (int i = 0; i < pop_size; i++) {
        double A = 2.0 * a * r1[i] - a;
        double C = 2.0 * r2[i];

        if (p_vec[i] < p_threshold) {
            if (std::abs(A) < 1.0) {
                // Encircling prey (exploitation)
                for (int j = 0; j < n_dim; j++) {
                    double D = std::abs(C * best_pos[j] - whales(i, j));
                    new_whales(i, j) = best_pos[j] - A * D;
                }
            } else {
                // Search for prey (exploration) using random whale
                int ri = rand_idx[i];  // 0-indexed
                for (int j = 0; j < n_dim; j++) {
                    double D = std::abs(C * whales(ri, j) - whales(i, j));
                    new_whales(i, j) = whales(ri, j) - A * D;
                }
            }
        } else {
            // Bubble-net attacking (spiral, Mirjalili 2016 correct formula)
            double l = l_vec[i];
            for (int j = 0; j < n_dim; j++) {
                double D = std::abs(best_pos[j] - whales(i, j));
                new_whales(i, j) = D * std::exp(b * l) * std::cos(2.0 * M_PI * l) + best_pos[j];
            }
        }

        // Clamp to bounds
        for (int j = 0; j < n_dim; j++) {
            if (new_whales(i, j) < lower_bounds[j]) new_whales(i, j) = lower_bounds[j];
            if (new_whales(i, j) > upper_bounds[j]) new_whales(i, j) = upper_bounds[j];
        }
    }

    return new_whales;
}

// PSO velocity and position update
// r1_flat / r2_flat are flattened pop_size*n_dim random matrices (row-major)
// [[Rcpp::export]]
List pso_update_positions_cpp(
        NumericMatrix positions,
        NumericMatrix velocities,
        NumericMatrix pbest,
        NumericVector gbest,
        NumericVector r1_flat,
        NumericVector r2_flat,
        double w,
        double c1,
        double c2,
        NumericVector lower_bounds,
        NumericVector upper_bounds,
        double v_max) {

    int pop_size = positions.nrow();
    int n_dim    = positions.ncol();
    NumericMatrix new_positions(pop_size, n_dim);
    NumericMatrix new_velocities(pop_size, n_dim);

    for (int i = 0; i < pop_size; i++) {
        for (int j = 0; j < n_dim; j++) {
            int idx = i * n_dim + j;
            double v = w * velocities(i, j)
                     + c1 * r1_flat[idx] * (pbest(i, j) - positions(i, j))
                     + c2 * r2_flat[idx] * (gbest[j]    - positions(i, j));
            // Clamp velocity
            if (v >  v_max) v =  v_max;
            if (v < -v_max) v = -v_max;
            new_velocities(i, j) = v;

            double x = positions(i, j) + v;
            // Clamp position to bounds
            if (x < lower_bounds[j]) x = lower_bounds[j];
            if (x > upper_bounds[j]) x = upper_bounds[j];
            new_positions(i, j) = x;
        }
    }

    return List::create(
        Named("positions")  = new_positions,
        Named("velocities") = new_velocities
    );
}

// Genetic Algorithm population evolution: one generation of
// tournament selection + BLX-alpha crossover + Gaussian mutation + elitism.
//
// @param population  pop_size x n_dim matrix (current generation)
// @param fitness     pop_size fitness values (GA maximises: larger = better)
// @param lower_bounds n_dim lower bounds
// @param upper_bounds n_dim upper bounds
// @param pcrossover  crossover probability (default 0.8)
// @param pmutation   mutation probability per gene (default 0.1)
// @param elitism_k   number of elite individuals copied unchanged (default 5)
// @param tournament_k tournament size for parent selection (default 3)
// @return new population matrix (pop_size x n_dim)
// [[Rcpp::export]]
NumericMatrix ga_evolve_population_cpp(
    NumericMatrix population,
    NumericVector fitness,
    NumericVector lower_bounds,
    NumericVector upper_bounds,
    double pcrossover  = 0.8,
    double pmutation   = 0.1,
    int    elitism_k   = 5,
    int    tournament_k = 3
) {
    int pop_size = population.nrow();
    int n_dim    = population.ncol();

    NumericMatrix new_pop(pop_size, n_dim);

    // ---- Elitism: identify top elitism_k individuals by fitness ----
    std::vector<int> order(pop_size);
    std::iota(order.begin(), order.end(), 0);
    int k = std::min(elitism_k, pop_size);
    std::partial_sort(order.begin(), order.begin() + k, order.end(),
        [&](int a, int b) { return fitness[a] > fitness[b]; });

    for (int i = 0; i < k; i++) {
        for (int j = 0; j < n_dim; j++) {
            new_pop(i, j) = population(order[i], j);
        }
    }

    // ---- Tournament selection (with replacement) ----
    auto tournament_select = [&]() -> int {
        int best = (int)(R::unif_rand() * pop_size);
        if (best >= pop_size) best = pop_size - 1;
        for (int t = 1; t < tournament_k; t++) {
            int challenger = (int)(R::unif_rand() * pop_size);
            if (challenger >= pop_size) challenger = pop_size - 1;
            if (fitness[challenger] > fitness[best]) best = challenger;
        }
        return best;
    };

    // ---- Generate offspring for non-elite slots ----
    for (int i = k; i < pop_size; i++) {
        int p1 = tournament_select();
        int p2 = tournament_select();

        for (int j = 0; j < n_dim; j++) {
            double x1 = population(p1, j);
            double x2 = population(p2, j);
            double child_val;

            // BLX-alpha crossover (alpha = 0.5)
            if (R::unif_rand() < pcrossover) {
                double lo  = std::min(x1, x2);
                double hi  = std::max(x1, x2);
                double d   = hi - lo;
                double ext_lo = std::max(lower_bounds[j], lo - 0.5 * d);
                double ext_hi = std::min(upper_bounds[j], hi + 0.5 * d);
                if (ext_hi > ext_lo) {
                    child_val = ext_lo + R::unif_rand() * (ext_hi - ext_lo);
                } else {
                    child_val = (ext_lo + ext_hi) * 0.5;
                }
            } else {
                child_val = R::unif_rand() < 0.5 ? x1 : x2;
            }

            // Gaussian mutation (sigma = 5% of range)
            if (R::unif_rand() < pmutation) {
                double range = upper_bounds[j] - lower_bounds[j];
                child_val += R::rnorm(0.0, range * 0.05);
                if (child_val < lower_bounds[j]) child_val = lower_bounds[j];
                if (child_val > upper_bounds[j]) child_val = upper_bounds[j];
            }

            new_pop(i, j) = child_val;
        }
    }

    return new_pop;
}
