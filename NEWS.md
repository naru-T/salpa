# salpa 0.0.2.2

## Bug fixes

- Restored `Depends: R (>= 4.1.0)`. The package uses the native pipe (`|>`).
- `URL` and `BugReports` now point at https://github.com/naru-T/salpa.
- README license section now says GPL-3, matching `DESCRIPTION` and `LICENSE.md`.
- README version heading updated to 0.0.2.2.
- Removed R CMD build metadata (`Packaged`, `Author`, `Maintainer`) from the source tree and cleared executable bits introduced when unpacking the 0.0.2.1 tarball.

# salpa 0.0.2.1

## Performance

- All four optimizers (GA, PSO, WOA, L-BFGS-B) are now faster than v0.0.1.8.0:
  GA 1.26x, PSO 1.24x, WOA 1.08x, L-BFGS-B 1.11x (Nikko 0.5m benchmark, 36 groups).
- GA/PSO/WOA revert to `eval_population_fitness()` (mclapply parallel); batch_fn approach
  was slower on real data.
- `positional_correction()` now accepts a `parallel` argument to control the number of
  workers used by population-based optimizers (defaults: GA=20, PSO=25, WOA=20, L-BFGS-B=1).

## Bug fixes

- Fixed `eval_population_fitness()` ignoring numeric `parallel=` values (was accidentally
  using `detectCores()-1` cores for any truthy value, causing inconsistent benchmarks).

# salpa 0.0.1.9.3

## Bug fixes

- Critical fix: `Rcpp::compileAttributes()` was not run after adding
  `ga_evolve_population_cpp`, causing the function to be missing from `RcppExports.R`
  and GA to silently return `best_x = 0` in all v0.0.1.9.x builds.

# salpa 0.0.1.8.0

## New features

- Replaced GA package dependency with a custom C++ GA implementation
  (`ga_evolve_population_cpp`: BLX-alpha crossover, tournament selection).
- Fork-safe parallel fitness evaluation via `terra::wrap()` / `terra::unwrap()`.

# salpa 0.0.1.7.0

## Bug fixes

- Fixed missing `parallel` parameter in `optimize_woa()`.
- Fixed WOA spiral formula to match Mirjalili (2016).
- Fixed adaptive `l` parameter (was fixed in [-1, 1]).

# salpa 0.0.1.6.1

## Bug fixes

- Fixed Euclidean distance formula (was computing sum instead of mean).
