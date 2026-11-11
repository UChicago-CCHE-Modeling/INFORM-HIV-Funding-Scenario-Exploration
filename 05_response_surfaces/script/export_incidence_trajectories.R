# =============================================================================
# Step 2b: absolute incidence trajectories for the incidence tables and ribbon plot.
# =============================================================================
# Reads  <private-dir>/out_summary/points_to_query.csv (set == "main": the elicited
#        PrEP component, ART component, both, marginal minima and maxima)
# Writes <private-dir>/out_summary/trajectory_summary.csv: one row per point x year
#        (0-10) with mean and 2.5/97.5% quantiles of absolute incidence (per 100 PY)
#        and of the IRR vs the (0,0) baseline. Point "baseline" is the surrogate at
#        (0,0). Reference points are round-number reductions (no survey input).
# Same seed and common random numbers as export_survey_points.R, so year 10
# (tick 11) and year 5 (tick 6) reproduce the IRRs in points_draws.csv exactly.
# Run: INFORM_PRIVATE_DIR=/path/outside/git /usr/local/bin/Rscript 05_response_surfaces/script/export_incidence_trajectories.R
# =============================================================================

args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grepl("^--file=", args)][1])
SCRIPT_DIR <- dirname(normalizePath(script_file))
STAGE_DIR  <- dirname(SCRIPT_DIR)
REPO_DIR   <- dirname(STAGE_DIR)

find_git_root <- function(p) {
  p <- normalizePath(p, mustWork = FALSE)
  repeat {
    if (file.exists(file.path(p, ".git"))) return(p)
    parent <- dirname(p)
    if (parent == p) return(NULL)
    p <- parent
  }
}
pdir <- Sys.getenv("INFORM_PRIVATE_DIR")
if (!nzchar(pdir)) stop("set INFORM_PRIVATE_DIR to a folder outside Git")
pdir <- normalizePath(pdir, mustWork = TRUE)
if (!is.null(find_git_root(pdir))) stop("REFUSED: ", pdir, " is inside a Git repository")
points_file <- file.path(pdir, "out_summary", "points_to_query.csv")
if (!file.exists(points_file)) stop("run survey_prepare_summary.py first")

setwd(file.path(REPO_DIR, "04_cost_mapping", "script"))
library(data.table)
library(hetGP)
load("../../03_intervention_scenario_surrogate/output/surrogate.Rdata")
source("../R/parameters.R")
source("../R/irr_common_random_numbers.R")
p <- cost_mapping_params()

# Round-number reference reductions spanning the surrogate domain (survey-free).
ref <- rbindlist(lapply(c(0.10, 0.25, 0.40), function(x) data.table(
  point = paste0(c("ref_vs_", "ref_prep_", "ref_both_"), round(100 * x)),
  r_prep = c(0, x, x), r_vs = c(x, 0, x))))
elic <- fread(points_file)[set == "main", .(point = paste0("elic_", name), r_prep, r_vs)]
pts <- rbind(ref, elic)
stopifnot(all(pts$r_prep <= 0.75), all(pts$r_vs <= 0.75))

set.seed(0)
d <- compute_scenario_trajectory_draws_crn(cbind(-pts$r_vs, -pts$r_prep), gp_incidence_fit,
                                           n_samples_per_checkpoint = p$n_samples_per_checkpoint,
                                           common_random_numbers = TRUE)
years <- seq_len(nrow(d$incidence[[1]])) - 1L
stopifnot(length(years) == 11L)

summ <- function(M, nm) {
  out <- data.table(rowMeans(M), apply(M, 1, quantile, 0.025), apply(M, 1, quantile, 0.975))
  setnames(out, paste0(nm, c("_mean", "_lo", "_hi")))
  out
}
base <- d$incidence[[1]] - d$diff[[1]]   # surrogate (0,0) incidence draws, shared under CRN
rows <- list(data.table(point = "baseline", r_prep = 0, r_vs = 0, year = years,
                        summ(base, "inc"), irr_mean = 1, irr_lo = 1, irr_hi = 1))
for (k in seq_len(nrow(pts)))
  rows[[k + 1]] <- data.table(pts[k], year = years, summ(d$incidence[[k]], "inc"), summ(d$irr[[k]], "irr"))
out <- rbindlist(rows)
fwrite(out, file.path(pdir, "out_summary", "trajectory_summary.csv"))
cat("wrote trajectory_summary.csv:", nrow(out), "rows,", nrow(pts) + 1L, "points x", length(years), "years\n")
