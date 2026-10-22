# =============================================================================
# Step 2: query the surrogate at the points prepared from the respondent file.
# =============================================================================
# Reads  <private-dir>/<out>/points_to_query.csv (from survey_prepare.py or survey_prepare_summary.py)
# Writes <private-dir>/<out>/points_draws.csv    (per-draw IRR at every point); <out> = INFORM_OUT_SUBDIR (default 'out')
# The private directory must be OUTSIDE any Git repository; this script refuses
# to run otherwise. No value from the respondent file is hard-coded here.
# Rows: surrogate row 11 = post-intervention year 10 (final); row 6 = year 5.
# Points with a coordinate above the surrogate support (0.75) get NA, never a
# clipped value. Constant z block (set.seed(0)) => identical draws across points.
# Run: INFORM_PRIVATE_DIR=/path/outside/git /usr/local/bin/Rscript 05_response_surfaces/script/export_survey_points.R
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
outsub <- Sys.getenv("INFORM_OUT_SUBDIR", "out")   # "out" (respondent mode) or "out_summary"
points_file <- file.path(pdir, outsub, "points_to_query.csv")
if (!file.exists(points_file)) stop("run survey_prepare.py first")

setwd(file.path(REPO_DIR, "04_cost_mapping", "script"))
library(data.table)
library(hetGP)
load("../../03_intervention_scenario_surrogate/output/surrogate.Rdata")
source("../R/parameters.R")
source("../R/irr_common_random_numbers.R")
p <- cost_mapping_params()

pts <- fread(points_file)
ROWS <- c(final = 11L, year5 = 6L)
R_MAX <- 0.75
NPER <- p$n_samples_per_checkpoint
stopifnot(nrow(gp_incidence_fit[[1]]$Kmat) == 11L, all(pts$r_prep >= 0), all(pts$r_vs >= 0))
ok <- pts$r_prep <= R_MAX + 1e-12 & pts$r_vs <= R_MAX + 1e-12

res <- list()
for (w in names(ROWS)) {
  set.seed(0)
  m <- compute_scenario_draws_crn(cbind(-pts$r_vs[ok], -pts$r_prep[ok]), gp_incidence_fit,
                                  n_samples_per_checkpoint = NPER, common_random_numbers = TRUE,
                                  tick = ROWS[[w]])$irr
  nd <- ncol(m)
  full <- matrix(NA_real_, nrow(pts), nd); full[ok, ] <- m
  for (i in seq_len(nrow(pts)))
    res[[length(res) + 1]] <- data.table(pts[i], when = w, row_index = ROWS[[w]], draw = seq_len(nd),
                                         checkpoint = ceiling(seq_len(nd) / NPER), irr = full[i, ])
}
out <- rbindlist(res)
stopifnot(all(is.finite(out$irr[out$r_prep <= R_MAX & out$r_vs <= R_MAX])))
fwrite(out, file.path(pdir, outsub, "points_draws.csv"))
cat("wrote points_draws.csv:", nrow(out), "rows;", sum(!ok), "point(s) outside surrogate support -> NA\n")
