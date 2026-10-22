#!/usr/bin/env bash
# Summary-statistics mode (median, min, max only). INFORM_PRIVATE_DIR (outside Git) must contain pre_mtg_survey_summary_gtz_20261019.csv.
set -euo pipefail
: "${INFORM_PRIVATE_DIR:?set INFORM_PRIVATE_DIR to a folder outside Git}"
cd "$(dirname "$0")/../.."
python3 05_response_surfaces/script/survey_prepare_summary.py
INFORM_OUT_SUBDIR=out_summary /usr/local/bin/Rscript 05_response_surfaces/script/export_survey_points.R
python3 05_response_surfaces/script/generate_summary_assets.py
python3 05_response_surfaces/script/acceptance_checks_summary.py
export PATH="$PATH:/Library/TeX/texbin"
O="$INFORM_PRIVATE_DIR/out_summary"; mkdir -p "$O/build" "$O/pdf"
( cd 05_response_surfaces && TEXINPUTS=".:$O:" latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir="$O/build" \
    proportional_funding_reductions_response_surfaces_summary.tex )
cp "$O/build/proportional_funding_reductions_response_surfaces_summary.pdf" "$O/pdf/"
echo "PDF: $O/pdf/proportional_funding_reductions_response_surfaces_summary.pdf"
