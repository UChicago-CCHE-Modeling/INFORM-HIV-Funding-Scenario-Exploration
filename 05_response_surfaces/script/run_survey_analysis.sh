#!/usr/bin/env bash
# Run the full respondent-level analysis. INFORM_PRIVATE_DIR must point to a folder OUTSIDE Git that holds
# config.json (see respondent_tools/config_template.json) and the respondent CSV it names.
set -euo pipefail
: "${INFORM_PRIVATE_DIR:?set INFORM_PRIVATE_DIR to a folder outside Git}"
cd "$(dirname "$0")/../.."
python3 05_response_surfaces/script/survey_prepare.py
/usr/local/bin/Rscript 05_response_surfaces/script/export_survey_points.R
python3 05_response_surfaces/script/generate_report_assets.py
python3 05_response_surfaces/script/acceptance_checks_survey.py
export PATH="$PATH:/Library/TeX/texbin"
mkdir -p "$INFORM_PRIVATE_DIR/out/build" "$INFORM_PRIVATE_DIR/out/pdf"
( cd 05_response_surfaces && TEXINPUTS=".:$INFORM_PRIVATE_DIR/out:" latexmk -pdf -interaction=nonstopmode -halt-on-error \
    -outdir="$INFORM_PRIVATE_DIR/out/build" proportional_funding_reductions_response_surfaces.tex )
cp "$INFORM_PRIVATE_DIR/out/build/proportional_funding_reductions_response_surfaces.pdf" "$INFORM_PRIVATE_DIR/out/pdf/"
echo "PDF: $INFORM_PRIVATE_DIR/out/pdf/proportional_funding_reductions_response_surfaces.pdf"
