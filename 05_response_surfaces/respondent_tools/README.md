# respondent_tools: local respondent-validation utility

Generic code only. It contains no respondent values and must be run on data kept **outside any Git repository**.

- `respondent_core.py`: privacy guard (`assert_outside_git`), config/CSV loading, per-respondent r = q(1-gamma), alternative r = q, componentwise medians, P1/P2/P3, and a read-only surface from the tracked survey-free grid.
- `respondent_validation.py`: writes one anonymous PDF per respondent. Each shows only that respondent's q, gamma, constructed r, an r-space plot, their P1/P2/P3 outcomes, and the question whether their reduction already accounted for alternative access. No network access (sockets are disabled); it refuses to read or write inside a Git repository.
- `config_template.json`: template with SYNTHETIC column names. Copy it next to your CSV (outside Git) and edit `input_csv`, `id_column`, the four column names, and `units` ("percent" or "proportion"). Set `data_status` to "real" for real data.
- `make_synthetic_data.py`: writes synthetic data and a matching config for testing. `test_respondent_tools.py`: tests (synthetic data only).

```sh
python3 respondent_tools/respondent_validation.py --config /path/outside/git/config.json --out-dir /path/outside/git/sheets [--anonymize]
```
`--anonymize` replaces ids with R01.. and writes `id_key.csv` (keep it private). Filenames use the id column, so supply anonymous ids.
