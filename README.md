# Enzyme Kinetics — Michaelis–Menten

A Python enzyme-kinetics library and CLI with a dependency-free browser workbench for initial-rate analysis.

## Browser workbench

The static app in [web/](web/) fits **Vmax** and **Km** by unweighted nonlinear least squares on the original velocity observations. It displays measured points and the fitted curve, R², predicted rates and residuals. It also calculates velocity under competitive, uncompetitive, pure noncompetitive and mixed inhibition, plus ideal-model IC50 relationships where defined. Results can be exported as CSV.

Enter a CSV file with `concentration,velocity` columns or use the bundled example. At least three distinct positive substrate concentrations and some nonzero measured rates are required. Keep concentration and velocity units consistent. An optional enzyme concentration supports `kcat = Vmax/[E]` and `kcat/Km`.

To run locally, serve the static directory (ES modules normally require HTTP rather than `file://`):

```bash
python -m http.server 8001 --directory web
# Open http://localhost:8001
```

A [GitHub Pages deployment workflow](.github/workflows/pages.yml) publishes `web/` from `master`. The repository must have **Settings → Pages → Source: GitHub Actions** enabled. The live link should only be published here once deployment and the application have been verified.

## Python installation and CLI

Python 3.10+ is required.

```bash
git clone https://github.com/abusuraihsakhri/enzyme-kinetics-michaelis-menten.git
cd enzyme-kinetics-michaelis-menten
python -m pip install -e .
python cli.py fit -i sample_kinetics.csv --json
python cli.py linearize -i sample_kinetics.csv
python cli.py inhibition --vmax 100 --km 10 --type competitive --inhibitor 5 --ki 2
python cli.py ic50 --vmax 100 --km 10 --substrate 5 --ki 2
python cli.py substrate-inhibition --vmax 100 --km 10 --ki-substrate 100
```

The original `sample.csv` contains *audit-task inputs*, not concentration/velocity pairs; `sample_kinetics.csv` is the appropriate example for fitting and linearizations. JSON input with `concentrations` and `velocities` arrays is also supported.

The core module implements Michaelis–Menten velocities, three reciprocal/linearized diagnostic methods (Lineweaver–Burk, Eadie–Hofstee, Hanes–Woolf), inhibition models, turnover number, catalytic efficiency, and substrate inhibition. Linearizations are **not used as the final nonlinear fit** because reciprocal transformations alter measurement-error weighting. No confidence intervals or automatic inhibition-mechanism classification are currently provided by the numerical engine.

## Optional audit API and container

The repository also retains a separate rule-based supervisor for synthetic audit tasks. Its fixed thresholds and enrichment heuristics are **not validated clinical decision rules or experimental enzyme-kinetics classifiers**.

```bash
python -m pip install -e ".[api]"
uvicorn agents.api:app --host 127.0.0.1 --port 8000
# API documentation: http://127.0.0.1:8000/docs
```

Endpoints include `/health`, `/metrics` (JSON), `/api/audit`, `/api/chat` (deterministic mock), and `/api/audit/logs`. The compatibility factory `enzyme_kinetics.server.create_app()` returns the same application.

For Docker Compose, set a strong, private `AUDIT_SECRET_KEY` in the environment before running `docker compose up --build`. No signing key is committed to the current Compose configuration. Prior deployed keys must be rotated because earlier repository history contained an example signing key.

## Privacy and limits

The browser workbench performs calculations entirely on the client; file import and CSV export do not send data to a server. Optional FastAPI endpoints receive submitted audit payloads and keep an in-process audit trail. The identifier-screening regexes are only basic safeguards and **do not establish de-identification, HIPAA compliance, or clinical suitability**. Do not submit patient-identifiable information.

The nonlinear fit assumes initial-rate measurements, an ordinary Michaelis–Menten mechanism and unweighted residuals. Confirm model adequacy and experimental assumptions before interpreting parameters. IC50 formulae represent ideal inhibition mechanisms and are not substitutes for fitting measured inhibitor dose–response curves.

## Development and tests

```bash
python -m pip install -e ".[api,dev]"
python -m pytest -q
python -m compileall -q agents enzyme_kinetics
node --test web/kinetics.test.mjs
```

The GitHub Actions workflow tests Python 3.10–3.13, package installation, numerical regressions, CLI commands and browser calculation modules. The static workbench targets modern browsers supporting ECMAScript modules, SVG, File and Blob APIs.

## License

[MIT](LICENSE).
