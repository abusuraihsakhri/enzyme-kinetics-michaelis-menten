# Enzyme Kinetics Michaelis Menten

> **Domain:** Computational Biology & AI Drug Discovery  
> **Reference Guidelines & Standards:** `wwPDB, IUPAC & CLSI Computational Guidelines`

<div align="center">

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12-3776AB.svg?logo=python&logoColor=white)
![FastAPI](https://img.shields.io/badge/FastAPI-0.111-009688.svg?logo=fastapi&logoColor=white)
![Audit Trail](https://img.shields.io/badge/Audit-HMAC--SHA256_Tamper--Evident-brightgreen.svg)
![Zero-PHI Guard](https://img.shields.io/badge/Guard-Zero--PHI_Outbound-blue.svg)
![Docker](https://img.shields.io/badge/Docker-Ready-2496ED.svg?logo=docker&logoColor=white)

</div>

---

## 📖 What It Does

**Enzyme Kinetics Michaelis Menten** is an advanced analytical and computational platform implementing:

- **Non-linear Michaelis-Menten Vmax/Km fitter** with Lineweaver-Burk, Eadie-Hofstee, and Hanes-Woolf linearization methods
- **Inhibitor classification** (competitive, uncompetitive, non-competitive, mixed)
- **IC50 calculation** using the Cheng-Prusoff equation
- **Substrate inhibition analysis** with optimal concentration determination
- **Enterprise security features** including Zero-PHI outbound guard and HMAC-SHA256 audit trail

---

## ⚙️ Key Capabilities & Algorithmic Modules

- **Deterministic Calculation Engine**: Strict compliance with standard reference formulations and thresholds.
- **Risk & Urgency Classification**: Multi-tier categorization with automated clinical/operational action recommendations.
- **Validation & Guardrails**: Rigorous input bounds checking and anomaly detection.
- **Multi-Agent Supervisor**: Distributed worker architecture for quality control, safety escalation, and protocol conformance.

---

## 💻 Installation

```bash
# Clone the repository
git clone https://github.com/abusuraihsakhri/enzyme-kinetics-michaelis-menten.git
cd enzyme-kinetics-michaelis-menten

# Install dependencies
pip install -e .

# Or install with optional dependencies
pip install fastapi uvicorn pydantic pytest
```

---

## 💻 CLI Quickstart & Usage

### 1. Fit Michaelis-Menten Parameters
```bash
python cli.py fit -i sample.csv
python cli.py fit -i data.json --enzyme-conc 1.0 --json
```

### 2. Linearization Analysis
```bash
python cli.py linearize -i sample.csv --json
```

### 3. Inhibition Analysis
```bash
python cli.py inhibition --vmax 100 --km 10 --type competitive --inhibitor 5 --ki 2
python cli.py inhibition --vmax 100 --km 10 --type mixed --inhibitor 5 --ki 2 --ki-prime 5
```

### 4. IC50 Calculation
```bash
python cli.py ic50 --vmax 100 --km 10 --substrate 5 --ki 2 --type competitive
```

### 5. Substrate Inhibition
```bash
python cli.py substrate-inhibition --vmax 100 --km 10 --ki-substrate 100
```

### 6. Supervisor Audit
```bash
python cli.py audit --task-id TASK-001 --primary-metric 25.5 --json
```

### 7. Chat Assistant
```bash
python cli.py chat "Explain enzyme kinetics"
```

### 8. Verify Audit Trail
```bash
python cli.py verify-audit
```

### Input Data Schema (for fit/linearize commands)

| Field | Description | Requirement |
|:------|:------------|:------------|
| `concentration` | Substrate concentration [S] | Required |
| `velocity` | Reaction velocity v | Required |

CSV format:
```csv
concentration,velocity
1.0,9.1
2.0,16.7
5.0,33.3
10.0,50.0
```

JSON format:
```json
{
  "concentrations": [1, 2, 5, 10, 20, 50, 100],
  "velocities": [9.1, 16.7, 33.3, 50.0, 66.7, 83.3, 90.9]
}
```

---

## 🛡️ Security & Enterprise Architecture

* **Zero-PHI Outbound Interceptor:** Active regex inspection blocking SSNs, MRNs, phone numbers, emails, and patient identifiers.
* **Tamper-Evident HMAC-SHA256 Audit Trail:** Chained, cryptographically signed logs for every evaluation and state transition.
* **Air-Gapped LLM Reasoning Adapter:** Agnostic integration for local Ollama instances (`llama3`, `mistral`), Claude 3.5 Sonnet, GPT-4o, and deterministic test mocks.
* **Active Learning Bayesian Calibration:** Dynamic tracker updating worker reliability weights and monitoring Brier calibration drift.
* **FastAPI & Prometheus Telemetry:** Exposes OpenAPI 3.1 REST endpoints and operational Prometheus metrics (`/metrics`).

### Security Configuration

Set the `AUDIT_SECRET_KEY` environment variable for production use:
```bash
export AUDIT_SECRET_KEY="your-secure-key-here"
```

---

## 🧪 Testing & Verification

Run the automated test suite:

```bash
pytest -v
```

Execute high-throughput batch simulation benchmarks:

```bash
python simulator.py 1000
```

---

## 🐳 Container Deployment

```bash
docker build -t enzyme-kinetics-michaelis-menten .
docker run -p 8000:8000 -e AUDIT_SECRET_KEY=your-key enzyme-kinetics-michaelis-menten
```

Or using Docker Compose:

```bash
docker-compose up -d
```

---

## 📁 Project Structure

```
enzyme-kinetics-michaelis-menten/
├── agents/                  # Multi-agent supervisor system
│   ├── api.py              # FastAPI REST endpoints
│   ├── base.py             # Security, PHI guard, audit trail
│   ├── models.py           # Pydantic data models
│   ├── supervisor.py       # Master orchestrator
│   ├── workers.py          # Specialized worker agents
│   ├── llm_factory.py      # LLM provider factory
│   ├── learning.py         # Bayesian calibration engine
│   ├── metrics.py          # Prometheus metrics
│   └── streamer.py         # WebSocket telemetry
├── enzyme_kinetics/         # Core enzyme kinetics engine
│   ├── engine.py           # Michaelis-Menten, linearization, inhibition
│   ├── cli.py              # Command-line interface
│   ├── models.py           # Data model exports
│   └── server.py           # FastAPI app factory
├── tests/                   # Test suite
├── web/                     # Operations console (HTML)
├── sample.csv               # Sample data
├── simulator.py             # High-throughput simulation
├── Dockerfile               # Container definition
└── docker-compose.yml       # Compose configuration
```
