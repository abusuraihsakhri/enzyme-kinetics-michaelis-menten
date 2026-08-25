# Enzyme Kinetics: Michaelis-Menten

Real implementations of enzyme kinetics analysis: Michaelis-Menten equation, linearization methods (Lineweaver-Burk, Eadie-Hofstee, Hanes-Woolf), inhibition models, IC50 calculation, and substrate inhibition.

## Features

- **Michaelis-Menten**: v = Vmax × [S] / (Km + [S])
- **Lineweaver-Burk**: 1/v = (Km/Vmax)×(1/[S]) + 1/Vmax
- **Eadie-Hofstee**: v = Vmax - Km×(v/[S])
- **Hanes-Woolf**: [S]/v = [S]/Vmax + Km/Vmax
- **kcat & Catalytic Efficiency**: kcat = Vmax/[E], efficiency = kcat/Km
- **Inhibition Models**: Competitive, uncompetitive, non-competitive, mixed
- **IC50**: Cheng-Prusoff equation for all inhibition types
- **Substrate Inhibition**: v = Vmax×[S] / (Km + [S]×(1 + [S]/Ki))

## Quick Start

```bash
# Fit Michaelis-Menten parameters from data
python cli.py fit -i data.json

# All linearization methods
python cli.py linearize -i data.json

# Inhibition analysis
python cli.py inhibition --vmax 100 --km 10 --type competitive --inhibitor 5 --ki 2

# IC50 calculation
python cli.py ic50 --vmax 100 --km 10 --substrate 5 --ki 2 --type competitive

# Substrate inhibition
python cli.py substrate-inhibition --vmax 100 --km 10 --ki-substrate 100
```

## Input Format

JSON file:
```json
{"concentrations": [1, 2, 5, 10, 20, 50, 100],
 "velocities": [9.1, 16.7, 33.3, 50.0, 66.7, 83.3, 90.9]}
```

Or CSV with `concentration` and `velocity` columns.

## Python API

```python
from enzyme_kinetics import (
    michaelis_menten, fit_michaelis_menten, full_kinetic_analysis,
    competitive_inhibition, calculate_ic50, calculate_kcat,
)

# Calculate velocity
v = michaelis_menten(s=10.0, vmax=100.0, km=10.0)

# Fit from data
params = fit_michaelis_menten(concentrations, velocities)
print(f"Vmax={params.vmax:.2f}, Km={params.km:.2f}")

# IC50
ic50 = calculate_ic50(vmax=100, km=10, s=5, ki=2, inhibition_type='competitive')
```

## Inhibition Models

| Type | Equation | Effect |
|------|----------|--------|
| Competitive | v = Vmax×[S] / (α×Km + [S]) | ↑ apparent Km |
| Uncompetitive | v = Vmax×[S] / (Km + α'×[S]) | ↓ both Km and Vmax |
| Non-competitive | v = Vmax×[S] / (α×Km + α'×[S]) | ↓ Vmax only |
| Mixed | v = Vmax×[S] / (α×Km + α'×[S]) | ↓ Vmax, changes Km |

## License

MIT License.
