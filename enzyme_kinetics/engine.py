"""
Enzyme Kinetics Engine: Michaelis-Menten, Lineweaver-Burk, Eadie-Hofstee,
Hanes-Woolf, inhibition models, IC50, and substrate inhibition.
"""
import math
from typing import Dict, Any, List, Optional, Tuple
from dataclasses import dataclass, field


# --- Data Models ---

@dataclass
class KineticParameters:
    """Result of kinetic parameter fitting."""
    vmax: float
    km: float
    kcat: Optional[float] = None
    catalytic_efficiency: Optional[float] = None
    enzyme_conc: Optional[float] = None
    r_squared: Optional[float] = None


@dataclass
class InhibitionResult:
    """Result of inhibition analysis."""
    inhibition_type: str  # competitive, uncompetitive, noncompetitive, mixed
    ki: float
    alpha: float = 1.0
    alpha_prime: float = 1.0
    vmax_apparent: Optional[float] = None
    km_apparent: Optional[float] = None


@dataclass
class SubstrateInhibitionResult:
    """Result of substrate inhibition analysis."""
    vmax: float
    km: float
    ki_substrate: float
    optimal_concentration: float
    max_velocity_at_optimal: float


@dataclass
class KineticFitResult:
    """Complete kinetic analysis result."""
    parameters: KineticParameters
    velocities: List[float]
    concentrations: List[float]
    predicted_velocities: List[float]
    residuals: List[float]
    lineweaver_burk: Optional[Dict[str, Any]] = None
    eadie_hofstee: Optional[Dict[str, Any]] = None
    hanes_woolf: Optional[Dict[str, Any]] = None


# --- Core Michaelis-Menten ---

def michaelis_menten(s: float, vmax: float, km: float) -> float:
    """Michaelis-Menten equation: v = Vmax * [S] / (Km + [S])
    
    Args:
        s: Substrate concentration
        vmax: Maximum velocity
        km: Michaelis constant
    Returns:
        Reaction velocity
    """
    if km + s == 0:
        return 0.0
    return vmax * s / (km + s)


def michaelis_menten_batch(concentrations: List[float], vmax: float, km: float) -> List[float]:
    """Calculate velocities for multiple substrate concentrations."""
    return [michaelis_menten(s, vmax, km) for s in concentrations]


# --- Linearization Methods ---

def lineweaver_burk(concentrations: List[float], velocities: List[float]) -> Dict[str, Any]:
    """Lineweaver-Burk (double reciprocal) plot: 1/v = (Km/Vmax)*(1/[S]) + 1/Vmax
    
    Linear regression on 1/[S] vs 1/v:
    - slope = Km/Vmax
    - intercept = 1/Vmax
    
    Returns dict with Vmax, Km, slope, intercept, r_squared.
    """
    inv_s = []
    inv_v = []
    for s, v in zip(concentrations, velocities):
        if s > 0 and v > 0:
            inv_s.append(1.0 / s)
            inv_v.append(1.0 / v)
    
    if len(inv_s) < 2:
        raise ValueError("Need at least 2 valid data points for Lineweaver-Burk")
    
    slope, intercept, r_sq = _linear_regression(inv_s, inv_v)
    
    if intercept == 0:
        raise ValueError("Lineweaver-Burk intercept is zero; cannot compute Vmax")
    
    vmax = 1.0 / intercept
    km = slope * vmax
    
    return {
        'vmax': vmax,
        'km': km,
        'slope': slope,
        'intercept': intercept,
        'r_squared': r_sq,
        'inv_s': inv_s,
        'inv_v': inv_v,
    }


def eadie_hofstee(concentrations: List[float], velocities: List[float]) -> Dict[str, Any]:
    """Eadie-Hofstee plot: v = Vmax - Km*(v/[S])
    
    Linear regression on v/[S] vs v:
    - slope = -Km
    - intercept = Vmax
    
    Returns dict with Vmax, Km, slope, intercept, r_squared.
    """
    v_over_s = []
    v_vals = []
    for s, v in zip(concentrations, velocities):
        if s > 0:
            v_over_s.append(v / s)
            v_vals.append(v)
    
    if len(v_over_s) < 2:
        raise ValueError("Need at least 2 valid data points for Eadie-Hofstee")
    
    # Regression: v = Vmax - Km * (v/[S])
    # y = v_vals, x = v_over_s
    slope, intercept, r_sq = _linear_regression(v_over_s, v_vals)
    
    km = -slope
    vmax = intercept
    
    return {
        'vmax': vmax,
        'km': km,
        'slope': slope,
        'intercept': intercept,
        'r_squared': r_sq,
        'v_over_s': v_over_s,
        'v': v_vals,
    }


def hanes_woolf(concentrations: List[float], velocities: List[float]) -> Dict[str, Any]:
    """Hanes-Woolf plot: [S]/v = [S]/Vmax + Km/Vmax
    
    Linear regression on [S] vs [S]/v:
    - slope = 1/Vmax
    - intercept = Km/Vmax
    
    Returns dict with Vmax, Km, slope, intercept, r_squared.
    """
    s_vals = []
    s_over_v = []
    for s, v in zip(concentrations, velocities):
        if v > 0:
            s_vals.append(s)
            s_over_v.append(s / v)
    
    if len(s_vals) < 2:
        raise ValueError("Need at least 2 valid data points for Hanes-Woolf")
    
    # Regression: [S]/v = (1/Vmax)*[S] + Km/Vmax
    slope, intercept, r_sq = _linear_regression(s_vals, s_over_v)
    
    if slope == 0:
        raise ValueError("Hanes-Woolf slope is zero; cannot compute Vmax")
    
    vmax = 1.0 / slope
    km = intercept * vmax
    
    return {
        'vmax': vmax,
        'km': km,
        'slope': slope,
        'intercept': intercept,
        'r_squared': r_sq,
        's': s_vals,
        's_over_v': s_over_v,
    }


# --- kcat and Catalytic Efficiency ---

def calculate_kcat(vmax: float, enzyme_conc: float) -> float:
    """Calculate turnover number: kcat = Vmax / [E]total
    
    Args:
        vmax: Maximum velocity (same units as enzyme_conc/time)
        enzyme_conc: Total enzyme concentration
    Returns:
        kcat (turnover number, 1/time)
    """
    if enzyme_conc <= 0:
        raise ValueError("Enzyme concentration must be positive")
    return vmax / enzyme_conc


def calculate_catalytic_efficiency(kcat: float, km: float) -> float:
    """Calculate catalytic efficiency: kcat/Km
    
    Args:
        kcat: Turnover number
        km: Michaelis constant
    Returns:
        Catalytic efficiency (1/(concentration*time))
    """
    if km <= 0:
        raise ValueError("Km must be positive")
    return kcat / km


# --- Inhibition Models ---

def competitive_inhibition(s: float, vmax: float, km: float,
                            inhibitor: float, ki: float) -> float:
    """Competitive inhibition: v = Vmax*[S] / (alpha*Km + [S])
    where alpha = 1 + [I]/Ki
    
    Inhibitor competes with substrate for active site.
    Increases apparent Km, Vmax unchanged.
    """
    alpha = 1.0 + inhibitor / ki
    return vmax * s / (alpha * km + s)


def uncompetitive_inhibition(s: float, vmax: float, km: float,
                              inhibitor: float, ki: float) -> float:
    """Uncompetitive inhibition: v = Vmax*[S] / (Km + alpha'*[S])
    where alpha' = 1 + [I]/Ki
    
    Inhibitor binds only to ES complex.
    Decreases both apparent Km and Vmax.
    """
    alpha_prime = 1.0 + inhibitor / ki
    return vmax * s / (km + alpha_prime * s)


def noncompetitive_inhibition(s: float, vmax: float, km: float,
                               inhibitor: float, ki: float,
                               ki_prime: Optional[float] = None) -> float:
    """Non-competitive inhibition: v = Vmax*[S] / (alpha*Km + alpha'*[S])
    where alpha = 1 + [I]/Ki, alpha' = 1 + [I]/Ki'
    
    Pure non-competitive: Ki = Ki' (alpha = alpha').
    Decreases Vmax, Km unchanged.
    """
    alpha = 1.0 + inhibitor / ki
    alpha_prime = 1.0 + inhibitor / (ki_prime if ki_prime else ki)
    return vmax * s / (alpha * km + alpha_prime * s)


def mixed_inhibition(s: float, vmax: float, km: float,
                      inhibitor: float, ki: float, ki_prime: float) -> float:
    """Mixed inhibition: v = Vmax*[S] / (alpha*Km + alpha'*[S])
    where alpha = 1 + [I]/Ki, alpha' = 1 + [I]/Ki'
    
    Inhibitor binds both E and ES with different affinities.
    Decreases Vmax and changes apparent Km.
    """
    alpha = 1.0 + inhibitor / ki
    alpha_prime = 1.0 + inhibitor / ki_prime
    return vmax * s / (alpha * km + alpha_prime * s)


def substrate_inhibition(s: float, vmax: float, km: float, ki_substrate: float) -> float:
    """Substrate inhibition: v = Vmax*[S] / (Km + [S]*(1 + [S]/Ki))
    
    At high substrate concentrations, velocity decreases.
    Optimal [S] = sqrt(Km * Ki)
    """
    denominator = km + s * (1.0 + s / ki_substrate)
    if denominator == 0:
        return 0.0
    return vmax * s / denominator


def substrate_inhibition_optimal_conc(km: float, ki_substrate: float) -> float:
    """Calculate optimal substrate concentration for substrate inhibition model.
    S_opt = sqrt(Km * Ki)
    """
    return math.sqrt(km * ki_substrate)


# --- IC50 Calculation ---

def calculate_ic50(vmax: float, km: float, s: float,
                    ki: float, inhibition_type: str = 'competitive') -> float:
    """Calculate IC50 (inhibitor concentration causing 50% inhibition).
    
    Uses the Cheng-Prusoff equation:
    - Competitive: IC50 = Ki * (1 + [S]/Km)
    - Uncompetitive: IC50 = Ki * (1 + Km/[S])
    - Non-competitive: IC50 = Ki (when Ki = Ki')
    
    Args:
        vmax: Maximum velocity
        km: Michaelis constant
        s: Substrate concentration at which IC50 is measured
        ki: Inhibition constant
        inhibition_type: Type of inhibition
    Returns:
        IC50 value
    """
    if inhibition_type == 'competitive':
        return ki * (1.0 + s / km)
    elif inhibition_type == 'uncompetitive':
        if s == 0:
            return float('inf')
        return ki * (1.0 + km / s)
    elif inhibition_type in ('noncompetitive', 'non-competitive'):
        return ki
    else:
        raise ValueError(f"Unknown inhibition type: {inhibition_type}")


def ic50_to_ki(ic50: float, s: float, km: float,
                inhibition_type: str = 'competitive') -> float:
    """Convert IC50 to Ki using the Cheng-Prusoff equation (inverse).
    
    Args:
        ic50: Measured IC50
        s: Substrate concentration
        km: Michaelis constant
        inhibition_type: Type of inhibition
    Returns:
        Ki value
    """
    if inhibition_type == 'competitive':
        return ic50 / (1.0 + s / km)
    elif inhibition_type == 'uncompetitive':
        if km == 0:
            return ic50
        return ic50 / (1.0 + km / s)
    elif inhibition_type in ('noncompetitive', 'non-competitive'):
        return ic50
    else:
        raise ValueError(f"Unknown inhibition type: {inhibition_type}")


# --- Parameter Fitting (Least Squares) ---

def fit_michaelis_menten(concentrations: List[float],
                          velocities: List[float]) -> KineticParameters:
    """Estimate Vmax and Km from experimental data using Lineweaver-Burk linearization.
    
    For production use, non-linear regression (Levenberg-Marquardt) would be preferred,
    but this provides a reasonable first estimate.
    """
    lb = lineweaver_burk(concentrations, velocities)
    return KineticParameters(
        vmax=lb['vmax'],
        km=lb['km'],
        r_squared=lb['r_squared'],
    )


def fit_with_kcat(concentrations: List[float], velocities: List[float],
                   enzyme_conc: float) -> KineticParameters:
    """Fit Michaelis-Menten and calculate kcat and catalytic efficiency."""
    params = fit_michaelis_menten(concentrations, velocities)
    params.enzyme_conc = enzyme_conc
    params.kcat = calculate_kcat(params.vmax, enzyme_conc)
    params.catalytic_efficiency = calculate_catalytic_efficiency(params.kcat, params.km)
    return params


def full_kinetic_analysis(concentrations: List[float],
                           velocities: List[float],
                           enzyme_conc: Optional[float] = None) -> KineticFitResult:
    """Run complete kinetic analysis with all linearization methods.
    
    Returns fitted parameters and all linearization results.
    """
    params = fit_michaelis_menten(concentrations, velocities)
    
    if enzyme_conc is not None:
        params.enzyme_conc = enzyme_conc
        params.kcat = calculate_kcat(params.vmax, enzyme_conc)
        params.catalytic_efficiency = calculate_catalytic_efficiency(params.kcat, params.km)
    
    predicted = michaelis_menten_batch(concentrations, params.vmax, params.km)
    residuals = [v - p for v, p in zip(velocities, predicted)]
    
    lb = lineweaver_burk(concentrations, velocities)
    eh = eadie_hofstee(concentrations, velocities)
    hw = hanes_woolf(concentrations, velocities)
    
    return KineticFitResult(
        parameters=params,
        velocities=list(velocities),
        concentrations=list(concentrations),
        predicted_velocities=predicted,
        residuals=residuals,
        lineweaver_burk=lb,
        eadie_hofstee=eh,
        hanes_woolf=hw,
    )


# --- Linear Regression Utility ---

def _linear_regression(x: List[float], y: List[float]) -> Tuple[float, float, float]:
    """Simple linear regression: y = slope*x + intercept.
    
    Returns (slope, intercept, r_squared).
    """
    n = len(x)
    if n < 2:
        raise ValueError("Need at least 2 points for regression")
    
    sum_x = sum(x)
    sum_y = sum(y)
    sum_xy = sum(xi * yi for xi, yi in zip(x, y))
    sum_x2 = sum(xi ** 2 for xi in x)
    sum_y2 = sum(yi ** 2 for yi in y)
    
    denom = n * sum_x2 - sum_x ** 2
    if abs(denom) < 1e-15:
        raise ValueError("Degenerate data: all x values are identical")
    
    slope = (n * sum_xy - sum_x * sum_y) / denom
    intercept = (sum_y - slope * sum_x) / n
    
    # R-squared
    ss_res = sum((yi - (slope * xi + intercept)) ** 2 for xi, yi in zip(x, y))
    mean_y = sum_y / n
    ss_tot = sum((yi - mean_y) ** 2 for yi in y)
    
    if ss_tot < 1e-15:
        r_squared = 1.0  # all y values identical
    else:
        r_squared = 1.0 - ss_res / ss_tot
    
    return slope, intercept, r_squared
