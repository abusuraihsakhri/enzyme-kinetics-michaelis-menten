"""
Real tests for Enzyme Kinetics: Michaelis-Menten, linearization, inhibition, IC50.
"""
import sys
import math
import json
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from enzyme_kinetics.engine import (
    michaelis_menten,
    michaelis_menten_batch,
    lineweaver_burk,
    eadie_hofstee,
    hanes_woolf,
    calculate_kcat,
    calculate_catalytic_efficiency,
    competitive_inhibition,
    uncompetitive_inhibition,
    noncompetitive_inhibition,
    mixed_inhibition,
    substrate_inhibition,
    substrate_inhibition_optimal_conc,
    calculate_ic50,
    ic50_to_ki,
    fit_michaelis_menten,
    fit_with_kcat,
    full_kinetic_analysis,
    _linear_regression,
)


# --- Test Michaelis-Menten Equation ---

def test_mm_at_km():
    """At [S] = Km, v should be Vmax/2."""
    v = michaelis_menten(10.0, 100.0, 10.0)
    assert abs(v - 50.0) < 0.001


def test_mm_low_substrate():
    """At very low [S], v ≈ Vmax*[S]/Km (first-order)."""
    v = michaelis_menten(0.1, 100.0, 10.0)
    expected = 100.0 * 0.1 / (10.0 + 0.1)
    assert abs(v - expected) < 0.001


def test_mm_high_substrate():
    """At very high [S], v approaches Vmax (zero-order)."""
    v = michaelis_menten(10000.0, 100.0, 10.0)
    assert abs(v - 100.0) < 0.1


def test_mm_zero_substrate():
    """At [S] = 0, v should be 0."""
    v = michaelis_menten(0.0, 100.0, 10.0)
    assert v == 0.0


def test_mm_batch():
    """Batch calculation should match individual calls."""
    concs = [1.0, 5.0, 10.0, 50.0, 100.0]
    velocities = michaelis_menten_batch(concs, 100.0, 10.0)
    assert len(velocities) == 5
    for s, v in zip(concs, velocities):
        assert abs(v - michaelis_menten(s, 100.0, 10.0)) < 0.001


# --- Test Linearization Methods ---

def _generate_mm_data(vmax=100.0, km=10.0, noise=0.0):
    """Generate Michaelis-Menten data for testing."""
    import random
    random.seed(42)
    concs = [1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0]
    velocities = []
    for s in concs:
        v = michaelis_menten(s, vmax, km)
        if noise > 0:
            v += random.gauss(0, noise)
        velocities.append(max(v, 0.001))
    return concs, velocities


def test_lineweaver_burk_recovery():
    """Lineweaver-Burk should recover known Vmax and Km from noiseless data."""
    concs, vels = _generate_mm_data(vmax=100.0, km=10.0)
    result = lineweaver_burk(concs, vels)
    assert abs(result['vmax'] - 100.0) < 1.0
    assert abs(result['km'] - 10.0) < 1.0
    assert result['r_squared'] > 0.99


def test_eadie_hofstee_recovery():
    """Eadie-Hofstee should recover known Vmax and Km."""
    concs, vels = _generate_mm_data(vmax=100.0, km=10.0)
    result = eadie_hofstee(concs, vels)
    assert abs(result['vmax'] - 100.0) < 1.0
    assert abs(result['km'] - 10.0) < 1.0
    assert result['r_squared'] > 0.99


def test_hanes_woolf_recovery():
    """Hanes-Woolf should recover known Vmax and Km."""
    concs, vels = _generate_mm_data(vmax=100.0, km=10.0)
    result = hanes_woolf(concs, vels)
    assert abs(result['vmax'] - 100.0) < 1.0
    assert abs(result['km'] - 10.0) < 1.0
    assert result['r_squared'] > 0.99


def test_all_linearizations_consistent():
    """All three linearization methods should give similar results."""
    concs, vels = _generate_mm_data(vmax=50.0, km=5.0)
    lb = lineweaver_burk(concs, vels)
    eh = eadie_hofstee(concs, vels)
    hw = hanes_woolf(concs, vels)
    
    # All should be within 10% of true values
    for result in [lb, eh, hw]:
        assert abs(result['vmax'] - 50.0) / 50.0 < 0.1
        assert abs(result['km'] - 5.0) / 5.0 < 0.1


# --- Test kcat and Catalytic Efficiency ---

def test_kcat():
    """kcat = Vmax / [E]total."""
    kcat = calculate_kcat(100.0, 1.0)
    assert kcat == 100.0


def test_kcat_different_conc():
    """kcat should scale with enzyme concentration."""
    kcat1 = calculate_kcat(100.0, 1.0)
    kcat2 = calculate_kcat(100.0, 0.5)
    assert kcat2 == 2 * kcat1


def test_kcat_zero_conc():
    """Zero enzyme concentration should raise error."""
    try:
        calculate_kcat(100.0, 0.0)
        assert False, "Should have raised ValueError"
    except ValueError:
        pass


def test_catalytic_efficiency():
    """Catalytic efficiency = kcat / Km."""
    eff = calculate_catalytic_efficiency(100.0, 10.0)
    assert abs(eff - 10.0) < 0.001


def test_catalytic_efficiency_zero_km():
    """Zero Km should raise error."""
    try:
        calculate_catalytic_efficiency(100.0, 0.0)
        assert False, "Should have raised ValueError"
    except ValueError:
        pass


# --- Test Inhibition Models ---

def test_competitive_no_inhibitor():
    """With [I] = 0, competitive inhibition should equal normal MM."""
    v_normal = michaelis_menten(10.0, 100.0, 10.0)
    v_comp = competitive_inhibition(10.0, 100.0, 10.0, 0.0, 5.0)
    assert abs(v_normal - v_comp) < 0.001


def test_competitive_increases_km():
    """Competitive inhibition should increase apparent Km."""
    v_normal = michaelis_menten(10.0, 100.0, 10.0)
    v_comp = competitive_inhibition(10.0, 100.0, 10.0, 10.0, 5.0)
    # alpha = 1 + 10/5 = 3, so apparent Km = 30
    # v = 100*10/(30+10) = 25
    assert v_comp < v_normal
    assert abs(v_comp - 25.0) < 0.001


def test_competitive_vmax_unchanged():
    """At infinite [S], competitive inhibition should reach same Vmax."""
    v_normal = michaelis_menten(100000.0, 100.0, 10.0)
    v_comp = competitive_inhibition(100000.0, 100.0, 10.0, 10.0, 5.0)
    assert abs(v_normal - v_comp) < 0.1


def test_uncompetitive_decreases_both():
    """Uncompetitive inhibition decreases both Vmax and Km."""
    v_normal = michaelis_menten(10.0, 100.0, 10.0)
    v_uncomp = uncompetitive_inhibition(10.0, 100.0, 10.0, 10.0, 5.0)
    assert v_uncomp < v_normal


def test_noncompetitive_decreases_vmax():
    """Non-competitive inhibition decreases Vmax."""
    v_normal = michaelis_menten(100000.0, 100.0, 10.0)
    v_noncomp = noncompetitive_inhibition(100000.0, 100.0, 10.0, 10.0, 5.0)
    assert v_noncomp < v_normal


def test_mixed_inhibition():
    """Mixed inhibition with Ki != Ki' should behave differently from pure noncompetitive."""
    v_mixed = mixed_inhibition(10.0, 100.0, 10.0, 10.0, 5.0, 10.0)
    v_noncomp = noncompetitive_inhibition(10.0, 100.0, 10.0, 10.0, 5.0, 5.0)
    # Different Ki' values should give different results
    assert v_mixed != v_noncomp


# --- Test Substrate Inhibition ---

def test_substrate_inhibition_optimal():
    """Optimal [S] = sqrt(Km * Ki)."""
    s_opt = substrate_inhibition_optimal_conc(10.0, 100.0)
    assert abs(s_opt - math.sqrt(1000.0)) < 0.001


def test_substrate_inhibition_decreases_at_high_s():
    """At very high [S], velocity should decrease."""
    v_low = substrate_inhibition(10.0, 100.0, 10.0, 100.0)
    v_high = substrate_inhibition(500.0, 100.0, 10.0, 100.0)
    assert v_high < v_low


def test_substrate_inhibition_at_optimal():
    """Velocity at optimal [S] should be maximum."""
    km, ki = 10.0, 100.0
    s_opt = substrate_inhibition_optimal_conc(km, ki)
    v_opt = substrate_inhibition(s_opt, 100.0, km, ki)
    
    # Check nearby points are lower
    v_below = substrate_inhibition(s_opt * 0.5, 100.0, km, ki)
    v_above = substrate_inhibition(s_opt * 2.0, 100.0, km, ki)
    assert v_opt >= v_below
    assert v_opt >= v_above


# --- Test IC50 ---

def test_ic50_competitive():
    """IC50 = Ki * (1 + [S]/Km) for competitive inhibition."""
    ic50 = calculate_ic50(100.0, 10.0, 5.0, 2.0, 'competitive')
    expected = 2.0 * (1.0 + 5.0 / 10.0)
    assert abs(ic50 - expected) < 0.001


def test_ic50_uncompetitive():
    """IC50 = Ki * (1 + Km/[S]) for uncompetitive inhibition."""
    ic50 = calculate_ic50(100.0, 10.0, 5.0, 2.0, 'uncompetitive')
    expected = 2.0 * (1.0 + 10.0 / 5.0)
    assert abs(ic50 - expected) < 0.001


def test_ic50_noncompetitive():
    """IC50 = Ki for non-competitive inhibition."""
    ic50 = calculate_ic50(100.0, 10.0, 5.0, 2.0, 'noncompetitive')
    assert abs(ic50 - 2.0) < 0.001


def test_ic50_to_ki_roundtrip():
    """IC50 -> Ki -> IC50 should be consistent."""
    ki = 5.0
    s = 10.0
    km = 20.0
    ic50 = calculate_ic50(100.0, km, s, ki, 'competitive')
    ki_back = ic50_to_ki(ic50, s, km, 'competitive')
    assert abs(ki_back - ki) < 0.001


# --- Test Parameter Fitting ---

def test_fit_mm():
    """fit_michaelis_menten should recover parameters."""
    concs, vels = _generate_mm_data(vmax=80.0, km=15.0)
    params = fit_michaelis_menten(concs, vels)
    assert abs(params.vmax - 80.0) / 80.0 < 0.1
    assert abs(params.km - 15.0) / 15.0 < 0.1


def test_fit_with_kcat():
    """fit_with_kcat should compute kcat and efficiency."""
    concs, vels = _generate_mm_data(vmax=100.0, km=10.0)
    params = fit_with_kcat(concs, vels, enzyme_conc=1.0)
    assert params.kcat is not None
    assert params.catalytic_efficiency is not None
    assert abs(params.kcat - params.vmax) < 1.0


def test_full_kinetic_analysis():
    """full_kinetic_analysis should return all linearization results."""
    concs, vels = _generate_mm_data(vmax=100.0, km=10.0)
    result = full_kinetic_analysis(concs, vels, enzyme_conc=1.0)
    assert result.parameters.vmax > 0
    assert result.parameters.km > 0
    assert result.parameters.kcat is not None
    assert result.lineweaver_burk is not None
    assert result.eadie_hofstee is not None
    assert result.hanes_woolf is not None
    assert len(result.velocities) == len(concs)
    assert len(result.predicted_velocities) == len(concs)
    assert len(result.residuals) == len(concs)


# --- Test Linear Regression ---

def test_linear_regression_perfect():
    """Perfect linear data should give R² = 1."""
    x = [1, 2, 3, 4, 5]
    y = [2, 4, 6, 8, 10]
    slope, intercept, r_sq = _linear_regression(x, y)
    assert abs(slope - 2.0) < 0.001
    assert abs(intercept) < 0.001
    assert abs(r_sq - 1.0) < 0.001


def test_linear_regression_constant():
    """Constant y should give slope ≈ 0."""
    x = [1, 2, 3, 4, 5]
    y = [5, 5, 5, 5, 5]
    slope, intercept, r_sq = _linear_regression(x, y)
    assert abs(slope) < 0.001
    assert abs(intercept - 5.0) < 0.001


# --- Test CLI ---

def test_cli_fit(tmp_path):
    """Test CLI fit command."""
    data = {'concentrations': [1, 2, 5, 10, 20, 50, 100],
            'velocities': [9.1, 16.7, 33.3, 50.0, 66.7, 83.3, 90.9]}
    f = tmp_path / 'data.json'
    f.write_text(json.dumps(data))
    
    from enzyme_kinetics.cli import main
    assert main(['fit', '-i', str(f)]) == 0


def test_cli_fit_json(tmp_path):
    """Test CLI fit with JSON output."""
    data = {'concentrations': [1, 2, 5, 10, 20, 50, 100],
            'velocities': [9.1, 16.7, 33.3, 50.0, 66.7, 83.3, 90.9]}
    f = tmp_path / 'data.json'
    f.write_text(json.dumps(data))
    
    from enzyme_kinetics.cli import main
    assert main(['fit', '-i', str(f), '--json']) == 0


def test_cli_ic50():
    """Test CLI IC50 command."""
    from enzyme_kinetics.cli import main
    assert main(['ic50', '--vmax', '100', '--km', '10', '--substrate', '5', '--ki', '2']) == 0


def test_cli_substrate_inhibition():
    """Test CLI substrate inhibition command."""
    from enzyme_kinetics.cli import main
    assert main(['substrate-inhibition', '--vmax', '100', '--km', '10', '--ki-substrate', '100']) == 0
