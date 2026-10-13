"""Regression tests for numerical fitting, validation, and audit-chain integrity."""
import math
import pytest

from enzyme_kinetics.engine import (
    fit_michaelis_menten, lineweaver_burk, michaelis_menten,
    competitive_inhibition, ic50_to_ki, calculate_ic50,
)
from agents.base import AuditTrail


def test_nonlinear_fit_minimizes_velocity_space_error():
    concentrations = [1, 2, 5, 10, 20, 50, 100]
    velocities = [
        michaelis_menten(s, 100.0, 10.0) + noise
        for s, noise in zip(concentrations, [-1, 1, -1, 1, 0, -1, 1])
    ]
    fitted = fit_michaelis_menten(concentrations, velocities)
    lb = lineweaver_burk(concentrations, velocities)
    nonlinear_sse = sum(
        (v - michaelis_menten(s, fitted.vmax, fitted.km)) ** 2
        for s, v in zip(concentrations, velocities)
    )
    linear_sse = sum(
        (v - michaelis_menten(s, lb["vmax"], lb["km"])) ** 2
        for s, v in zip(concentrations, velocities)
    )
    assert nonlinear_sse <= linear_sse + 1e-8
    assert fitted.vmax == pytest.approx(100, rel=0.06)
    assert fitted.km == pytest.approx(10, rel=0.15)
    assert 0.0 < fitted.r_squared <= 1.0


@pytest.mark.parametrize("x,y", [
    ([1, 2, 3], [1, 2]),
    ([1, 2, float("nan")], [1, 2, 3]),
    ([1, 2, -3], [1, 2, 3]),
    ([1, 1, 1], [1, 2, 3]),
    ([1, 2, 3], [0, 0, 0]),
])
def test_reject_invalid_fitting_inputs(x, y):
    with pytest.raises(ValueError):
        fit_michaelis_menten(x, y)


@pytest.mark.parametrize("s,vmax,km", [
    (-1, 10, 1), (1, -10, 1), (1, 10, 0), (1, 10, float("inf"))
])
def test_reject_nonphysical_rate_parameters(s, vmax, km):
    with pytest.raises(ValueError):
        michaelis_menten(s, vmax, km)


def test_reject_zero_inhibitor_constant():
    with pytest.raises(ValueError):
        competitive_inhibition(1, 10, 1, 1, 0)


def test_reject_ic50_without_substrate():
    with pytest.raises(ValueError):
        calculate_ic50(100, 10, 0, 2)


def test_reject_inverse_ic50_without_substrate():
    with pytest.raises(ValueError):
        ic50_to_ki(5, 0, 1)


def test_hmac_detects_data_tampering():
    audit = AuditTrail("test-secret")
    audit.log("tester", "scientist", "FIT", {"count": 7})
    assert audit.verify_integrity()
    audit.logs[0]["payload_hash"] = "tampered"
    assert not audit.verify_integrity()


def test_hmac_detects_modified_signature():
    audit = AuditTrail("test-secret")
    audit.log("tester", "scientist", "FIT", {"count": 7})
    audit.logs[0]["current_hash"] = "00" * 32
    assert not audit.verify_integrity()


def test_audit_snapshot_is_not_mutable_internal_state():
    audit = AuditTrail("test-secret")
    audit.log("tester", "scientist", "FIT", {"count": 7})
    snapshot = audit.get_trail()
    snapshot[0]["actor"] = "intruder"
    assert audit.verify_integrity()
    assert audit.get_trail()[0]["actor"] == "tester"


def test_legacy_app_factory_is_usable():
    from enzyme_kinetics.server import create_app
    app = create_app()
    assert any(route.path == "/health" for route in app.routes)
    assert any(route.path == "/api/audit" for route in app.routes)
