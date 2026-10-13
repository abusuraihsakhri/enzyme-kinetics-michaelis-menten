"""Compatibility application factory for the maintained FastAPI service.

Install with `pip install ".[api]"` before calling `create_app`.
"""


def create_app():
    """Return the canonical FastAPI application used by the container."""
    from agents.api import app
    return app
