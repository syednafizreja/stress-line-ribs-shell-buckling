import os

from cecore.env import setup_ccx


def _call():
    try:
        setup_ccx()
    except FileNotFoundError:   # CI has no CalculiX; the PATH logic runs before the lookup
        pass


def test_setup_ccx_is_idempotent():
    _call()
    n = len(os.environ["PATH"])
    for _ in range(50):
        _call()
    assert len(os.environ["PATH"]) == n
