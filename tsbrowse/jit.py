import functools
import logging
import os

import numba

logger = logging.getLogger(__name__)

_DISABLE_NUMBA = os.environ.get("TSBROWSE_DISABLE_NUMBA", "0")

try:
    ENABLE_NUMBA = {"0": True, "1": False}[_DISABLE_NUMBA]
except KeyError as e:  # pragma: no cover
    raise KeyError(
        "Environment variable 'TSBROWSE_DISABLE_NUMBA' must be '0' or '1'"
    ) from e

# We will mostly be using disable numba for debugging and running tests for
# coverage, so raise a loud warning in case this is being used accidentally.

if not ENABLE_NUMBA:
    logger.warning(
        "numba globally disabled for tsbrowse; performance will be drastically"
        " reduced."
    )


DEFAULT_NUMBA_ARGS = {
    "nopython": True,
    "cache": True,
}


def numba_jit(**numba_kwargs):
    def _numba_jit(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            return func(*args, **kwargs)  # pragma: no cover

        if ENABLE_NUMBA:  # pragma: no cover
            combined_kwargs = {**DEFAULT_NUMBA_ARGS, **numba_kwargs}
            return numba.jit(**combined_kwargs)(func)
        else:
            return func

    return _numba_jit


def numba_jitclass(spec):
    def _numba_jitclass(cls):
        if ENABLE_NUMBA:  # pragma: no cover
            return numba.experimental.jitclass(spec)(cls)
        else:
            return cls

    return _numba_jitclass
