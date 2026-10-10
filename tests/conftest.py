"""Test session configuration.

sbmlutils resolves the annotations of a model with the registry of
identifiers.org, which pymetadata downloads into its cache on first use. The
tests run in parallel, so on a machine without the cache every worker which
builds an annotated model, and every snippet of the documentation which does,
would download the registry and replace the cache file at the same time; on
Windows a replace of a file which another process reads fails
(matthiaskoenig/pymetadata#102). The registry is therefore loaded once in the
controller of the session, before the workers start, and the workers only read
the cache.

The SBML model of the toy model package (`tests/pkpdsim_toy`) is created in the
controller of the session as well, before the workers start; the tests only read
it.
"""

from contextlib import suppress

import pytest
from pymetadata.webservices.registry import get_registry
from pymetadata.webservices.webservice import WebserviceError

from pkpdsim_toy import MODEL_BASE_PATH
from pkpdsim_toy.factory import create_toy_models


def pytest_configure(config: pytest.Config) -> None:
    """Load the registry of identifiers.org and create the toy model."""
    if hasattr(config, "workerinput"):
        # a worker of pytest-xdist, the controller has prepared the session
        return
    # offline without a cache the tests which resolve annotations fail with
    # this error themselves, the others still run
    with suppress(WebserviceError):
        get_registry()
    create_toy_models(MODEL_BASE_PATH)
