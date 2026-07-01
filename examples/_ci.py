"""Support for running the examples quickly in CI.

The CI examples job runs every example end-to-end.  GitHub Actions sets the
``CI`` environment variable, in which case the longer-running examples swap
their full, human-meaningful evaluation counts for small ones that finish
quickly.  When run locally, the full counts are used.
"""
import os


def scaled(full, ci=500):
    """Return `full`, or the smaller `ci` value when running in CI."""
    return ci if os.getenv("CI") else full
