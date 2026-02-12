"""Disk cache for isistools.

Uses diskcache (SQLite-backed) to cache expensive I/O results like
parsed footprints, control networks, and campt coordinate conversions.
Cache is keyed by file path + mtime so entries auto-invalidate when
source files change.
"""

from __future__ import annotations

from pathlib import Path

import diskcache

_CACHE_DIR = Path.home() / ".cache" / "isistools"


def get_cache() -> diskcache.Cache:
    """Return the shared isistools disk cache."""
    return diskcache.Cache(str(_CACHE_DIR))
