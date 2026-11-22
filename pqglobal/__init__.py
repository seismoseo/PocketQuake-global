"""PocketQuake-global — from an earthquake catalog anywhere to a relocation-summary notebook.

The same pipeline as PocketQuake (picking → HYPOINVERSE → HypoDD dt.ct/dt.cc → SKHASH → notebook
→ report, plus --augment), with the stations, picks and waveforms fetched through ObsPy FDSN
clients instead of the KMA NECIS portal:

    catalog (USGS ComCat, UTC)
        → the network's own picks for those events (FDSN event service)
        → the stations it picked, within a radius, with a 3-component sensor at >= 100 sps
        → FDSN event waveforms into the cluster's kma_waveforms/<event>/FDSN.{v,a}/
        → PocketQuake's orchestrator and its relocation pipeline

PocketQuake is the submodule external/PocketQuake (pinned release); it runs the pipeline in its
own submodule external/PocketQuake/external/korea-cluster-relocation. Its NECIS downloader is
never imported.
"""

__version__ = "0.2.0"

import os
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
POCKETQUAKE_DIR = os.path.join(ROOT, "external", "PocketQuake")

# Tool locations (EQNET_DIR, SKHASH_DIR, ...) from <repo>/.env, as pocketquake.sh sources its own;
# variables already set in the environment win.
_ENV = os.path.join(ROOT, ".env")
if os.path.exists(_ENV):
    with open(_ENV) as _fh:
        for _line in _fh:
            _line = _line.split("#", 1)[0].strip()
            if "=" in _line:
                _k, _v = _line.split("=", 1)
                os.environ.setdefault(_k.strip(), _v.strip().strip('"').strip("'"))

if POCKETQUAKE_DIR not in sys.path:
    sys.path.insert(0, POCKETQUAKE_DIR)
if not os.path.isdir(os.path.join(POCKETQUAKE_DIR, "pocketquake")):
    raise ImportError(f"PocketQuake submodule missing at {POCKETQUAKE_DIR}: run "
                      "`git submodule update --init --recursive --depth 1`")

from pocketquake import EQCYCLE_DIR  # noqa: E402  (the pipeline checkout PocketQuake runs)

if EQCYCLE_DIR not in sys.path:
    sys.path.insert(0, EQCYCLE_DIR)
