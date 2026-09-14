"""
ag_auth.py  -  single place the AlphaGenome API key is loaded from.

THE KEY IS NOT IN THIS REPOSITORY AND MUST NEVER BE.

This project tree is covered by daily OneDrive snapshots and a private git
repository. Anything written here leaves the machine. A previous incident in this
portfolio began with credentials sitting in an RTF file and ended in rotating
AWS, Anthropic and JWT secrets, so this is a repeat risk, not a hypothetical one.

The key therefore lives at ~/.alphagenome/api_key, mode 600, outside every synced
directory. That matches the convention already used for the trading keys in
~/.t212.

Resolution order:
    1. environment variable ALPHAGENOME_API_KEY
    2. ~/.alphagenome/api_key

Usage:
    from ag_auth import get_key, make_client
    client = make_client()
"""

from __future__ import annotations

import os
from pathlib import Path

KEY_PATH = Path.home() / ".alphagenome" / "api_key"
ENV_VAR = "ALPHAGENOME_API_KEY"


class MissingKey(RuntimeError):
    pass


def get_key() -> str:
    """Return the API key, or raise with instructions rather than a stack trace."""
    env = os.environ.get(ENV_VAR, "").strip()
    if env:
        return env

    if KEY_PATH.exists():
        key = KEY_PATH.read_text().strip()
        if key:
            mode = KEY_PATH.stat().st_mode & 0o777
            if mode & 0o077:
                raise MissingKey(
                    f"{KEY_PATH} is mode {mode:o}, readable beyond the owner. "
                    f"Run: chmod 600 {KEY_PATH}")
            return key

    raise MissingKey(
        "No AlphaGenome API key found.\n"
        f"  Expected at {KEY_PATH} (mode 600), or in ${ENV_VAR}.\n"
        "  Obtain one from https://deepmind.google.com/science/alphagenome\n"
        "  using a PERSONAL Google account. Workspace and Enterprise accounts\n"
        "  cannot generate a key. See docs/GETTING_MODEL_ACCESS.md.\n"
        "  Do NOT place the key anywhere inside this repository."
    )


def masked() -> str:
    """Safe-to-log representation. Use this in any print or log line."""
    try:
        k = get_key()
    except MissingKey:
        return "<no key>"
    return f"{k[:6]}...{k[-4:]} ({len(k)} chars)"


def make_client():
    """Construct an AlphaGenome client, with a clear error if the package is absent."""
    try:
        from alphagenome.models import dna_client
    except ImportError as exc:  # noqa: BLE001
        raise MissingKey(
            "The alphagenome package is not installed. Install with:\n"
            "  git clone https://github.com/google-deepmind/alphagenome.git\n"
            "  pip install ./alphagenome"
        ) from exc
    return dna_client.create(get_key())


if __name__ == "__main__":
    print("key:", masked())
    print("source:", "env" if os.environ.get(ENV_VAR) else str(KEY_PATH))
