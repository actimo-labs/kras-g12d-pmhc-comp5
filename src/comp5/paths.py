"""Path display helper.

`Path.relative_to(ROOT)` raises ValueError when the path is relative while ROOT
is absolute, which turned a successful run into a traceback three separate
times during Stage 1. Always print through `rel`.
"""
from pathlib import Path


def rel(p, root: Path) -> str:
    try:
        return str(Path(p).resolve().relative_to(root))
    except (ValueError, OSError):
        return str(p)
