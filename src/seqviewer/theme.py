"""How a rendered page bridges into a host application's theme.

Kept apart from any one viewer because every viewer this package grows will write
pages into the same application and should read the same stored preference.
"""

from __future__ import annotations

from dataclasses import dataclass

__all__ = ["Theme"]


@dataclass(frozen=True)
class Theme:
    """Names the page uses to bridge into a host application's theme.

    A page written with the defaults is self-contained.  An application that
    already stores a light/dark preference passes its own names so the pileup
    follows the same setting as the rest of its output.
    """

    storage_key: str = "seqviewer-theme"
    css_prefix: str = "cv"
    style_id: str = "cv-theme-bridge"
    script_id: str = "cv-theme-sync"
