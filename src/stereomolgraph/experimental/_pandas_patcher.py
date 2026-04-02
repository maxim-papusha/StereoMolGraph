# pyright: standard

# I am fully written by ChatGPT. But i am not performence critical code :)

from __future__ import annotations

import html
import logging
from collections.abc import Callable

try:
    import pandas as pd
except Exception:  # pragma: no cover - pandas is optional
    pd = None

from stereomolgraph import (
    CondensedReactionGraph,
    MolGraph,
    StereoCondensedReactionGraph,
    StereoMolGraph,
)
from stereomolgraph.ipython import View2D

log = logging.getLogger(__name__)

molRepresentation = "svg"
molJustify = "center"

_render_images_in_all_dataframes = True
_original_df_repr_html = None
_original_series_repr_html = None

_MOL_TYPES = (
    MolGraph,
    StereoMolGraph,
    CondensedReactionGraph,
    StereoCondensedReactionGraph,
)


def _get_renderer_for_frame(frame) -> str:
    renderer = getattr(frame, "_smg_renderer", None)
    if renderer:
        return renderer
    return "image" if _render_images_in_all_dataframes else "text"


def _column_contains_mol(series, max_scan: int = 50) -> bool:
    if pd is None:
        return False
    if getattr(series, "dtype", None) != "object":
        return False
    seen = 0
    for value in series.values:
        if pd.isna(value):
            continue
        if isinstance(value, _MOL_TYPES):
            return True
        seen += 1
        if seen >= max_scan:
            break
    return False


def _format_as_svg(value) -> str:
    if value is None:
        return ""
    view = View2D()
    width = int(view.width)
    height = int(view.height)
    try:
        svg = view.svg(value).replace("svg:", "")
    except Exception as exc:
        log.warning("Failed to render molecule SVG: %s", exc)
        svg = ""
    return (
        f'<div style="width: {width}px; height: {height}px; '
        f'text-align: {molJustify};" data-content="stereomolgraph/molecule">'
        f"{svg}</div>"
    )


def _format_as_text(value) -> str:
    return html.escape(str(value))


def PrintAsImageString(value):
    if molRepresentation.lower() == "svg":
        return _format_as_svg(value)
    return _format_as_text(value)


def _get_mol_formatter(renderer: str) -> Callable[[object], str]:
    if renderer == "image":
        return PrintAsImageString
    return _format_as_text


def _get_formatters_for_frame(frame, renderer: str):
    if pd is None:
        return {}
    formatters = {}
    for col in frame.columns:
        series = frame[col]
        if _column_contains_mol(series):
            formatters[col] = _get_mol_formatter(renderer)
    return formatters


def _index_level_contains_mol(
    index, level: int | None = None, max_scan: int = 50
) -> bool:
    if pd is None:
        return False
    if isinstance(index, pd.MultiIndex):
        values = index.get_level_values(level)
    else:
        values = index
    seen = 0
    for value in values:
        if pd.isna(value):
            continue
        if isinstance(value, _MOL_TYPES):
            return True
        seen += 1
        if seen >= max_scan:
            break
    return False


def _get_index_formatters(index, renderer: str):
    if pd is None:
        return {}
    formatter = _get_mol_formatter(renderer)
    if isinstance(index, pd.MultiIndex):
        level_formatters = {}
        for level in range(index.nlevels):
            if _index_level_contains_mol(index, level=level):
                level_formatters[level] = formatter
        return level_formatters
    if _index_level_contains_mol(index):
        return formatter
    return {}


def _repr_html_dataframe(self) -> str:
    renderer = _get_renderer_for_frame(self)
    formatters = _get_formatters_for_frame(self, renderer)
    index_formatters = _get_index_formatters(self.index, renderer)
    if renderer == "text" or (not formatters and not index_formatters):
        return _original_df_repr_html(self)
    styler = self.style
    if formatters:
        styler = styler.format(formatters, escape=None)
    if index_formatters:
        styler = styler.format_index(index_formatters, axis=0, escape=None)
    return styler.to_html()


def _repr_html_series(self) -> str:
    renderer = _get_renderer_for_frame(self)
    if not _column_contains_mol(self):
        return _original_series_repr_html(self)
    formatters = {self.name: _get_mol_formatter(renderer)}
    return self.to_frame().to_html(escape=False, formatters=formatters)


def patchPandas() -> None:
    global _original_df_repr_html, _original_series_repr_html
    if pd is None:
        log.warning("pandas is not installed; skipping pandas patch")
        return
    if _original_df_repr_html is None:
        _original_df_repr_html = pd.DataFrame._repr_html_
        pd.DataFrame._repr_html_ = _repr_html_dataframe
    if _original_series_repr_html is None and hasattr(pd.Series, "_repr_html_"):
        _original_series_repr_html = pd.Series._repr_html_
        pd.Series._repr_html_ = _repr_html_series


def unpatchPandas() -> None:
    global _original_df_repr_html, _original_series_repr_html
    if pd is None:
        return
    if _original_df_repr_html is not None:
        pd.DataFrame._repr_html_ = _original_df_repr_html
        _original_df_repr_html = None
    if _original_series_repr_html is not None and hasattr(pd.Series, "_repr_html_"):
        pd.Series._repr_html_ = _original_series_repr_html
        _original_series_repr_html = None


def renderImagesInAllDataFrames(images: bool = True) -> None:
    global _render_images_in_all_dataframes
    _render_images_in_all_dataframes = bool(images)


def changeMoleculeRendering(frame, renderer: str = "image") -> None:
    if renderer not in {"image", "text"}:
        raise ValueError("renderer must be 'image' or 'text'")
    setattr(frame, "_smg_renderer", renderer)


__all__ = [
    "PrintAsImageString",
    "changeMoleculeRendering",
    "molJustify",
    "molRepresentation",
    "molSize",
    "patchPandas",
    "renderImagesInAllDataFrames",
    "unpatchPandas",
]
