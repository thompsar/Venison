"""Tab rail icons: 24×24 stroke glyphs drawn for Venison, colored by currentColor."""

_PATHS = {
    "project": '<path d="M4 6.5A1.5 1.5 0 0 1 5.5 5h4l2 2.5h7A1.5 1.5 0 0 1 20 9v8.5a1.5 1.5 0 0 1-1.5 1.5h-13A1.5 1.5 0 0 1 4 17.5z"/>',
    "signal": '<path d="M3 12h2.5l2-6 3 12 3-9 2 5 1.5-2H21"/>',
    "background": '<path d="M4 5v14h16"/><path d="M6 7c2 7 5 9 13 9.5"/><path d="M6 15.5h13" stroke-dasharray="1.5 2"/>',
    "lcurve": '<path d="M4 4v16h16"/><path d="M7 6c0 6 1 9 4 10.5s6 1 9 1"/><circle cx="10" cy="15" r="1.6"/>',
    "validation": '<path d="M4 6c3 0 3 9 6 9s3-5 5-5 2 4 5 4"/><path d="M8 19h8" /><path d="M10 15v4"/>',
    "stability": '<path d="M4 20V4"/><path d="M4 20h16"/><path d="M6 18c2-9 4-11 6-6l1 2"/><path d="M15 6v12" stroke-dasharray="2 2"/>',
    "gaussians": '<path d="M3 19h18"/><path d="M4 19c3 0 3.5-12 6-12s3 12 6 12"/><path d="M12 19c2 0 2.5-6 4.5-6s2.5 6 4 6"/>',
    "montecarlo": '<circle cx="7" cy="15" r="1"/><circle cx="10" cy="10" r="1"/><circle cx="12.5" cy="13" r="1"/><circle cx="15" cy="8" r="1"/><circle cx="17" cy="12" r="1"/><circle cx="9" cy="17.5" r="1"/><path d="M4 4v16h16"/>',
    "compare": '<path d="M3 19h18"/><path d="M4 18c3 0 3-11 6-11s3 11 6 11"/><path d="M8 18c3 0 3-7 6-7s3 7 6 7" stroke-dasharray="2.2 2"/>',
    "results": '<rect x="4" y="4" width="16" height="16" rx="2"/><path d="M4 9.5h16M10 9.5V20"/>',
    "console": '<rect x="3.5" y="5" width="17" height="14" rx="2"/><path d="M7.5 10l2.5 2-2.5 2M12 15h4"/>',
}


def icon_svg(name: str, size: int = 24) -> str:
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{size}" height="{size}" '
        f'viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.6" '
        f'stroke-linecap="round" stroke-linejoin="round" '
        f'style="display:block;pointer-events:none">{_PATHS[name]}</svg>'
    )
