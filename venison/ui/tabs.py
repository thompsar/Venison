"""Rail order. Analysis tabs follow the notebook's workflow order."""

from venison.ui.background import BackgroundTab
from venison.ui.console import ConsoleTab
from venison.ui.gaussians import GaussiansTab
from venison.ui.monte_carlo import MonteCarloTab
from venison.ui.project_tab import ProjectTab
from venison.ui.raw_signal import RawSignalTab
from venison.ui.results import ResultsTab
from venison.ui.stability import StabilityTab
from venison.ui.tikhonov import TikhonovTab
from venison.ui.validation import ValidationTab

TAB_CLASSES = [
    ProjectTab,
    RawSignalTab,
    BackgroundTab,
    TikhonovTab,
    ValidationTab,
    StabilityTab,
    GaussiansTab,
    MonteCarloTab,
    ResultsTab,
    ConsoleTab,  # last: pinned to the bottom of the rail
]
