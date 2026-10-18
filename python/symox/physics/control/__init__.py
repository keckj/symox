"""Linear control models, interconnections, discretization and response plots."""
from .lti import *
from .margins import phase_margin,gain_margin
from .control_plots import *
from . import lti,control_plots
__all__=lti.__all__+['phase_margin','gain_margin']+control_plots.__all__
