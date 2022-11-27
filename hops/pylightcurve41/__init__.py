import os

__version__ = open(os.path.join(os.path.abspath(os.path.dirname(__file__)), '__version__.txt')).read()
__author__ = open(os.path.join(os.path.abspath(os.path.dirname(__file__)), '__author__.txt')).read()
__author_email__ = open(os.path.join(os.path.abspath(os.path.dirname(__file__)), '__author_email__.txt')).read()
__description__ = open(os.path.join(os.path.abspath(os.path.dirname(__file__)), '__description__.txt')).read()
__url__ = open(os.path.join(os.path.abspath(os.path.dirname(__file__)), '__url__.txt')).read()

from .analysis.distributions import *
from .analysis.numerical_integration import *
from .analysis.optimisation import *
from .analysis.gaussian import *
from .analysis.curve_fit import *
from .analysis.stats import *

from .catalogues.catalogues import *
from .catalogues.simbad import *

from .images.find_stars import *

from .models.exoplanet import *
from .models.exoplanet_lc import *

from .processes.counter import *
from .processes.files import *

from .spacetime.angles import *
from .spacetime.targets import *
from .spacetime.times import *
from .spacetime.observing import *

from .errors import *
