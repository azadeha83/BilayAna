import os
from . import analysis
from . import files
from . import log
from .definitions import lipidmolecules, structure_formats
from .systeminfo import SysInfo
os.environ["GMX_MAXBACKUP"] = "-1"

LOGGER = log.LOGGER