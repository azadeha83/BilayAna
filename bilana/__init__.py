import os
from . import analysis
from . import files
from . import log
from . import command_line as cmd
from .definitions import lipidmolecules, structure_formats
from .systeminfo import SysInfo
from .common import exec_gromacs
os.environ["GMX_MAXBACKUP"] = "-1"

LOGGER = log.LOGGER