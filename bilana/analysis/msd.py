'''
    This module should contain a class that automates gromacs tool gmx msd
'''
from ..systeminfo import SysInfo
from ..common import exec_gromacs

class MSD(SysInfo):
    def __init__(self, inputfilename="inputfile"):
        super().__init__(inputfilename)