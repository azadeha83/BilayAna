''' Gather information about system
    1. Reads in an inputfile
    2. Contains information of:
        - Creates all needed folders, if not present
        - Paths of mdfiles
        - System size: Number of molecules and particles
        - Types of lipids
        - MD times to consider
        - Reads information about which lipid lies in which leaflet
        - Conversion dictionaries: resid to type and index to resid
'''

import os
import re
import numpy as np
import MDAnalysis as mda
from . import log
from .definitions import lipidmolecules
from .definitions.structure_formats import REGEXP_GRO

LOGGER = log.LOGGER
INPUTFILENAME = 'inputfile'


class SysInfo():
    ''' Gather all relevant information about the system to analyse
    Info about:
        -Paths to MD files
        -Output paths (where to put generated files to)
        -Important dictionaries:
            index to resid --> maps index of an atom to corresponding resid
            resid to lipid --> maps resid to lipid type
            res to leaflet --> maps resid to leaflet index
        -number of lipids in system
        -number of atoms of lipids in system
        -
    '''
    LOGGER = LOGGER

    def __init__(self, inputfilename="inputfile", load_univ=True):
        LOGGER.debug("Initialize")
        self.startres = 1 # Initialize value
        self.system_info = self.read_infofile(inputfilename)
        cwd = os.getcwd()
        os.makedirs(cwd+'/datafiles/',   exist_ok=True)
        os.makedirs(cwd+'/indexfiles/',  exist_ok=True)
        os.makedirs(cwd+'/tempfiles/',   exist_ok=True)
        os.makedirs(cwd+'/energyfiles/', exist_ok=True)

        # ''' general system information '''
        self.system       = self.system_info['System']
        self.temperature  = self.system_info['Temperature']
        self.cutoff       = float(self.system_info['cutoff'])
        self.times        = [int(x) for x in self.system_info['Timeframe'].split(',')] # Start,End,step
        if len(self.times) < 3:
            self.times.append(None)
        self.molecules    = sorted([x.strip() for x in\
                                self.system_info['Lipidmolecules'].split(',')])
        if 'CHOL' in self.molecules:
            self.molecules.append('CHL1')
            self.molecules.remove('CHOL')
        if 'WSC1' in self.molecules:
            self.molecules.remove('WSC1')
            self.molecules += lipidmolecules.PROTEINS

        # ''' absolute_ paths to  md-files  '''
        self.mdfilepath = self.system_info['mdfiles']
        self.trjpath    = '{}/md_trj/{}_{}.trr'.format(self.mdfilepath, self.system, self.temperature)
        self.gropath    = '{}/initial_coords/{}.gro'.format(self.mdfilepath, self.system)
        self.toppath    = '{}/psf/{}.top'.format(self.mdfilepath, self.system)
        self.tprpath    = '{}/tpr/{}_{}.tpr'.format(self.mdfilepath, self.system, self.temperature)
        self.edrpath    = '{}/enr/{}_{}.edr'.format(self.mdfilepath, self.system, self.temperature)

        self.trjpath_whole = '{}/md_trj/{}_{}_whole.xtc'.format(self.mdfilepath, self.system, self.temperature)

        # Check if all paths exist
        for fpath in [self.gropath, self.toppath, self.tprpath]: self.check_file_exists(fpath)
        try:
            self.check_file_exists(self.trjpath)
        except FileNotFoundError:
            self.trjpath = self.trjpath.replace(".trr", ".xtc")
            self.check_file_exists(self.trjpath)

        # ''' outputpaths (specify _absolute_ paths! '''
        self.indexpath  = "{}/indexfiles/".format(cwd)
        self.datapath   = "{}/datafiles/".format(cwd)
        self.temppath   = "{}/tempfiles/".format(cwd)
        self.energypath = "{}/energyfiles/".format(cwd)

        # ''' Dictionaries and info '''
        self.index_to_resid, self.resid_to_lipid = self.index_conversion_dict()
        self.system_size, self.number_of_lipids, self.RESIDS =\
            self.determine_systemsize_and_number_of_lipids()
        self.res_to_leaflet = self.assign_res_to_leaflet()

        if load_univ:
            if os.path.isfile(self.trjpath_whole):
                LOGGER.info("Loading pbc whole trajectory into universe")
                self.universe   = mda.Universe(self.gropath, self.trjpath_whole)
            else:
                LOGGER.info("Reading raw traj")
                self.universe   = mda.Universe(self.gropath, self.trjpath)

            #''' Time information '''
            self.t_end_real = int(self.universe.trajectory[-1].time)
            if self.times[2] is None or self.times[2] < int(self.universe.trajectory.dt):
                LOGGER.warning("Time step in input file smaller than actual time step read from MDAnalysis. Using MDA.dt")
                self.dt         = int(self.universe.trajectory.dt)
            else:
                self.dt         = self.times[2]
            if int(self.times[1]) < self.t_end_real:
                self.t_end = int(self.times[1])
            else:
                self.t_end = int(self.t_end_real)
            self.t_start = int(self.times[0])

            self.RESNAMES          = list(set(self.universe.residues.resnames))
            [self.RESNAMES.remove(solvent) for  solvent in lipidmolecules.SOLVENTS if solvent in self.RESNAMES]

        # Set constants
        self.NUMBEROFMOLECULES = self.number_of_lipids
        self.MOLRANGE          = self.RESIDS

    def read_infofile(self, inputfname):
        ''' Reads the inputfile. Caution!
            Input arguments are not checked for validity - double check yourself
        '''
        system_info = {}
        with open(inputfname,"r") as inputf:
            # Creates a list like [[system,dppc_chol],[temperature,290]]
            regex = re.compile(r'^([\w,\d,\s]*):([\w,\d,\s, \., /]*)#*.*$')
            for line in inputf:
                match = regex.match(line)
                if match is not None:
                    key = match.group(1).strip()
                    item = match.group(2).strip().replace(" ", "").replace("\n", "")
                    system_info[key] = item
        return system_info

    def info(self):
        ''' Print out information on system '''
        outstr = []
        for key, val in self.system_info.items():
            outstr.append("{: <25}{: <20}\n".format(key+':', val))
        outstr.append('\n')
        outstr.append('{: <25}{: <20}\n'.format('Total number of atoms:', self.system_size))
        outstr.append('{: <25}{: <20}\n'.format('Number of lipids:', self.number_of_lipids))
        outstr.append('{: <25}{: <20}\n'.format('Residue types found:', ' '.join(self.RESNAMES)))
        outstr.append('\n')
        outstr = ''.join(outstr)
        print(outstr)
        return outstr

    def index_conversion_dict(self):
        ''' returns a dictionary for conversion
            from index to resid as well as resid to molecule
        '''
        grofile = self.gropath
        in2res = {}
        res2mol = {}
        firstres = True
        with open(grofile,"r") as fgro:
            # get rid of first 2 lines:
            fgro.readline()
            fgro.readline()
            for line in fgro:
                regmatch = REGEXP_GRO.match(line)
                if regmatch:
                    grps = regmatch.groups()
                    resid = int(grps[0].strip())
                    lipid = grps[1].strip()
                    ind = int(grps[3].strip())
                    if firstres:
                        firstres = False
                        self.startres = resid
                in2res[ind] = resid
                res2mol[resid] = lipid
        return in2res, res2mol

    def determine_systemsize_and_number_of_lipids(self):
        ''' returns number of atoms in system and number of lipid molecules '''
        grofile = self.gropath
        number_of_lipids = 0
        with open(grofile,"r") as fgro:
            resids = []
            lipids_found = []
            fgro.readline() # Gro headline
            system_size = int(fgro.readline()) # second line is systemsize info
            for line in fgro:
                regmatch = REGEXP_GRO.match(line)
                if regmatch:
                    grps = regmatch.groups()
                    resid = int(grps[0].strip())
                    resname = grps[1].strip()
                    if resname in self.molecules:
                        lipids_found.append(resname)
                        resids.append(resid)
            LOGGER.debug("Read gro file.")
            lipids_found = list(set(lipids_found))
            resids = list(set(resids))
            number_of_lipids = len(resids)
            if not number_of_lipids:
                raise ValueError("No lipid found! Name of mols wrong?")
            if len(lipids_found) != len(self.molecules):
                LOGGER.warning("Not all lipids, given in the input file, found in structure file! %s", grofile)
        LOGGER.debug("Output: %s %s %s", system_size, number_of_lipids, resids)
        return system_size, number_of_lipids, resids

    def assign_res_to_leaflet(self):
        '''
            Reads file leaflet_assignment created with function
            mainanalysis.create_leaflet_assignment_file
            and returns dictionary res:leaflet_ind
        '''
        outdict = {}
        try:
            with open("leaflet_assignment.dat", "r") as lfile:
                lfile.readline()
                for line in lfile:
                    cols = line.split()
                    res = int(cols[0])
                    leaflet = int(cols[1])
                    outdict[res] = leaflet
        except FileNotFoundError:
            LOGGER.warning('File "leaflet_assignment.dat" does not exist.\n'
                  'Consider creating it using mainanalysis.create_leaflet_assignment_file()')
        return outdict

    def check_file_exists(self, fpath):
        ''' Checks existence of mdfiles '''
        if os.path.isfile(fpath):
            return
        else:
            raise FileNotFoundError("File does not exist {}".format(fpath))

    def check_structurefile_format(self):
        ''' Check wether coordinate file from simulation has correct format and atomnaming '''
