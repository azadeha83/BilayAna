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
import pandas as pd
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
            resid_to_lipid[resid] --> str(resname)
            index_to_resid[ndx] --> int(resid)
            res_to_leaflet[resid] --> int(0 | 1)
        -number of lipids in system
        -number of atoms of lipids in system

    '''
    LOGGER = LOGGER

    def __init__(self, inputfilename="inputfile", load_univ=True):
        LOGGER.debug("Initialize")
        self.startres = 1 # Initialize value
        self.system_info = self.read_infofile(inputfilename)
        cwd = os.getcwd()

        # ''' general system information '''
        self.system       = self.system_info['System']
        self.temperature  = self.system_info['Temperature']
        self.cutoff       = float(self.system_info['cutoff'])
        self.times        = [int(x) for x in self.system_info['Timeframe'].split(',')] # Start,End,step
        self.reference_atom_selection = self.system_info['refatomselection']
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
        self.PL_molecules = [mol for mol in self.molecules if not (lipidmolecules.is_protein(mol) or lipidmolecules.is_sterol(mol))]

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


        # ''' Dictionaries and info '''
        self.index_to_resid, self.resid_to_lipid = self.index_conversion_dict()
        self.system_size, self.number_of_lipids, self.resids, self.protein_resids, self.protein_resnames =\
            self.determine_systemsize_and_number_of_lipids()
        self.res_to_leaflet = self.assign_res_to_leaflet()

        # ''' Initialize MDAnalysis Universe '''
        if load_univ:
            if os.path.isfile(self.trjpath_whole):
                LOGGER.info("Loading pbc whole trajectory into universe")
                self.universe = mda.Universe(self.gropath, self.trjpath_whole)
            else:
                LOGGER.info("Reading raw traj")
                self.trjpath_whole = ''
                self.universe = mda.Universe(self.gropath, self.trjpath)

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

            self.RESNAMES     = list(set(self.universe.residues.resnames))
            self.SOLVENT      = [solv for solv in self.RESNAMES if solv in ["SOL", "TIP3"]]
            if len(self.SOLVENT) > 1:
                raise ValueError("Multiple water models found {}".format(self.SOLVENT))
            else:
                self.SOLVENT = self.SOLVENT[0]
            [self.RESNAMES.remove(solvent) for  solvent in lipidmolecules.SOLVENTS if solvent in self.RESNAMES]
        else:
            #''' Time information '''
            self.t_start = int(self.times[0])
            self.t_end = int(self.times[1])
            self.dt         = self.times[2]

        # ''' Set constants '''
        self.NUMBEROFMOLECULES = self.number_of_lipids
        self.MOLRANGE          = self.resids
        self.MOLRANGE_PROT     = self.protein_resids

        # ''' Do protein stuff '''
        if self.MOLRANGE_PROT:
            self.res_to_leaflet_prot, self.res_to_region_prot = self.assign_res_to_leaflet_prot()

        # ''' Create dirs and store paths '''
        os.makedirs(cwd+'/datafiles/',   exist_ok=True)
        os.makedirs(cwd+'/indexfiles/',  exist_ok=True)
        os.makedirs(cwd+'/tempfiles/',   exist_ok=True)
        os.makedirs(cwd+'/energyfiles/', exist_ok=True)

        self.indexpath  = "{}/indexfiles/".format(cwd)
        self.datapath   = "{}/datafiles/".format(cwd)
        self.temppath   = "{}/tempfiles/".format(cwd)
        self.energypath = "{}/energyfiles/".format(cwd)

        # ''' If tail reference is used, do following '''
        self.tail = self.system_info.get("tail", None)
        if self.tail is not None:
            LOGGER.warn("The tail implementation is not well tested and will not work for most of the tools included in this package")
            self.recalc_resid_dicts()

    def read_infofile(self, inputfname):
        ''' Reads the inputfile. Caution!
            Input arguments are not checked for validity - double check yourself
        '''
        system_info = {}
        with open(inputfname,"r") as inputf:
            # Creates a list like [[system,dppc_chol],[temperature,290]]
            #regex = re.compile(r'^([\w,\d,\s]*):([\w,\d,\s, \(, \), \., /]*)#*.*$')
            regex = re.compile(r'^([\w,\d,\s]*):([\w, \d, \s, \(, \), \., /]*)\s*#*.*$')
            for line in inputf:
                match = regex.match(line)
                if match is not None:
                    key = match.group(1).strip()
                    item = match.group(2).strip().replace("\n", "")
                    if key != "refatomselection":
                        item = item.replace(" ", "")
                    system_info[key] = item
        return system_info

    def info(self):
        ''' Print out information on system '''
        outstr = []
        for key, val in self.system_info.items():
            outstr.append("{: <25}{: <20}\n".format(key+':', val))
        outstr.append('\n')
        outstr.append('{: <25} {}\n'.format('String for reference selection:', self.reference_atom_selection))
        outstr.append('{: <25}{: <20}\n'.format('Total number of atoms:', self.system_size))
        outstr.append('{: <25}{: <20}\n'.format('Number of lipids:', self.number_of_lipids))
        outstr.append('{: <25}{: <20}\n'.format('Residue types found:', ' '.join(self.RESNAMES)))
        if self.protein_resnames is not None:
            outstr.append('{: <25}\n'.format('Protein found:') )
            outstr.append('{: <25}{: <20}\n'.format('... with resnames:', ' '.join(self.protein_resnames)) )
            outstr.append('{: <25}{: <20}\n'.format('... and resids:', ' '.join([str(i) for i in self.protein_resids])) )
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
        protein_resids   = None
        protein_resnames = None
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
                        if lipidmolecules.is_protein(resname):
                            if protein_resids is None:
                                protein_resids   = []
                                protein_resnames = []
                            protein_resids.append(resid)
                            protein_resnames.append(resname)
                        else:
                            lipids_found.append(resname)
                            resids.append(resid)
            LOGGER.debug("Read gro file.")
            lipids_found = list(set(lipids_found))
            if protein_resids is not None:
                protein_resnames = list( set(protein_resnames))
                protein_resids = list( set(protein_resids) )
            resids = list(set(resids))
            number_of_lipids = len(resids)
            if not number_of_lipids:
                raise ValueError("No lipid found! Name of mols wrong?")
            if len(lipids_found) != len(self.molecules):
                LOGGER.warning("Not all lipids, given in the input file, found in structure file! %s", grofile)
        LOGGER.debug("Output: %s %s %s", system_size, number_of_lipids, resids)
        return system_size, number_of_lipids, resids, protein_resids, protein_resnames

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

    def assign_res_to_leaflet_prot(self, inputfilename="leaflet_assignment_prot.csv"):
        try:
            dat = pd.read_csv(inputfilename)
        except ValueError:
            #except FileNotFoundError:
            LOGGER.warning('File "leaflet_assignment.dat" does not exist.\n'
                  'Consider creating it using mainanalysis.create_leaflet_assignment_file()')
            return None, None
        resid_to_leaflet = dat.filter(["resid", "leaflet"]).set_index("resid").T.to_dict()
        resid_to_region  = dat.filter(["resid", "region"]).set_index("resid").T.to_dict()
        return resid_to_leaflet, resid_to_region

    def check_file_exists(self, fpath):
        ''' Checks existence of mdfiles '''
        if os.path.isfile(fpath):
            return
        else:
            raise FileNotFoundError("File does not exist {}".format(fpath))

    def check_structurefile_format(self):
        ''' Check wether coordinate file from simulation has correct format and atomnaming '''

    def recalc_resid_dicts(self, refpos_per_resid=2):
        ''' This function had to be included for the tail implementation
            as now there are more than one reference position per resid

            Attributes that are changed:
                MOLRANGE = resids
                resid_to_lipid
                res_to_leaflet
                index_to_resid

        '''
        added_resids = 0

        old_to_new_resid = {}
        new_resids = []
        new_resid_to_lipid = {}
        new_res_to_leaflet = {}
        new_index_to_resid = {}

        for res in self.resids:
            resname = self.resid_to_lipid[res]
            if self.res_to_leaflet:
                leaflet = self.res_to_leaflet[res]
            else:
                leaflet = ""

            if resname in ["CHL1", "ERG"]: # Because this resname still only have one ref position per resid
                nres = res + added_resids
                old_to_new_resid[res] = [nres]
                new_resids.append(nres)
                new_res_to_leaflet[nres] = leaflet
                new_resid_to_lipid[nres] = resname

            else:
                new_res = [ (res + i + added_resids) for i in range(refpos_per_resid) ]

                old_to_new_resid[res] = new_res

                for nres in new_res:
                    new_resids.append(nres)
                    new_resid_to_lipid[nres] = resname
                    new_res_to_leaflet[nres] = leaflet

                added_resids += refpos_per_resid - 1
                LOGGER.debug("<new res> %s : <old res> %s", new_res, res)
                LOGGER.debug("added residues: %s", added_resids)

        # get new index_to_resid
        for ind, res in self.index_to_resid.items():
            if res not in self.resids:
                continue
            for nres in old_to_new_resid[res]:
                new_index_to_resid[ind] = nres

        # overwrite old dicts
        self.MOLRANGE = self.resids = new_resids
        self.resid_to_lipid = new_resid_to_lipid
        self.res_to_leaflet = new_res_to_leaflet
        self.index_to_resid = new_index_to_resid
