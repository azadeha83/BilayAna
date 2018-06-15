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
from bilana import lipidmolecules
from bilana import common as com
import warnings
inputfilename_default = 'inputfile' 

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

    NUMBEROFMOLECULES = 'all'

    def __init__(self, inputfile):
        self.system_info = self.read_infofile(inputfile)
        cwd = os.getcwd()
        os.makedirs(cwd+'/datafiles', exist_ok=True)
        os.makedirs(cwd+'/indexfiles', exist_ok=True)
        os.makedirs(cwd+'/tempfiles', exist_ok=True)
        os.makedirs(cwd+'/energyfiles', exist_ok=True)
        # ''' general system information '''
        self.system = self.system_info['System']
        self.temperature = self.system_info['Temperature']
        self.cutoff = float(self.system_info['cutoff'])
        self.molecules = [x.upper().strip() for x in self.system_info['Lipidmolecules'].split(',')] # Lipid molecules in system
        if 'CHOL' in self.molecules:    # This must be declared immediately !!!
            self.molecules.append('CHL1')
            self.molecules.remove('CHOL')
        self.times = [x for x in self.system_info['Timeframe'].split(',')] # Start,End,step
        # ''' absolute_ paths to  md-files  '''
        self.mdfilepath = self.system_info['mdfiles']
        self.trjpath = '{}/md_trj/{}_{}.trr'.format(self.mdfilepath, self.system, self.temperature)
        self.gropath = '{}/initial_coords/{}.gro'.format(self.mdfilepath, self.system)
        self.toppath = '{}/psf/{}.top'.format(self.mdfilepath, self.system)
        self.tprpath = '{}/tpr/{}_{}.tpr'.format(self.mdfilepath, self.system, self.temperature)
        # ''' outputpaths (specify _absolute_ paths! '''
        self.indexpath = "{}/indexfiles".format(cwd)
        self.datapath = "{}/datafiles".format(cwd)
        self.temppath = "{}/tempfiles".format(cwd)
        self.energypath = "{}/energyfiles".format(cwd)
        # ''' Dictionaries and info '''
        self.index_to_resid, self.resid_to_lipid = self.index_conversion_dict()
        self.system_size, self.number_of_lipids = self.determine_systemsize_and_number_of_lipids()
        self.res_to_leaflet = self.assign_res_to_leaflet()
        #''' Time information '''
        if self.times[1] == 'inf':
            self.t_end = 1000000000 # Ugly but works
        else:
            self.t_end = int(self.times[1])
        self.t_start = int(self.times[0])
        self.dt = int(self.times[2])
        #####
        print('Total number of atoms: {}\nNumber of lipids: {}\n\n'
              .format(self.system_size, self.number_of_lipids))
        if self.NUMBEROFMOLECULES == 'all':
            self.NUMBEROFMOLECULES = self.number_of_lipids

    def read_infofile(self, inputfname):
        ''' Reads the inputfile. Caution!
            Input arguments are not checked for validity - double check yourself
        '''
        system_info = {}
        with open(inputfname,"r") as inputf:
            # Creates a list like [[system,dppc_chol],[temperature,290]]
            filecontent = [x.split(': ') for x in [y.strip('\n') for y in inputf.readlines()]] 
            for info in filecontent:
                if len(info) == 2 and not '#' in info[0]:
                    system_info.update({info[0].strip(): info[1].strip()})
                    print("{} : {}".format(info[0].strip(), info[1].strip()),"\n")  
        return system_info

    def index_conversion_dict(self): 
        ''' returns a dictionary for conversion
            from index to resid as well as resid to molecule
        '''
        grofile = self.gropath
        in2res = {}
        res2mol = {}
        with open(grofile,"r") as fgro:
            # get rid of first 2 lines:
            fgro.readline()
            fgro.readline() 
            for line in fgro:
                regmatch = com.GRO_format.regexp.match(line)
                if regmatch:
                    grps = regmatch.groups()
                    resid = int(grps[0].strip())
                    lipid = grps[1].strip()
                    ind = int(grps[3].strip())
                #resid = int(float(line[:5].strip()))
                #lipid = line[5:9]
                #ind = int(float(line[15:20].strip()))
                if lipid in lipidmolecules.described_molecules:
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
                regmatch = com.GRO_format.regexp.match(line)
                if regmatch:
                    grps = regmatch.groups()
                    resid = int(grps[0].strip())
                    resname = grps[1].strip()
                    if resname in self.molecules:
                        lipids_found.append(resname)
                        resids.append(resid)
            lipids_found = list(set(lipids_found))
            resids = list(set(resids))
            number_of_lipids = len(resids)
            if len(lipids_found) != len(self.molecules):
                warnings.warn("Not all lipids, given in the input file, found in structure file!") 
            return system_size, number_of_lipids
            #lines = fgro.readlines()
            #del lines[-1]
            #for item in lines:
            #    lipid = item[5:9]
            #    resid = item[:5].strip()
            #    if lipid in self.molecules:
            #        resids += [resid]
            #resids = list(set(resids))
            #resids = [int(x) for x in resids]
            #if max(resids) == len(resids):
            #    number_of_lipids = len(resids)
            #else:
            #    print('Something went wrong! Some resids are missing...')
        return system_size, number_of_lipids

    def assign_res_to_leaflet(self):
        '''
            Reads file leaflet_assignment created with function
            mainanalysis.create_leaflet_assignment_file
            and returns dictionary res:leaflet_ind
        '''
        outdict = dict()
        try:
            with open("leaflet_assignment.dat", "r") as lfile:
                lfile.readline()
                for line in lfile:
                    cols = line.split()
                    res = int(cols[0])
                    leaflet = int(cols[1])
                    outdict[res] = leaflet
        except FileNotFoundError:
            warnings.warn('File "leaflet_assignment.dat" does not exist.\n'
                  'Consider creating it using mainanalysis.create_leaflet_assignment_file()')
        return outdict


