''' Gather information about system '''

import sys
import os 

inputfilename_default = 'inputfile' 

class SysInfo():
    ''' Gather all relevant information about the system to analyse '''

    NUMBEROFMOLECULES = 'all' #'all'

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
        self.molecules = [x.upper() for x in self.system_info['Lipidmolecules'].split(',')] # Lipid molecules in system
        if 'CHOL' in self.molecules:
            self.molecules.append('CHL1')
            self.molecules.remove('CHOL')
        self.times = [x for x in self.system_info['Timeframe'].split(',')] # Start,End,step
        # '''_absolute_ paths to  md-files  '''
        self.mdfilepath = self.system_info['mdfiles']
        self.trjpath = '{}/md_trj/{}_{}.trr'.format(self.mdfilepath, self.system, self.temperature)
        self.gropath = '{}/initial_coords/{}.gro'.format(self.mdfilepath, self.system)
        self.toppath = '{}/psf/{}.top'.format(self.mdfilepath, self.system)
        self.tprpath = '{}/tpr/{}_{}'.format(self.mdfilepath, self.system, self.temperature)
        # ''' outputpaths (specify _absolute_ paths! '''
        self.indexpath = "{}/indexfiles".format(cwd)
        self.datapath = "{}/datafiles".format(cwd)
        self.temppath = "{}/tempfiles".format(cwd)
        self.energypath = "{}/energyfiles".format(cwd)
        # ''' Dictionaries '''
        self.index_to_resid, self.resid_to_lipid = self.index_conversion_dict()
        self.system_size, self.number_of_lipids = self.determine_systemsize_and_number_of_lipids()
        #''' '''
        if self.times[1] == 'end':
            raise ValueError("Not yet implemented")
        else:
            self.t_end = int(self.times[1])
            self.t_start = int(self.times[0])
            self.dt = int(self.times[2])
        print('Total number of atoms: {}\nNumber of lipids: {}\n\n'.format(self.system_size, self.number_of_lipids))
        if self.NUMBEROFMOLECULES == 'all':
            self.NUMBEROFMOLECULES = self.number_of_lipids

    def read_infofile(self, inputfname):
        ''' Reads the inputfile. Caution!
        Input arguments are not checked for validity - double check yourself'''
        system_info = {}
        with open(inputfname,"r") as inputf:
            # Creates a list like [[system,dppc_chol],[temperature,290]]
            filecontent = [x.split(': ') for x in [y.strip('\n') for y in inputf.readlines()]] 
            for info in filecontent:
                if len(info) == 2 and not '#' in info[0]:
                    system_info.update({info[0]: info[1]})
                    print("{} : {}".format(info[0], info[1]),"\n")  
        return system_info

    def index_conversion_dict(self): 
        ''' returns a dictionary for conversion from index to resid as well as resid to molecule'''
        grofile = self.gropath
        in2res = {}
        res2mol = {}
        with open(grofile,"r") as fgro:
            # get rid of first 2 lines:
            fgro.readline()
            fgro.readline() 
            for lines in fgro:
                resid = int(float(lines[:5].strip()))
                lipid = lines[5:9]
                atom = lines[10:15].strip()
                ind = int(float(lines[15:20].strip()))
                if lipid == 'DPPC' or lipid == 'DUPC' and atom == 'P':
                    in2res.update({ind:resid})
                    res2mol.update({resid:lipid})
                elif lipid == 'CHL1' and atom == 'O3':
                    in2res.update({ind:resid})
                    res2mol.update({resid:lipid})
        return in2res, res2mol

    def determine_systemsize_and_number_of_lipids(self):
        ''' returns number of atoms in system and number of lipid molecules '''
        grofile = self.gropath
        number_of_lipids = 0
        with open(grofile,"r") as fgro:
            resids = []
            fgro.readline() # get rid of first line
            system_size = int(fgro.readline()) # second line is systemsize info
            lines = fgro.readlines()
            del lines[-1]
            for item in lines:
                lipid = item[5:9]
                resid = item[:5].strip()
                if lipid in self.molecules:
                    resids += [resid]
            resids = list(set(resids))
            resids = [int(x) for x in resids]
            if max(resids) == len(resids):
                number_of_lipids = len(resids)
            else:
                print('Something went wrong! Some resids are missing...')
        return system_size, number_of_lipids

    #def determine_traj_length(self):
    #    trjpath=self.trjpath
    #    gmxarglist=[gmx_exec,'check','-f',trjpath]
    #    out,err=self.exec_gromacs(gmxarglist)
    #    err=err.split()
    #    endtime=err[err.index(b'Step')+1]
    #   return int(endtime.decode())

if os.path.isfile(inputfilename_default):
    inp = inputfilename_default
else:
    try:
        inp = sys.argv[1]
    except IndexError:
        inp = input('Specify (relative) path to inputfile\n')
try:
    mysystem = SysInfo(inp)
except FileNotFoundError:
    print("Inputfile not found.")
    raise
