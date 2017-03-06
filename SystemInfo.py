


class system_info(type=lipid_systems):
    
    def __init__(self, inputfile):
        print("\n\n _____Initialize analysis______ \n\n")
        self.system_info={} # Dict to extract infos from system info file
        '''    ################################## '''
        with open(inputfile,"r") as inp:
            filecontent=[x.split(': ') for x in [y.strip('\n') for y in inp.readlines()]] #Creates a list like [[system,dppc_chol],[temperature,290]]
            for info in filecontent:
                if len(info)==2 and not '#' in info[0]:
                    self.system_info.update({info[0]:info[1]})
                    print("{} : {}".format(info[0],info[1]),"\n")
            sys.stdout.flush()
        '''   ##################################'''
        '''create all needed folders '''
        cwd=os.getcwd()
        try:
            os.mkdir('datafiles')
        except OSError as err:
            pass
        try:
            os.mkdir('indexfiles')
        except OSError as err:
            pass
        try:
            os.mkdir('tempfiles')
        except OSError as err:
            pass
        try:
            os.mkdir('energyfiles')
        except OSError as err:
            pass
            
        '''_absolute_ paths to  md-files  '''
        self.trjpath='{}/md_trj/{}_{}.trr'.format(self.mdfilepath,self.system,self.temperature)
        self.gropath='{}/initial_coords/{}.gro'.format(self.mdfilepath,self.system)
        self.toppath='{}/psf/{}.top'.format(self.mdfilepath,self.system)
        self.tprpath='{}/tpr/{}_{}'.format(self.mdfilepath,self.system,self.temperature)
        
        ''' outputpaths (specify _absolute_ paths! '''
        self.mdfilepath=self.system_info['mdfiles']
        self.indexpath="{}/indexfiles".format(cwd)
        self.datapath="{}/datafiles".format(cwd)
        self.temppath="{}/tempfiles".format(cwd)
        self.energypath="{}/energyfiles".format(cwd)

        
        ''' general system information '''
        self.system=self.system_info['System']
        self.temperature=self.system_info['Temperature']
        self.molecules=[x.upper() for x in self.system_info['Lipidmolecules'].split(',')]   #Lipid molecules in system
        if 'CHOL' in self.molecules:
            self.molecules.append('CHL1')
            self.molecules.remove('CHOL')
        self.times=[self.system_info['Timeframe'].split(',')[0],self.system_info['Timeframe'].split(',')[1],self.system_info['Timeframe'].split(',')[2]] #Start,End,step
        if self.times[1]=='end':
            print("Not yet implemented")
        else:
            self.t_end=int(self.times[1])
            self.t_start=int(self.times[0])
            self.dt=int(self.times[2])


        
        
        self.index_to_resid,self.resid_to_lipid=self.index_conversion_dict()
        self.system_size,self.number_of_lipids=self.determine_systemsize_and_number_of_lipids()
        global g_number_of_lipids
        if g_number_of_lipids == 'all':
            g_number_of_lipids=self.number_of_lipids
        print('Total number of atoms: {}\nNumber of lipids: {}\n\n'.format(self.system_size,self.number_of_lipids))
        sys.stdout.flush()
        
    #########################################################################################################
        
        
        

    
    def index_conversion_dict(self): 
        ''' returns a dictionary for conversion from index to resid as well as resid to molecule'''
        grofile=self.gropath
        in2res={}
        res2mol={}
        with open(grofile,"r") as fgro:
            fgro.readline(); fgro.readline() #get rid of first 2 lines
            for lines in fgro:
                resid=int(float(lines[:5].strip()))
                lipid=lines[5:9]
                atom=lines[10:15].strip()
                ind=int(float(lines[15:20].strip()))
                if lipid=='DPPC' or lipid=='DUPC' and atom=='P':
                    in2res.update({ind:resid})
                    res2mol.update({resid:lipid})
                elif lipid=='CHL1' and atom=='O3':
                    in2res.update({ind:resid})
                    res2mol.update({resid:lipid})
        return in2res,res2mol

    
    def determine_systemsize_and_number_of_lipids(self):
        ''' returns number of atoms in system and number of lipid molecules '''
        grofile=self.gropath
        number_of_lipids=0
        with open(grofile,"r") as fgro:
            resids=[]
            fgro.readline()                     #get rid of first line
            system_size=int(fgro.readline())    #second line is systemsize info
            lines=fgro.readlines()
            del lines[-1]
            for item in lines:
                lipid=item[5:9]
                resid=item[:5].strip()
                if lipid in self.molecules:
                    resids+=[resid]
            resids=list(set(resids))
            resids=[int(x) for x in resids]
            if max(resids) == len(resids):
                number_of_lipids=len(resids)
            else:
                print('Something went wrong: Not all lipids found.')
        sys.stdout.flush()
        return system_size,number_of_lipids

    def create_selectionfile_neighborsearch(self,resid):
        filename="{}/neighbors_of_residue{}".format(self.temppath,resid)
        with open(filename,"w") as selection:
            print("host =  resid {} and (name P O3);\n\
                allOAtoms = resname CHL1 and name O3 and not host;\n\
                allPAtoms = resname DPPC DUPC and name P and not host;\n\
                neibOs = allOAtoms and within 1.0 of host;\n\
                neibPs = allPAtoms and within 1.0 of host;\n\
                neibs = neibPs or neibOs;\n\
                neibs;".format(resid),file=selection)
        return filename

    def find_all_neighbors(self):
        ''' Returns a list of all neighbors being in the cutoff distance at least once in the trajectory. Neighborfile is !required! and is output of "determine_neighbors()" '''
        neighborfile="neighbor_info"
        neibdict={}
        with open(neighborfile,"r") as neibmap:
            fileheader=neibmap.readline() 
            for line in neibmap:
                cols=line.split()
                resid=int(cols[0])
                time=float(cols[1])
                if resid not in neibdict.keys():
                    neibdict.update({resid:{}})
                try:
                    neiblist=[int(x) for x in cols[3].split(',')]
                    neibdict[resid].update({time:neiblist})
                except IndexError:
                    neibdict[resid].update({time:[]})
                    print("No neighbors of residue {} at time {}.".format(cols[0],cols[1])) 
        sys.stdout.flush()
        return neibdict

    def determine_traj_length(self):
        trjpath=self.trjpath
        gmxarglist=[gmx_exec,'check','-f',trjpath]
        out,err=self.exec_gromacs(gmxarglist)
        err=err.split()
        endtime=err[err.index(b'Step')+1]
        return int(endtime.decode())
   
