''' 
    #######    Program to analyse a membrane system  ########

    This will create files in the current folder.
    
    An input file with system information and for energy calculation a raw-.mdp-file is required.
    
    UTF-8 encoded
    
    Specify global variables:
     self.NUMBEROFPARTICLES: the last lipid molecule to calculate, if all lipids are to be calculated set self.NUMBEROFPARTICLES='all'.
     gmx_exec: path to the name of Gromacs binary.
     If you don't want backup files: Set gromacs env variable to -1 
     self.DENOMINATOR: for energy calculation, to split the neighbors of host in different groups.
    
'''

##################################################

#from datetime import datetime
import os
import subprocess
import sys
import re
from time import localtime, strftime
import numpy as np


gmx_exec = 'gmx' #'gmx_5.0.4_mpi'
os.environ["GMX_MAXBACKUP"] = "-1"


#def display_runtime(func):
#    starttime=datetime.now()
#    print("Job started at:",strftime("%c:", localtime()))
#    #func
#    print("Finished job at:",strftime("%c:", localtime()))
#    finishtime=datetime.now()
#    print("It took ",finishtime-starttime," to finish the job.")


###################################################################
#def print_to_file(str, file): # Defined for Python2.7 compatibility
#    file.write(str + '\n')
#def print(str, end):
##    print(end+str)
###################################################


class AnalyseLipidSystem:
    DENOMINATOR=40 #62
    NUMBEROFPARTICLES='all' #'all'

    
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
#        '''   ##################################'''



        
#        '''create all needed folders '''
        cwd=os.getcwd()
        try:
            os.mkdir('datafiles')
        except OSError:
            pass
        try:
            os.mkdir('indexfiles')
        except OSError:
            pass
        try:
            os.mkdir('tempfiles')
        except OSError:
            pass
        try:
            os.mkdir('energyfiles')
        except OSError:
            pass
            
            
#        ''' outputpaths (specify _absolute_ paths! '''
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


        
#        '''_absolute_ paths to  md-files  '''
        self.trjpath='{}/md_trj/{}_{}.trr'.format(self.mdfilepath,self.system,self.temperature)
        self.gropath='{}/initial_coords/{}.gro'.format(self.mdfilepath,self.system)
        self.toppath='{}/psf/{}.top'.format(self.mdfilepath,self.system)
        self.tprpath='{}/tpr/{}_{}'.format(self.mdfilepath,self.system,self.temperature)
        

        ''' atom selections  '''
        #self.tail_atoms_of={\
        #                'DPPC':[['C21','C22','C23','C24','C25', 'C26','C27', 'C28','C29', 'C210','C211', 'C212','C213', 'C214','C215', 'C216'],\
        #                    ['C31','C32','C33', 'C34','C35', 'C36','C37', 'C38','C39', 'C310','C311', 'C312','C313', 'C314','C315', 'C316']],\
        #                'CHL1':[[]],\
        #                'DUPC':[['C21','C22','C23','C24','C25', 'C26','C27', 'C28','C29', 'C210','C211', 'C212','C213', 'C214','C215', 'C216','C217','C218'],\
        #                      ['C31','C32','C33', 'C34','C35', 'C36','C37', 'C38','C39', 'C310','C311', 'C312','C313', 'C314','C315', 'C316','C317','C318']]}
           
        self.scd_tail_atoms_of={\
                        'DPPC':[['C22','C24','C26', 'C28', 'C210', 'C212', 'C214', 'C216'],\
                            ['C32', 'C34', 'C36','C38', 'C310', 'C312', 'C314', 'C316']],\
                        'CHL1':[['C3','C17']],\
                        'DUPC':[['C21','C23','C25', 'C27','C29','C211', 'C213', 'C215','C217'],\
                            ['C31','C33','C35', 'C37','C39','C311', 'C313', 'C315','C317']]}
        
        ##### For energy groups
        headhydr=['H13A', 'H13B', 'H13C', 'H14A', 'H14B', 'H14C','H15A', 'H15B', 'H15C', 'H12A', 'H12B','H11A', 'H11B','HA', 'HB', 'HS', 'HX', 'HY']
        tailcarbonsdppc=[['C22','C23','C24','C25','C26','C27', 'C28','C29', 'C210','C211', 'C212','C213', 'C214','C215', 'C216'],\
                    ['C32','C33','C34','C35','C36','C37', 'C38','C39', 'C310','C311', 'C312','C313', 'C314','C315', 'C316']]
        tailhydrdppc=[['H2R', 'H2S','H3R', 'H3S', 'H4R', 'H4S', 'H5R', 'H5S', 'H6R', 'H6S','H7R','H7S','H8R', 'H8S', 'H9R', 'H9S', 'H10R', 'H10S','H11R', 'H11S',\
           'H12R', 'H12S','H13R', 'H13S','H14R', 'H14S', 'H15R', 'H15S','H16R', 'H16S', 'H16T'],\
                  ['H2X', 'H2Y','H3X', 'H3Y','H4X', 'H4Y','H5X', 'H5Y', 'H6X', 'H6Y','H7X', 'H7Y','H8X', 'H8Y','H9X', 'H9Y','H10X', 'H10Y','H11X', 'H11Y',\
          'H12X', 'H12Y', 'H13X', 'H13Y','H14X', 'H14Y','H15X', 'H15Y','H16X', 'H16Y','H16Z']]
        tailcarbonsdupc=[tailcarbonsdppc[0]+['C217','C218'],tailcarbonsdppc[1]+['C317','C318']]
        tailhydrdupc=[tailhydrdppc[0].copy()+['H17X','H17Y','H18X','H18Y','H18Z'],tailhydrdppc[1].copy()+['H17X','H17Y','H18X','H18Y','H18Z']]
        tailhydrdupc[0].remove('H16T')
        tailhydrdupc[1].remove('H16Z')
        #####
        
        
        self.head_atoms_of={\
                        'DPPC':['C11','C12','C13','C14','C15','N','P','O11','O12','O13','O14','C1','C2','O21','C21','O22','C3','O31','O32','C31']+headhydr,\
                        'DUPC':['C11','C12','C13','C14','C15','N','P','O11','O12','O13','O14','C1','C2','O21','C21','O22','C3','O31','O32','C31']+headhydr,\
                        'CHL1':['all']}
        
        
        self.tail_atoms_of={\
                            'DPPC':[tailcarbonsdppc[0],tailhydrdppc[0],tailcarbonsdppc[1],tailhydrdppc[1]],\
                            'DUPC':[tailcarbonsdupc[0],tailhydrdupc[0],tailcarbonsdupc[1],tailhydrdupc[1]],\
                            'CHL1':['all']}
            
        #print(self.tail_atoms_of['DUPC'])           
        self.central_atom_of={\
                                'DPPC':'P',\
                                'CHL1':'O3',\
                                'DUPC':'P',\
                                'DOPC':'P'}

        #########################################################################################################
        
        self.index_to_resid,self.resid_to_lipid=self.index_conversion_dict()
        self.system_size,self.number_of_lipids=self.determine_systemsize_and_number_of_lipids()
        if self.NUMBEROFPARTICLES == 'all':
            self.NUMBEROFPARTICLES=self.number_of_lipids
        print('Total number of atoms: {}\nNumber of lipids: {}\n\n ___________________'.format(self.system_size,self.number_of_lipids))
        sys.stdout.flush()
        
        #########################################################################################################
        
    


    ############################################################################################################################
    ############################################################################################################################
    ############################################################################################################################

#    ''' preparation tools '''

    def exec_gromacs(self,cmd,inp_str=None): 
        '''arglist (cmd) is list of arguments like "['gmx cmd','-f','tprfile','-e','en.edr']".        --> inp_str must be byte b' ' !'''
        if inp_str is None:
            proc = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            out, err = proc.communicate()
        else:
            proc = subprocess.Popen(cmd,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            out, err = proc.communicate(inp_str)
        proc.wait()
        proc.stdout.close()
        proc.stderr.close()
        return out, err
    
    
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


    def find_all_neighbors(self,verbose='off'):
        ''' Returns a list of all neighbors being in the cutoff distance at least once in the trajectory. Neighborfile is !required! and is output of "determine_neighbors()" '''
        neighborfile="neighbor_info"
        #neib=[0,[]]
        #neiblist=[[]]
        neibdict={}
        with open(neighborfile,"r") as neibmap:
            neibmap.readline() 
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
                    if verbose=='on':
                        print("No neighbors of residue {} at time {}.".format(cols[0],cols[1])) 
                #if neib[0]==int(cols[0]) and int(cols[2])!=0:
                #    neib[1]+=[int(x) for x in cols[3].split(",")]
                #    resindex=int(cols[0])
                #    neiblist[resindex][1]+=neib[1]
                #    neiblist[resindex][1]=list(set(neiblist[resindex][1]))
                #    if neiblist[resindex][0] in neiblist[resindex][1]:
                #        print("Residue {} is its own neighbor.".format(cols[0]))
                #        break
                #elif int(cols[2])!=0:
                #    neib=[0,[]]
                #    neib[0]=int(cols[0])
                #    neib[1]+=[int(x) for x in cols[3].split(",")]
                #    neiblist+=[neib]
                #elif int(cols[2])==0:
                #    print("No neighbors of residue {} at time {}.".format(cols[0],cols[1]))
                #else:
                #    print("Something went wrong on line: \n '{}'".format(line))
        sys.stdout.flush()
        return neibdict

    def create_indexfile(self):
        ''' Creates an indexfile containing all indices of atom of each residue in system (resid_X) and all indices of all atoms in system.   '''
        print("\n_____Creating index file____\n")
        #OUTPUT IS:    resindex_all.ndx     | in the cwd!
        resindex_all=open("resindex_all.ndx","w")
        
        for i in range(1,self.number_of_lipids+1):
        #for i in range(1,2):
            print("Working on residue {}".format(i),end='\r')
            selectionfile=self.temppath+'/tmp_selectionfile'
            with open(selectionfile,"w") as sf:
                lipidtype=self.resid_to_lipid[i]
                if lipidtype  != 'CHL1': 
                    tailhalf12_l=[] ### to get half the tails
                    tailhalf22_l=[]
                    for molpart in self.tail_atoms_of[lipidtype]:
                        tailhalf12_l.extend(molpart[:len(molpart)//2])
                        tailhalf22_l.extend(molpart[len(molpart)//2:])
                    tailhalf12=' '.join(tailhalf12_l)
                    tailhalf22=' '.join(tailhalf22_l)
                    #lentailbyfour=len(self.tail_atoms_of[lipidtype])//4                 ### Problem for the half tail t12/t22:
                    #tailhalf12=' '.join(self.tail_atoms_of[lipidtype][:lentailbyfour]+self.tail_atoms_of[lipidtype][2*lentailbyfour:3*lentailbyfour])
                    #tailhalf22=' '.join(self.tail_atoms_of[lipidtype][lentailbyfour:2*lentailbyfour]+self.tail_atoms_of[lipidtype][3*lentailbyfour:])
                    tailatoms=[x for index in self.tail_atoms_of[lipidtype] for x in index] ##unpacking
                    headatoms=self.head_atoms_of[lipidtype]
                    
                    
                    selectionlist = [''.join(["resid_",str(i),"=resid ",str(i)," and resname ",lipidtype,";\n"])]
                    selectionlist += [''.join(["resid_h_",str(i),"=resid ",str(i)," and resname ",lipidtype," and name ",' '.join(headatoms),";\n"])]
                    selectionlist += [''.join(["resid_t_",str(i),"=resid ",str(i)," and resname ",lipidtype," and name ",' '.join(tailatoms),";\n"])]
                    selectionlist += [''.join(["resid_t12_",str(i),"=resid ",str(i)," and resname ",lipidtype," and name ",tailhalf12 ,";\n"])]
                    selectionlist += [''.join(["resid_t22_",str(i),"=resid ",str(i)," and resname ",lipidtype," and name ",tailhalf22,";\n"])]
                    selectionlist += [''.join(["resid_",str(i),";\n","resid_h_",str(i),";\n","resid_t_",str(i),";\n","resid_t12_",str(i),";\n","resid_t22_",str(i),";\n"])]
                    selectionstring=''.join(selectionlist)
                    sf.write(selectionstring)
                elif self.resid_to_lipid[i] == 'CHL1':
                    selectionstring=''.join(["resid_",str(i),"=resid ",str(i)," and resname CHL1;\nresid_",str(i),';'])
                    sf.write(selectionstring)
                    
            outputindex=self.indexpath+"/resid_"+str(i)+".ndx"
            gmx_select_arglist=[gmx_exec,'select','-s',self.gropath,'-sf',selectionfile,'-on',outputindex]
            out,err=self.exec_gromacs(gmx_select_arglist)
            
            with open("gmx_select.log","w") as logfile: 
                logfile.write(err.decode()); logfile.write(150*'_'); logfile.write(out.decode()); logfile.write(150*'_')
            ##append resid.ndx to resindex_all.ndx
            with open(outputindex,"r") as output_index:
                filecontent=output_index.readlines()
                resindex_all.write(''.join(filecontent)+'\n\n')
            
        ### To have whole system indices in one group    
        make_ndx_output=self.temppath+'/make_ndx_system.ndx' 
        gmx_make_ndx_arglist=[gmx_exec,'make_ndx','-f',self.gropath,'-o',make_ndx_output] 
        inp_str=b'keep 0\nq\n'  
        out,err=self.exec_gromacs(gmx_make_ndx_arglist,inp_str) 
        
        with open("gmx_make_ndx.log","w") as logfile, open(make_ndx_output,"r") as output:
            logfile.write(err.decode()); logfile.write(150*'_'); logfile.write(out.decode()); logfile.write(150*'_')
            filecontent=output.readlines()
            resindex_all.write(''.join(filecontent)+"\n\n")
        resindex_all.close() 

    def determine_traj_length(self):
        trjpath=self.trjpath
        gmxarglist=[gmx_exec,'check','-f',trjpath]
        out,err=self.exec_gromacs(gmxarglist)
        err=err.split()
        print(out,err)
        endtime=err[err.index(b'Step')+1]
        return int(endtime.decode())
   
   

#     def get_xy_angle(self,res,moltype,time,getcoords=None):
#         res=str(res)
#         tot_xyangle=0
#         if getcoords==None:
#             getcoords=self.trajectory_to_gro(time)[0]
#         for i in range(len(self.scd_tail_atoms_of[moltype])):
#             tailcoords+=([getcoords[(moltype,str(time),str(res),x)] for x in self.scd_tail_atoms_of[moltype][i]])
#         #tailcoords=[getcoords[''.join([lipidmolecule,time,str(res),x])] for x in self.scd_tail_atoms_of[lipidmolecule][0]]+[getcoords[''.join([lipidmolecule,time,str(res),x])] for x in self.scd_tail_atoms_of[lipidmolecule][1]]          #Saves coordinates of residue [tail1]+[tail2]
#         headcoords=np.array(getcoords[(moltype,str(time),str(res),self.central_atom_of[moltype])])
#         geocenter=self.calculate_geometriccenter(tailcoords)
#         tiltangle=np.arccos(np.dot((geocenter-headcoords),[1,0,0])/np.linalg.norm(geocenter-headcoords))
#         return (tiltangle*180/np.pi)
     
     
    def calculate_geometriccenter(self,coordinateinput):
        geocenter=[0,0,0]
        for atomcoords in coordinateinput:
            for dimension in atomcoords:
                geocenter[atomcoords.index(dimension)]+=dimension
        geocenter=[x/len(coordinateinput) for x in geocenter]
        return np.array(geocenter)
        
    ############################################################################################################################
    ############################################################################################################################
    ############################################################################################################################
    ############################################################################################################################

#    ''' ____________________________________________ main calculations tools __________________________________________'''
    
    
    
    def trajectory_to_gro(self,overwrite='off',atomlist=None,lipids='all'):
        try:
            os.mkdir(self.datapath+'/grofiles')
        except OSError:
            pass
        getcoords={}
        print('Converting trajectory-file to structure-file...\n')
        if lipids!='all':
            molecules=[lipids]
        else:
            if "CHOL" in self.molecules:
                molecules_new=self.molecules.copy()
                molecules_new.remove("CHOL")
                #print(molecules_new)
            else:
                molecules_new=self.molecules
            molecules=molecules_new
            
        for lipidmolecule in molecules:
            grofile_output=self.datapath+'/grofiles'+'/calc_scd_for'+str(lipidmolecule)+'.gro'
            sys.stdout.flush()
            
            if atomlist==None:
                atomlist=self.tail_atoms_of[lipidmolecule]+['P','O3']
            
            print(strftime("%H:%M:%S :", localtime()),"Processing {} ...".format(lipidmolecule))
            inp_str=str(lipidmolecule).encode()
            if os.path.isfile(grofile_output) and overwrite=='off':
                pass
            else:
                gmx_traj_arglist = [gmx_exec,'trjconv','-s',self.tprpath, '-f',self.trjpath,\
                                        '-o',grofile_output,\
                                        '-b', str(self.t_start), '-e', str(self.t_end),\
                                        '-dt', str(self.dt),\
                                        '-pbc', 'whole',]
                out,err=self.exec_gromacs(gmx_traj_arglist,inp_str)      #Creates a gro file containing {timeframe:resid:atoms:xyz}
                with open("gmx_traj.log","w") as logfile:
                    logfile.write(err.decode()); logfile.write(150*'_'); logfile.write(out.decode()); logfile.write(150*'_')
            with open(grofile_output,"r") as grofile:
                print(strftime("%H:%M:%S :", localtime()),"...read data output...")
                regexp=re.compile(r'[\s]*\d+'+lipidmolecule)
                for line in grofile:
                    if 't=' in line:
                        time=float(line[line.index('t=')+2:].strip())   #to get rid of obsolete decimals
                        print("...at time {}".format(time),end="\r")
                        if float(self.t_end)<time:
                            print("breaking at",time)
                            break
                    #print("Match",regexp.match(line),line[:15],end='\r')
                    #print("Time match is",float(self.t_start)<=time,end='\r')
                    if float(self.t_start)<=time and regexp.match(line)!=None:
                        #print("Reading data at",time,end='\r')
                        atom=line[9:15].strip()
                        lipidtype=line[5:9]
                        if atom not in atomlist and lipidtype not in molecules: continue
                        resid=line[:5].strip()
                        coordinates=[float(x) for x in line[20:44].split()]
                        keytuple=(str(lipidtype),time,int(resid),str(atom))
                        getcoords.update({keytuple:coordinates})
        return getcoords,time-1000.0        

     
    def scd_of_res(self,res,lipidmolecule,time,calculation_scheme='off',getcoords=None):#, neibstraightness=0,neiblist=None):
        time=float(time)
#         def calculate_distance(atomcoords1,atomcoords2):    #expect numpy arrays as np.array([x,y,z])
#             diffvector=atomcoords2-atomcoords1
#             distance=np.linalg.norm(diffvector)
#             return distance          
#         def calculate_tilt(lipidmolecule,time,res,getcoords):
#             tailcoords=[]
#             for i in range(len(self.scd_tail_atoms_of[lipidmolecule])):
#                 tailcoords+=([getcoords[''.join([lipidmolecule,time,str(res),x])] for x in self.scd_tail_atoms_of[lipidmolecule][i]])
#             headcoords=np.array(getcoords[''.join([lipidmolecule,time,str(res),self.central_atom_of[lipidmolecule]])])
#             geocenter=self.calculate_geometriccenter(tailcoords)
#             tiltangle=np.arccos(np.dot((geocenter-headcoords),[0,0,1])/np.linalg.norm(geocenter-headcoords))
#             if tiltangle>np.pi/2:
#                 tiltangle=abs(tiltangle-np.pi)
#             return tiltangle
             
        if getcoords==None:
            getcoords=self.trajectory_to_gro()[0]
            
        if calculation_scheme=='off':
            scds_of_atoms=[]
            scds_of_tails=[]
            #scds_of_tails_corrected=[]            
            for tail in self.scd_tail_atoms_of[lipidmolecule]:
                for atomindex in range(len(tail)-1): ### -1 because last res is not taken (taking index len() implies "range+1") 
                    atm1,atm2=tail[atomindex],tail[atomindex+1]    ### Attention: In the tail list only Scd-specific (every 2nd) atom is included!! Thus: atomindex-atomindex+1
                    coords_atm1,coords_atm2=np.array(getcoords[''.join([lipidmolecule,time,res,atm1])]),np.array(getcoords[''.join([lipidmolecule,time,res,atm2])])
                    diffvector=coords_atm1-coords_atm2
                    normdiffvector=np.linalg.norm(diffvector)
                    cos=np.dot(diffvector,[0,0,1])/normdiffvector
                    scds_of_atoms += [0.5 * (3 * cos**2 - 1)]
                scds_of_tails+=[sum(scds_of_atoms)/len(scds_of_atoms)] 
            totalscd=sum(scds_of_tails)/len(scds_of_tails)
            return (totalscd, )            
#         elif calculation_scheme=='Length_Based':
#             tiltangle=calculate_tilt(lipidmolecule, time, res, getcoords)
#             scds_of_atoms=[]
#             scds_of_tails=[]
#             scds_of_tails_corrected=[]            
#             for tail in self.scd_tail_atoms_of[lipidmolecule]:
#                 
#                 #1st calculate end to end length and initiate reallength
#                 startcoords=np.array([getcoords[(lipidmolecule,time,res,tail[0])]])
#                 endcoords=np.array([getcoords[(lipidmolecule,time,res,tail[-1])]])
#                 diffvector=endcoords-startcoords
#                 end_to_end_length,reallength=calculate_distance(startcoords,endcoords),0
#                 #################################
#                 
#                 for atomindex in range(len(tail)-1): ### -1 because last res is not taken (taking index len() implies "range+1") 
#                     atm1,atm2=tail[atomindex],tail[atomindex+1]    ### Attention: In the tail list only Scd-specific (every 2nd) atom is included!! Thus: atomindex-atomindex+1
#                     coords_atm1,coords_atm2=np.array(getcoords[(lipidmolecule,time,res,atm1)][0]),np.array(getcoords[(lipidmolecule,time,res,atm2)][0])
#                     diffvector=coords_atm1-coords_atm2
#                     normdiffvector=np.linalg.norm(diffvector)
#                     cos=np.dot(diffvector,[0,0,1])/normdiffvector
#                     scds_of_atoms += [0.5 * (3 * cos**2 - 1)]
#                     reallength+=calculate_distance(coords_atm1,coords_atm2)
#                 #########
#                 scds_of_tails+=[sum(scds_of_atoms)/len(scds_of_atoms)] 
#                 straightness=end_to_end_length/reallength            
#             
#                 if straightness>=0.90:
#                     order_angle=np.arccos(((2*scds_of_tails[-1]+1)/3)**0.5)
#                     new_cos=np.cos(order_angle-tiltangle)
#                     scds_of_tails[-1]=0.5 * (3 * new_cos**2 - 1)   
#             totalscd=sum(scds_of_tails)/len(scds_of_tails)
#             return (totalscd, tiltangle*180/np.pi, straightness)
#         
# 
#         elif calculation_scheme=='New':
#             
#             tiltangle=calculate_tilt(lipidmolecule, time, res, getcoords)
#             scds_of_atoms=[]
#             scds_of_tails=[]
#             scds_of_tails_corrected=[]            
#             for tail in self.scd_tail_atoms_of[lipidmolecule]:
#                 
#                 #1st calculate end to end length and initiate reallength
#                 startcoords=np.array([getcoords[(lipidmolecule,time,str(res),tail[0])]])
#                 endcoords=np.array([getcoords[(lipidmolecule,time,str(res),tail[-1])]])
#                 diffvector=endcoords-startcoords
#                 end_to_end_length,reallength=calculate_distance(startcoords,endcoords),0
#                 #################################
#                 
#                 for atomindex in range(len(tail)-1): ### -1 because last res is not taken (taking index len() implies "range+1") 
#                     atm1,atm2=tail[atomindex],tail[atomindex+1]    ### Attention: In the tail list only Scd-specific (every 2nd) atom is included!! Thus: atomindex-atomindex+1
#                     coords_atm1,coords_atm2=np.array(getcoords[(lipidmolecule,time,str(res),atm1)]),np.array(getcoords[(lipidmolecule,time,str(res),atm2)])
#                     diffvector=coords_atm1-coords_atm2
#                     normdiffvector=np.linalg.norm(diffvector)
#                     cos=np.dot(diffvector,[0,0,1])/normdiffvector
#                     scds_of_atoms += [0.5 * (3 * cos**2 - 1)]
#                     reallength+=calculate_distance(coords_atm1,coords_atm2)
#                 #########
#                 scds_of_tails+=[sum(scds_of_atoms)/len(scds_of_atoms)] 
#                 straightness=end_to_end_length/reallength
#                 if straightness>=0.95:
#                     neighbor_tilts=[]
#                     neighbors=neiblist[res][float(time)]
#                     for neib in neighbors:
#                         if neib in neibstraightness:
#                             neib_straight,neibtilt=neibstraightness[neib]
#                         else:
#                             neibtype=self.resid_to_lipid[int(neib)]
#                             scd,neibtilt,neib_straight=self.scd_of_res('Length_Based',neib,neibtype,time,getcoords)
#                             neibstraightness.update({neib:(neib_straight,neibtilt)})
#                         if neib_straight>=0.95:
#                             #straight_neighbors+=neib_straight
#                             #xyangle=self.get_xy_angle(neib,neibtype,time,getcoords)
#                             neighbor_tilts+=[neibtilt]
#                     if len(neighbor_tilts)!=0 and 0.75*len(neighbor_tilts)<=len(neighbors): #and abs(avg_tilt-tiltangle*180/np.pi)<=maxdiff:
#                         avg_tilt=sum(neighbor_tilts)/len(neighbor_tilts) 
#                         order_angle=np.arccos(((2*scds_of_tails[-1]+1)/3)**0.5)
#                         new_cos=np.cos(order_angle-avg_tilt*np.pi/180)
#                         scds_of_tails_corrected+=[0.5 * (3 * new_cos**2 - 1)]
#             if len(scds_of_tails) != len(scds_of_tails_corrected):
#                 totalscd=sum(scds_of_tails)/len(scds_of_tails)
#             else:
#                 totalscd=sum(scds_of_tails_corrected)/len(scds_of_tails_corrected)
#             return (totalscd, tiltangle*180/np.pi, straightness)        
        
              
    def create_scdfile(self,include_tilt='off',separate='on',maxdiff=0.0):           
        print("\n_____Extracting Scd values____\n\nTilt inclusion: {}\n".format(include_tilt)) 
        #resid2lipid=self.index_conversion_dict()[1] 
        getcoords,time=self.trajectory_to_gro()
        endtime=int(float(time))
        #neiblist=self.find_all_neighbors()
        if include_tilt!='off' and maxdiff==0.0:
            scd_outputfile=''.join(["scd_distribution",include_tilt,".dat"])
        elif include_tilt!='off' and maxdiff!=0.0:
            scd_outputfile=''.join(["scd_distribution",include_tilt,str(maxdiff),".dat"])
        else:
            scd_outputfile='scd_distribution.dat'
        with open(scd_outputfile,"w") as scdfile:
            print("Time \t Residue \t Scd",file=scdfile)
            print(strftime("%H:%M:%S :", localtime()),"Write to file...")
            #endtime=self.determine_traj_length()
            for t in range(self.t_start,endtime+1,self.dt):                       #the "time" variable should be the last declared time, thus the last read frame!
                time=float(t)
                print("at time {} ...".format(t),end="\r") 
                sys.stdout.flush()
                #neibstraightness={}
                for res in range(1,self.NUMBEROFPARTICLES+1):
                    lipidmolecule=self.resid_to_lipid[res]
                    totalscd = self.scd_of_res(res, lipidmolecule, time, getcoords)[0]
                    print("{} \t {} \t {} \t {}".format(time,res,lipidmolecule,totalscd),file=scdfile)
        if separate == 'on':
            self.create_scd_histogram(scd_outputfile)
                
    ############################################################################################################################
        
    def determine_neighbors(self,overwrite=True):
        print("\n____Determining neighbors____\n")
        try:
            os.mkdir(self.datapath+'/neighborfiles')
        except OSError:
            pass
        with open("neighbor_info","w") as outfile:
            outfile.write('Resid \t Time \t Number_of_neighbors \t List_of_Neighbors \n')
            for residue in range(1,self.NUMBEROFPARTICLES+1):
                print(". . . Working on residue: {} . . .".format(residue),end="\r")
                sys.stdout.flush()
                selectionfile=self.create_selectionfile_neighborsearch(residue)
                indexoutput=self.indexpath+'/neighbors_of_residue{}.ndx'.format(residue)
                datafileoutput=self.datapath+'/neighborfiles'+'/neighbors_of_residue{}.dat'.format(residue)
                if os.path.isfile(datafileoutput) and overwrite==False:
                    print("Neighbor file of residue {} already exists. Skipping.".format(residue))
                else:
                    cmdlist=[gmx_exec,'select','-s',self.tprpath,'-f',self.trjpath,'-sf',selectionfile,'-on',indexoutput,'-oi',datafileoutput,\
                                '-b',str(self.t_start), '-e', str(self.t_end),\
                                '-dt', str(self.dt),\
                             ]
                    out,err=self.exec_gromacs(cmdlist)
                    with open("gmx_select.log","w") as logfile:
                        logfile.write(err.decode()); logfile.write(150*'_'); logfile.write(out.decode()); logfile.write(150*'_')
                with open(datafileoutput,"r") as datfile:
                    for line in datfile:
                        cols=line.split()
                        time=cols.pop(0)
                        nneibs=cols.pop(0)
                        neibindeces=[int(x) for x in cols]
                        neibresid=[self.index_to_resid[x] for x in neibindeces]
                        residlist=','.join([str(x) for x in neibresid])
                        print('{} \t {} \t {} \t {}'.format(residue,time,nneibs,residlist), file=outfile)
                        
                        
    def get_neighbor_of(self,hostres,time):
        'returns list of resids of neighbors of host resid; Inputs hostres: Resid host, time: Time as integer'
        time=float(time)
        with open("neighbor_info","r") as ninfo:
            ninfo.readline()
            for line in ninfo:
                cols=line.split()
                res=str(cols[0])
                t=cols[1]   
                if res==str(hostres) and t==''.join([str(time),'00']):
                    try:
                        return cols[3].split(',')
                    except IndexError:
                        return []
        print(time,"I should never get here...")
    
    def get_res_info(self,info,infofile,res,time):
        'Returns info to residue, infos depend on file input'
        residues=self.get_neighbor_of(res,time)+[str(res)]
        outputdict={}
        time=''.join([str(time),'.0'])
        with open(infofile,"r") as infofile:
            if info=='scd':
                infofile.readline()
                for line in infofile:
                    cols=line.split()
                    t=cols[0]
                    resid=cols[1]
                    scd=cols[3]
                    if t==time and resid in residues:
                        #outputdict.update({''.join([time,'_',resid]):scd})
                        outputdict.update({resid:float(scd)})
                        residues.remove(resid)
                    if len(residues)==0:
                        return outputdict
            if info=='energy':
                pass
                      
    def rotate_by(self,coordinatelist,axisvector,angle=180):
        pass

    def mirror_molecules(self):
        pass
                                
    def calc_avgstructure(self,lipid='DPPC'):
        ncholset=set([])
        file_created={}
        neiblist=self.find_all_neighbors()
        try:
            lipidatoms=[atom for l in self.tail_atoms_of[lipid] for atom in l]+self.head_atoms_of[lipid]
        except IndexError:
            lipidatoms=[self.tail_atoms_of[lipid]]+[self.head_atoms_of[lipid]]
        getcoords,endtime=self.trajectory_to_gro(atomlist=lipidatoms,lipids=lipid) #(lipidtype,time,res,atom)
        ncholdict={}
        outnames={}
        serial={}
        for res in range(1,self.number_of_lipids+1):
            lipidtype=self.resid_to_lipid[res]
            ncholdict.update({res:{}})
            if lipidtype!=lipid:
                continue
            for time in range(self.t_start,int(endtime)+1,1000):
                lipid_time_coords=[]
                time=float(time)
                neighbors=neiblist[res][time]
                nchol=[self.resid_to_lipid[neib] for neib in neighbors].count('CHL1')
                ncholset=set(list(ncholset)+[nchol])
                if nchol not in file_created.keys(): file_created.update({nchol:False})
                if nchol not in serial.keys(): serial.update({nchol:1})
                if nchol not in ncholdict[res].keys():
                    ncholdict[res].update({nchol:{}})
                    outnames.update({nchol:"avg_structure_"+lipid+str(nchol)+".gro"})
                    
                if not file_created[nchol]:
                    f=open(outnames[nchol],"w")
                    f.write("Title\nXXXX\n")
                    f.close()
                    file_created[nchol]=True
                for atom in lipidatoms:
                    atomcoords=getcoords[(lipidtype,time,res,atom)]
                    lipid_time_coords.append(atomcoords)
                if time==self.t_start:
                    init_geocenter=np.array(self.calculate_geometriccenter(lipid_time_coords))    
                geocenter=np.array(self.calculate_geometriccenter(lipid_time_coords))    
                deltageocenter=init_geocenter-geocenter
                lipid_time_coords=list(np.array(lipid_time_coords)+deltageocenter)    
                
                for atom in lipidatoms:
                    try:
                        ncholdict[res][nchol][atom]+=[lipid_time_coords[lipidatoms.index(atom)]]
                    except KeyError:
                        ncholdict[res][nchol].update({atom:[lipid_time_coords[lipidatoms.index(atom)]]})
                        
            for Nneib in ncholdict[res].keys():
                with open(outnames[Nneib],"a") as outf:
                    for atom in lipidatoms:
                        x,y,z=np.mean(ncholdict[res][Nneib][atom],axis=0)
                        print("{:5d}{:5s}{:5s}{:5d}{:8.3f}{:8.3f}{:8.3f}{:8.4f}{:8.4f}{:8.4f}".format(res,lipid,atom,serial[Nneib],x,y,z,0.0,0.0,0.0),file=outf) ### COORDINATES IN ANGSTROM FOR PDB!!!
                        serial[Nneib]+=1
        for Nneib in sorted(ncholset):
            out,err=self.exec_gromacs(['sed','-i','s/XXXX/'+str(serial[Nneib]-1)+'/g',outnames[Nneib]])
            print(out,err)
            with open(outnames[Nneib],"a") as outf:
                outf.write("   1 1 1")
        
        #totserial=1
        #with open("avg_structure_"+total+".gro","w") as outf:
        #    outf.write("Title\nXXXX")
        #    for atom in lipidatoms:
        #        for Nneib in sorted(ncholset):
        #                for res in range(1,self.number_of_lipids+1):
        #                    if ncholdict[res][Nneib][atom][2]<=3.0:
        #                        x,y,z=np.mean(ncholdict[res][Nneib][atom],axis=0)
        #                        print("{:5d}{:5s}{:5s}{:5d}{:8.3f}{:8.3f}{:8.3f}{:8.4f}{:8.4f}{:8.4f}".format(res,lipid,atom,serial[Nneib],x,y,z,0.0,0.0,0.0),file=outf) ### COORDINATES IN ANGSTROM FOR PDB!!!
        #                        totserial+=1
        #        
            #    #ncholdict[res].append(nchol)
            #for Nneib in ncholdict[res].keys():
            #    #ncholset=set(ncholdict[res])
            ##for Nneib in ncholset:
            #    if Nneib not in lowestvar.keys():
            #        lowestvar.update({Nneib:(0,0,0,-1,0)})
            #    #occurrence=ncholdict[res].count(Nneib)/(endtime/1000.0)
            #    occurrence=0
            #    ncholset=0
            #    #print(occurrence)
            #    ncholdictar=np.array(ncholdict[res][Nneib])
            #    
            #    #ncholdictar=np.array(ncholdict[res])
            #    ncholavg=np.mean(ncholdictar)
            #    ncholvar=np.var(ncholdictar)
            #    ncholstd=np.std(ncholdictar)
            #    #if lowestvar[Nneib][3]==-1 or ncholvar<lowestvar[Nneib][3] and occurrence>=0.9:
            #    #    lowestvar[Nneib]=(res,ncholset,ncholavg,ncholvar,ncholvar)
            #    if lowestvar[Nneib][3]==-1 or abs(ncholvar-Nneib)<=abs(lowestvar[Nneib][2]-Nneib): #and occurrence>=lowestvar[Nneib][5]:
            #        lowestvar[Nneib]=(res,ncholset,ncholavg,ncholvar,ncholstd,occurrence)
            ##print(res,ncholset,ncholavg,ncholvar,ncholvar)
            #print(ncholdict[res])
        
        #for Nneib in lowestvar.keys():
        #    pass
        #    #print(Nneib,lowestvar[Nneib])
        #    #print(allN(lowestvar[Nneib][0],neiblist))
        #    #print("\n\n")
        
    ############################################################################################################################

    def calculate_energies(self,resindex_all,mdp_raw,overwrite=True,startres=1,endres=-1,parts='complete'):
        neiblist=self.find_all_neighbors()
        if endres==-1:
            endres=self.number_of_lipids
        if parts=='complete':
            molparts=["resid_"]
            parts=''
            denominator=self.DENOMINATOR
        elif parts=='head-tail':
            molparts=["resid_h_","resid_t_"]
            denominator=int(self.DENOMINATOR/2)
        elif parts=='head-tailhalfs':
            molparts=["resid_h_","resid_t12_","resid_t22_"]
            denominator=int(self.DENOMINATOR/4)
        
        
        print('\n____Rerunning MD for energyfiles, yielding xvgtables with relevant energies.____\nCaution mdp-file must not have energy_grps indicated!\n')
        print('\n Calculating for energygroups:',molparts)
        for res in range(startres,endres+1): ### Structure: For each res->get neibs->Divide neiblist in smaller groups->Create resid selection for .mdp->Make calculation   
            print('\n',strftime("%H:%M:%S :", localtime()),'Working on lipid '+str(res)+'...')
            all_N_of_res=list(set([neibs for t in neiblist[res].keys() for neibs in neiblist[res][t]]))
            n_neibs=len(all_N_of_res)
            if n_neibs % denominator == 0:
                number_of_groupfragments=(n_neibs//denominator)
            else:
                number_of_groupfragments=(n_neibs//denominator)+1
            print("Needing {} energy run(s)".format(number_of_groupfragments)) 
           
            for groupfragment in range(number_of_groupfragments):
                print("On fragment",groupfragment)
                sys.stdout.flush()
                g_energy_output=self.energypath+'/xvgtables/energies_residue'+str(res)+'_'+str(groupfragment)+parts+'.xvg'
                if os.path.isfile(g_energy_output) and overwrite==False:
                    print("Xvgtable for lipid {} part {} already exists. Will skip this calculation.".format(res,groupfragment))
                    continue
                     
                groupblockstart=groupfragment*denominator
                groupblockend=(groupfragment+1)*denominator
                

                


                
                energygroup_indeces=[res]+all_N_of_res[groupblockstart:groupblockend]
                energygroup_list=[]
                for index in energygroup_indeces:
                    if self.resid_to_lipid[index]=='CHL1':
                        energygroup_list.append(''.join(["resid_",str(index)]))
                    else:
                        for part in molparts:
                            energygroup_list.append(''.join([part,str(index)]))
                energygroup_string=' '.join(energygroup_list)        
                #energygroup_strings=' '.join([part+str(x) for x in energygroup_indeces])
                
                
                Etypes=["Coul-SR:","LJ-SR:"]
                energyselection=[]
                for interaction in Etypes:
                    counterhost=0 #for the cholesterol as it has just 1 molpart
                    for parthost in molparts:
                        if self.resid_to_lipid[res]=='CHL1' and counterhost==0:
                            parthost="resid_"
                            counterhost+=1
                        elif self.resid_to_lipid[res]=='CHL1' and counterhost!=0:
                            continue
                        for neib in all_N_of_res[groupblockstart:groupblockend]:
                            counterneib=0
                            for partneib in molparts:
                                if self.resid_to_lipid[neib]=='CHL1' and counterneib==0:
                                    partneib='resid_'
                                    counterneib+=1
                                elif self.resid_to_lipid[neib]=='CHL1' and counterneib!=0:
                                    continue
                                energyselection.append(''.join([interaction,parthost,str(res),"-",partneib,str(neib)]))
                all_relev_energies='\n'.join(energyselection+['\n'])
                
                #select_energies_coulomb='\n'.join(["Coul-SR:resid_"+str(res)+"-resid_"+str(x) for x in all_N_of_res[groupblockstart:groupblockend]])
                #select_energies_LJ='\n'.join(["LJ-SR:resid_"+str(res)+"-resid_"+str(x) for x in all_N_of_res[groupblockstart:groupblockend]])
                #select_all_relevant_energies=select_energies_coulomb+"\n"+select_energies_LJ
                
                groupfragment=str(groupfragment) 


                #Create Mdpfile:
                try:
                    os.mkdir(self.energypath+'/mdpfiles')
                except OSError:
                    pass
                mdpfile=self.energypath+'/mdpfiles/energy_mdp_recalc_resid'+str(res)+'_'+groupfragment+parts+'.mdp'
                with open(mdpfile,"w") as mdpfile_rerun, open(mdp_raw,"r") as mdpfile_raw:
                    mdp_raw_content=mdpfile_raw.readlines()
                    energygrpline='energygrps\t\t\t='+energygroup_string+'\n'
                    mdp_raw_content.append(energygrpline)
                    mdpfile_rerun.write('\n'.join(mdp_raw_content)+'\n')
                
                
                #Create TPRFILE with GROMPP:
                print(strftime("%H:%M:%S :", localtime()),'...Creating .tpr-file...')
                try:
                    os.mkdir(self.energypath+'/tprfiles')
                except OSError:
                    pass
                tprfile_energyrerun=self.energypath+'/tprfiles/mdrerun_resid'+str(res)+'_'+groupfragment+parts+'.tpr'
                mdpoutfile=self.energypath+'/mdpfiles'+'/mdrerun_resid'+str(res)+'_'+groupfragment+parts+'.mdp'
                grompp_arglist=[gmx_exec,'grompp','-f',mdpfile,'-p',self.toppath,'-c',self.gropath,'-o',tprfile_energyrerun,'-n',resindex_all,'-po',mdpoutfile]
                out,err=self.exec_gromacs(grompp_arglist)
                with open("gmx_grompp.log","a") as logfile:
                    logfile.write(err.decode()); logfile.write(150*'_'); logfile.write(out.decode()); logfile.write(150*'_')
                                
                
                #Create ENERGYFILE with mdrun -rerun:
                print(strftime("%H:%M:%S :", localtime()),'...Rerunning trajectory for energy calculation...')
                try:
                    os.mkdir(self.energypath+'/edrfiles')
                except OSError:
                    pass
                try:
                    os.mkdir(self.energypath+'/logfiles')
                except OSError:
                    pass
                energyfile_output=self.energypath+'/edrfiles/energyfile_resid'+str(res)+'_'+groupfragment+parts+'.edr'
                logoutput_file=self.energypath+'/logfiles'+'/mdrerun_resid'+str(res)+parts+'.log'
                trajoutput="traj"+str(res)+'_'+groupfragment+parts+'.trr'
                mdrun_arglist=[gmx_exec,'mdrun','-s',tprfile_energyrerun,'-rerun',self.trjpath,'-e',energyfile_output,'-o',trajoutput,'-g',logoutput_file]#,'-nt','8']
                out,err=self.exec_gromacs(mdrun_arglist)
                with open("gmx_mdrun.log","a") as logfile:
                    logfile.write(err.decode()); logfile.write(150*'_'); logfile.write(out.decode()); logfile.write(150*'_')
                #os.remove(trajoutput)

       

                #Create XVG-TABLE with all relevant energies:
                print(strftime("%H:%M:%S :", localtime()),'...Extracting all relevant energies from .edr file...')
                try:
                    os.mkdir(self.energypath+'/xvgtables')
                except OSError:
                    pass
                g_energy_output=self.energypath+'/xvgtables/energies_residue'+str(res)+'_'+groupfragment+parts+'.xvg'
                g_energy_arglist=[gmx_exec,'energy','-f',energyfile_output,'-s',tprfile_energyrerun,'-o',g_energy_output]
                inp_str=all_relev_energies.encode()
                out,err=self.exec_gromacs(g_energy_arglist,inp_str)
                with open("gmx_energy.log","a") as logfile:
                    logfile.write(err.decode()); logfile.write(150*'_'); logfile.write(out.decode()); logfile.write(150*'_')
    

    def write_energyfile(self,parts='complete'):
        
        if parts=='complete':
            molparts=[""]
            parts=''
            all_energies='all_energies'
            interactions=['']
            denominator=self.DENOMINATOR
        elif parts=='head-tail':
            molparts=["h_","t_"]
            denominator=int(self.DENOMINATOR/2)
            interactions=['head-tail','head-head','tail-tail']
            all_energies='all_energies_headtail.dat'
        elif parts=='head-tailhalfs':
            molparts=["h_","t12_","t22_"]
            denominator=int(self.DENOMINATOR/4)
            interactions=['head-tail12','tail12-tail12','head-tail22','tail22-tail22']
            all_energies='all_energies_headtailhalfs.dat'
        
        
        print("______________Writing energy file____________\n")
        neiblist=self.find_all_neighbors()
        for inter in interactions[-1:]:
            with open(all_energies,"w") as energyoutput:
                energyoutput.write("{: <10}{: <10}{: <10}{: <20}{: <20}{: <20}\n".format("Time","Host","Neighbor","VdW","Coul","Etot"))
                for resid in range(1,self.NUMBEROFPARTICLES+1):
                    print("Working on residue {}".format(resid),end="\r")
                    all_N_of_res=list(set([neibs for t in neiblist[resid].keys() for neibs in neiblist[resid][t]]))
                    n_neibs=len(all_N_of_res)
                    if n_neibs % denominator == 0:
                        number_of_groupfragments=(n_neibs//denominator)
                    else:
                        number_of_groupfragments=(n_neibs//denominator)+1
                    #print(number_of_groupfragments)
                    for part in range(number_of_groupfragments):
                        groupblockstart=part*denominator
                        groupblockend=(part+1)*denominator
                        neighbors_part_are=all_N_of_res[groupblockstart:groupblockend]
                        energypath=self.energypath+'/xvgtables/energies_residue'+str(resid)+'_'+str(part)+parts+'.xvg'
                        with open(energypath,"r") as energyfile:
                            res_to_row={}
                            for energyline in energyfile: #folderlayout is: time Coul_resHost_resNeib LJ_resHost_resNeib ...
                                energyline_cols=energyline.split()
                                if '@ s' in energyline:                     #creating a dict to know which column(energies) belong to which residue
                                    row=int(energyline_cols[1][1:])+1                 #time is at row 0 !
                                    neib=energyline_cols[3].split("resid_")[2][:-1]
                                    host=energyline_cols[3].split("resid_")[1][:-1]
                                    energytype=energyline_cols[3].split("-")[0][1:]
                                    #print(host,neib)
                                    res_to_row.update({(energytype,host,neib):row})
                                if '@' not in energyline and '#' not in energyline:     #pick correct energies from energyfile and print
                                    time=float(energyline_cols[0])
                                    for x in neighbors_part_are:
                                        counterhost=0
                                        for parthost in molparts:
                                            if self.resid_to_lipid[resid]=='CHL1' and counterhost==0:
                                                parthost=''
                                                counterhost+=1
                                            elif self.resid_to_lipid[resid]=='CHL1' and counterhost!=0:
                                                continue
                                            counterneib=0
                                            for partneib in molparts:
                                                if self.resid_to_lipid[x]=='CHL1' and counterneib==0:
                                                    partneib=''
                                                    counterneib+=1
                                                elif self.resid_to_lipid[x]=='CHL1' and counterneib!=0:
                                                    continue
                                                
                                                if parthost[:-1] == '':
                                                    interhost='w'
                                                else:
                                                    interhost=parthost[:-1]
                                                if partneib[:-1] == '':
                                                    interneib='w'
                                                else:
                                                    interneib=partneib[:-1]
                                                inter=''.join([interhost,'_',interneib])
                                                vdw=energyline_cols[res_to_row[('LJ',parthost+str(resid),partneib+str(x))]]
                                                coul=energyline_cols[res_to_row[('Coul',parthost+str(resid),partneib+str(x))]]
                                                Etot=float(vdw)+float(coul)
                                                print("{: <10}{: <10}{: <10}{: <20}{: <20}{: <20}{: <20.5f}".format(time,resid,x,inter,vdw,coul,Etot), file=energyoutput)
       

    ############################################################################################################################
           
    def radialdistribution(self):
        print("\n_____Calculating radial distribution function ____\n")
        try:
            os.mkdir(self.datapath+'/rdf')
        except OSError:
            pass
        selectstringP='name P;\nname P;\nname O3;'
        selectstringO3='name O3;\nname P;\nname O3;'
        selectionfileP=self.temppath+'/selrdfP'
        selectionfileO3=self.temppath+'/selrdfO3'
        outputfileP=self.datapath+'/rdf'+'/rdfhostP.xvg'
        outputfileO3=self.datapath+'/rdf'+'/rdfhostO3.xvg'
        print("...preparing selectionfiles...")
        with open(selectionfileP,"w") as sfP,open(selectionfileO3,"w") as sfO3:
            sfP.write(selectstringP)
            sfO3.write(selectstringO3)
        g_rdf_arglistP=[gmx_exec,'-nobackup','rdf','-f',self.trjpath,'-s',self.tprpath,'-sf',selectionfileP,'-o',outputfileP,'-xy']
        print("...first selection...")
        out,err=self.exec_gromacs(g_rdf_arglistP)
        with open("gmx_rdfP.log","w") as logfile:
            logfile.write(err.decode()); logfile.write(150*'_'); logfile.write(out.decode()); logfile.write(150*'_')
        g_rdf_arglistO3=[gmx_exec,'-nobackup','rdf','-f',self.trjpath,'-s',self.tprpath,'-sf',selectionfileO3,'-o',outputfileO3,'-xy']
        print("...second selection...")
        out,err=self.exec_gromacs(g_rdf_arglistO3)
        with open("gmx_rdfO3.log","w") as logfile:
            logfile.write(err.decode()); logfile.write(150*'_'); logfile.write(out.decode()); logfile.write(150*'_')
         
         
            
    def calc_gofr(self):
        print("\n_____Calculating radial distribution function ____\n")
        os.makedirs(self.datapath+'/rdf',exist_ok=True)
        makeselection=1
        return makeselection
        
    ############################################################################################################################
    
    def distance_distr(self):
        ''' ''' 
        def calculate_distance(atomcoords1,atomcoords2):    #expect numpy arrays as np.array([x,y,z])
            diffvector=atomcoords2-atomcoords1
            distance=np.linalg.norm(diffvector)
            return distance          
        getcoords,time=self.trajectory_to_gro()
        endtime=int(float(time))
        neiblist=self.find_all_neighbors()
        distances={}
        print(strftime("%H:%M:%S :", localtime()),"Write to file...")
        for t in range(self.t_start,endtime+1,self.dt):                       #the "time" variable should be the last declared time, thus the last read frame!
            time=float(t)
            print("at time {} ...".format(t),end="\r") 
            sys.stdout.flush()
            for res in range(1,self.NUMBEROFPARTICLES+1):
                lipidmolecule=self.resid_to_lipid[res]
                if lipidmolecule!='DPPC':
                    continue
                neighbors=neiblist[res][time]
                nchol=[self.resid_to_lipid[i] for i in neighbors].count('CHL1')
                if nchol not in distances.keys():
                    distances.update({nchol:{}})
                for tail in self.scd_tail_atoms_of[lipidmolecule]:
                    for atom in tail: 
                        atomcoords=np.array(getcoords[(lipidmolecule,time,res,atom)])
                        if atom == tail[0]:
                            firstcoords=np.array(getcoords[(lipidmolecule,time,res,atom)])
                        distance=calculate_distance(firstcoords, atomcoords)
                        try:                        
                            distances[nchol][atom].append(distance)
                        except KeyError:
                            distances[nchol].update({atom:[distance]})
        
        with open("Carbondistance_distribution.dat","w") as outf:
            for nchol in distances.keys():
                for atom in distances[nchol].keys():
                    meandistance=np.mean(distances[nchol][atom])
                    print("{: <10}{: <10}{: <10}".format(str(nchol)+str(atom[1:2]),atom[2:],meandistance),file=outf)
                
                          
                          
    def calc_averagedistance(self,distancefile):
        distdata={}
        neiblist=self.find_all_neighbors()
        with open(distancefile,"r") as df:
            df.readline()
            for line in df:
                col=line.split()
                time=float(col[0])
                if float(time)<float(self.t_start):
                    continue
                host=int(col[1])
                print(host,time,end='\r')
                neib=int(col[2])
                if self.resid_to_lipid[host]!='DPPC' and self.resid_to_lipid[neib]!='DPPC':
                    continue
                distance=float(col[3])
                try:
                    neighbors=neiblist[host][time]
                    nchol=[self.resid_to_lipid[int(i)] for i in neighbors].count('CHL1')
                except KeyError:
                    print('Attention: Neighbor file seems not to be complete time "{}" is missing'.format(time))
                try:
                    distdata[nchol]+=[distance]
                except KeyError:
                    distdata.update({nchol:[distance]})
        with open("average_distance_DPPC.dat",'w') as outputf:
            for key in distdata:
                avg=sum(distdata[key])/len(distdata[key])
                print("{0: <5} {1: <20}".format(key,avg),file=outputf) 
                              
    def calculate_distance(self):
        print("\n____Calculating distances____\n")
        neiblist=self.find_all_neighbors()
        try:
            os.mkdir(self.datapath+'/distcalc')
        except OSError:
            pass
        for i in range(1,self.NUMBEROFPARTICLES+1):
            all_N_of_res=list(set([neibs for t in neiblist[i].keys() for neibs in neiblist[i][t]]))
            sys.stdout.flush()
            print("Working on lipid ",i,'...',end='\r')
            neighbors=[str(x) for x in all_N_of_res]
            host=str(i)
            selectstring='resid '+host+' and name P O3;\n'
            for x in neighbors:
                selectstring+=''.join(['resid ',x,' and name P O3;\n'])
            selectionfile=''.join([self.temppath,'/seldist'])
            outputfile=''.join([self.datapath,'/distcalc','/dist_to_hostresid',str(i),'.xvg'])
            with open(selectionfile,"w") as sf:
                sf.write(selectstring)    
            g_pairdist_arglist=[gmx_exec,'-nobackup','pairdist','-f',self.trjpath,'-s',self.tprpath,'-sf',selectionfile,'-o',outputfile]
            out,err=self.exec_gromacs(g_pairdist_arglist)
            with open("gmx_pairdist.log","w") as logfile:
                logfile.write(err.decode()); logfile.write(150*'_'); logfile.write(out.decode()); logfile.write(150*'_')
        with open("all_distances.dat","w") as all_distances:
            # print("Time \t Host \t Neib \t distance/nm",file=all_distances)
            print("Time \t Host \t Neib \t distance/nm",file=all_distances)
            for resid in range(1,self.NUMBEROFPARTICLES+1):
                sys.stdout.flush()
                print("Working on residue {}".format(resid),end="\r")
                distancefile=self.datapath+'/distcalc/dist_to_hostresid'+str(resid)+'.xvg'
                with open(distancefile,"r") as dfile:
                    res_to_row={}
                    for distanceline in dfile: #folderlayout is: time distance_resHost_resNeib ...
                        distanceline_cols=distanceline.split()
                        if '@ s' in distanceline:                     #creating a dict to know which column(energies) belong to which residue
                            rownumber=int(distanceline_cols[1][1:])+1                 #time is at row 0 !
                            resnumber=distanceline_cols[4]
                            res_to_row.update({resnumber:rownumber})
                        if '@' not in distanceline and '#' not in distanceline:     #pick correct energies from energyfile and print
                            time_distanceline=float(distanceline_cols[0])
                            neighbors_are=neiblist[resid][float(time_distanceline)]
                            for x in neighbors_are:
                                distance=distanceline_cols[res_to_row[str(x)]]
                                print("{} \t {} \t {} \t {}".format(time_distanceline,resid,x,distance),file=all_distances)

#     def create_scdmap(self,scd_distribution,overwrite=False):
#         grofile_output=self.temppath+'/calc_scd_for'+str(lipidmolecule)+'.gro'
#         if os.path.isfile(grofile_output) and overwrite==False:
#             print("Writing trajectory as .gro")
#             inp_str=str(lipidmolecule).encode()
#             gmx_traj_arglist = [gmx_exec,'trjconv','-s',self.tprpath, '-f',self.trjpath,\
#                                 '-o',grofile_output,\
#                                 '-b', str(self.t_start), '-e', str(self.t_end),\
#                                 '-dt', str(self.dt),\
#                                 '-pbc', 'whole',\
#                                 ]
#             out,err=self.exec_gromacs(gmx_traj_arglist,inp_str)      #Creates a gro file containing {timeframe:resid:atoms:xyz}
#             with open("gmx_traj.log","w") as logfile:
#                 logfile.write(err.decode()); logfile.write(150*'_'); logfile.write(out.decode()); logfile.write(150*'_')
#             with open(grofile_output,"r") as grofile:
#                 print(strftime("%H:%M:%S :", localtime()),"...read data output...")
#                 for line in grofile:
#                     if 't=' in line:
#                         time=str(float(line[line.index('t=')+2:].strip()))   #to get rid of obsolete decimals
#                         print("...at time {}".format(time),end="\r")
#                     for atom in list_of_needed_atoms:
#                         if str(atom) in line:
#                             resid=line[:5].strip()
#                             lipidtype=line[5:9]
#                             #resid2lipid.update({int(resid):lipidtype})
#                             coordinates=[float(x) for x in line[20:44].split()]
#                             keystring=lipidtype+time+resid+atom
#         return



    ############################################################################################################################
    ############################################################################################################################
    ############################################################################################################################

#    ''' __________________________________________ Create Inputfiles ____________________________________________________ ''' 


    def create_Eofr_input(self,energyfile,distancefile):
        time_pair_to_E={}
        time_pair_to_r={}
        with open(energyfile,"r") as efile, open(distancefile,"r") as dfile:
            efile.readline()
            dfile.readline()
            for line in efile:
                cols=line.split()
                time=cols[0]
                respair=cols[1]+'-'+cols[2]
                Etot=cols[5]
                keystring=time+'_'+respair
                time_pair_to_E.update({keystring:Etot})
            for line in dfile:
                cols=line.split()
                time=cols[0]
                respair=cols[1]+'-'+cols[2]
                distance=cols[3]
                keystring=time+'_'+respair
                time_pair_to_r.update({keystring:distance})
        with open("Eofr_dat","w") as outfile:
            print("{: <10} {: <10} {: <10} {: < 10} {: <20}\n".format("Time","Host","Neib","Distance","Etot"),file=outfile)
            endtime=int(float(time))
            neiblist=self.find_all_neighbors()
            for i in range(1,self.NUMBEROFPARTICLES+1):
                for t in range(self.t_start,endtime+1,self.dt):
                    time=str(float(t))
                    neighbors_are=neiblist[i][float(t)]
                    for neib in neighbors_are:
                        respair=str(i)+'-'+str(neib)
                        Etot=time_pair_to_E[time+'_'+respair]
                        dist=time_pair_to_r[time+'_'+respair]
                        print("{: <10} {: <10} {: <10} {: < 10} {: <20.5f}".format(time,i,neib,dist,Etot),file=outfile)


    def create_EofScd_input(self,energyfile,scdfile,parts='complete'):
        print("______________Creating EofScd input file____________\n")
        if parts=='complete':
            #interactions=['']
            parts=''
        elif parts=='head-tail':
            #interactions=['head-tail','head-head','tail-tail']
            interactionskey=['h_h','h_t','t_t','h_w','t_w','w_w']
        elif parts=='head-tailhalfs':
            #interactions=['head-tail12','tail12-tail12','head-tail22','tail22-tail22']
            interactionskey=['h_h','h_t12','t12_t12','h_t22','t22_t22','h_w','t12_w','t22_w','w_w']
        neiblist=self.find_all_neighbors()
        timetoenergy={}
        timetoscd={}
        with open(energyfile,"r") as efile, open(scdfile,"r") as sfile:
            efile.readline()
            sfile.readline()
            if parts=='complete':
                for line in efile:
                    cols=line.split()
                    if int(cols[1])>int(cols[2]):
                        continue
                    time=float(cols[0])
                    endtime1=time
                    respair=(int(cols[1]),int(cols[2]))
                    #print(time,respair)
                    Etot=float(cols[5])
                    VDW=float(cols[4])
                    COUL=float(cols[3])
                    timetoenergy.update({(time,respair):(Etot,VDW,COUL)})
            else:
                for line in efile:
                    cols=line.split()
                    if int(cols[1])>int(cols[2]):
                        continue
                    time=float(cols[0])
                    endtime1=time
                    respair=(int(cols[1]),int(cols[2]))
                    print(time,respair,end='\r')
                    Etot=float(cols[6])
                    VDW=float(cols[5])
                    COUL=float(cols[4])
                    inttype=cols[3]
                    timetoenergy.update({(time,respair,inttype):(Etot,VDW,COUL)})
            for line in sfile:
                cols=line.split()
                time=float(cols[0])
                endtime2=time
                res=int(cols[1])
                scd=float(cols[3])
                timetoscd.update({(time,res):scd})
        lipidinteractions=[]
        for lipid1 in self.molecules:
            for lipid2 in self.molecules:
                if self.molecules.index(lipid2)>self.molecules.index(lipid1):
                    break
                lipidinteractions.append(''.join([lipid2,'_',lipid1]))
        outputfiles=[(''.join(['Eofscd',interaction,energyfile[13:]]),interaction) for interaction in lipidinteractions]
        endtime=int(min(endtime1,endtime2))
        for outf in outputfiles:
            with open(outf[0],"w") as out:
                print("{: <10}{: <10}{: <20}{: <10}{: <20}{: <20}{: <20}{: <20}{: <20}{: <20}{: <10}".format("Time","Host","Host_Scd","Neib","Neib_Scd","DeltaScd","AvgScd","Etot","Evdw","Ecoul", "NChol"),file=out)
                for host in range(1,self.NUMBEROFPARTICLES+1):
                    type_host=self.resid_to_lipid[host]
                    for t in range(self.t_start,endtime+1,self.dt):
                        t=float(t)
                        print("Working on residue {} at  {}".format(host,t),end="\r")
                        neighbors=neiblist[host][float(t)]
                        nchol=[self.resid_to_lipid[neib] for neib in neighbors].count('CHL1')
                        for neib in neighbors:
                            type_neib=self.resid_to_lipid[neib]
                            interaction=(''.join([type_host,'_',type_neib]),''.join([type_neib,'_',type_host]))
                            if neib<host or interaction[0]!=outf[1] or interaction[1]!=outf[1]:
                                continue
                            respair=(host,neib)
                            #type_neib=self.resid_to_lipid[int(neib)]
                            #type_pair=type_host+'_'+type_neib
                            scd_host=timetoscd[(t,host)]
                            scd_neib=timetoscd[(t,neib)]
                            delta_scd=abs(scd_host-scd_neib)
                            avg_scd=(scd_host+scd_neib)/2
                            if parts=='complete':
                                Etot,VDW,COUL=timetoenergy[(t,respair)]
                                print("{: <10}{: <10}{: <20.5f}{: <10}{: <20.5f}{: <20.5f}{: <20.5f}{: <20.5f}{: <20.5f}{: <20.5f} {}".format(t,host,scd_host,neib,scd_neib,delta_scd,avg_scd,float(Etot),float(VDW),float(COUL),nchol),file=out)
                            else:
                                for inter in interactionskey:
                                    if (t,respair,inter) in timetoenergy.keys():
                                        Etot,VDW,COUL=timetoenergy[(t,respair,inter)]
                                    else:
                                        continue
                                    print("{: <10}{: <10}{: <20.5f}{: <10}{: <20.5f}{: <15}{: <20.5f}{: <20.5f}{: <20.5f}{: <20.5f}{: <20.5f} {}".format(t,host,scd_host,neib,scd_neib,inter,delta_scd,avg_scd,float(Etot),float(VDW),float(COUL),nchol),file=out)
                            


    def create_NofScd_input(self,scdfile,neighborfile):
        time_resid_to_scd={}
        neighbors_of_host={}
        hosts_without_neib=[]
        with open(scdfile,"r") as sfile, open(neighborfile,"r") as nfile:
            sfile.readline()
            nfile.readline()
            for line in sfile:
                cols=line.split()
                time=cols[0]
                res=cols[1]
                scd=cols[3]
                keystring=time+'_'+res
                time_resid_to_scd.update({keystring:scd})
            for line in nfile:
                cols=line.split()
                time=str(float(cols[1]))
                host=cols[0]
                if int(cols[2])==0:
                    print(host,"has no neighbors at time",time)
                    hosts_without_neib+=[time+'_'+host]
                    continue
                neighbors=cols[3]
                keystring=time+'_'+host
                neighbors_of_host.update({keystring:neighbors})
        with open("Nofscd_dppc.dat","w") as outfile_dppc, open("Nofscd_chol.dat","w") as outfile_chol,open("Nofscd_dupc.dat","w") as outfile_dupc:
            print("{: <10} {: <10} {: <20} {: <20} {: <10} {: <10} {: <10}".format("Time","Host","Host_Scd","N Neighbor","N Chol","N DPPC","N DUPC"),file=outfile_dppc)
            print("{: <10} {: <10} {: <20} {: <20} {: <10} {: <10} {: <10}".format("Time","Host","Host_Scd","N Neighbor","N Chol","N DPPC","N DUPC"),file=outfile_chol)
            print("{: <10} {: <10} {: <20} {: <20} {: <10} {: <10} {: <10}".format("Time","Host","Host_Scd","N Neighbor","N Chol","N DPPC","N DUPC"),file=outfile_dupc)
            endtime=int(float(time))
            for i in range(1,self.NUMBEROFPARTICLES+1):
                restype=self.resid_to_lipid[i]
                print("Working on residue {} ".format(i),end="\r")
                for t in range(self.t_start,endtime+1,self.dt):
                    time=str(t)+'.0'
                    if time+'_'+str(i) in hosts_without_neib:
                        continue
                    n_neibs=len(neighbors_of_host[time+'_'+str(i)].split(','))
                    neibindexlist=neighbors_of_host[time+'_'+str(i)].split(',')
                    neibtypelist=[self.resid_to_lipid[int(resid)] for resid in neibindexlist]
                    nchol=neibtypelist.count('CHL1')
                    ndppc=neibtypelist.count('DPPC')
                    ndupc=neibtypelist.count('DUPC')
                    scd_host=float(time_resid_to_scd[time+'_'+str(i)])
                    if restype=='DPPC':
                        print("{: <10}{: <10}{: <20.5f}{: <20.5f}{: <10}{: <10}{: <10}".format(time,i,scd_host,float(n_neibs),nchol,ndppc,ndupc),file=outfile_dppc)
                    elif restype=='CHL1':
                        print("{: <10}{: <10}{: <20.5f}{: <20.5f}{: <10}{: <10}{: <10}".format(time,i,scd_host,float(n_neibs),nchol,ndppc,ndupc),file=outfile_chol)
                    elif restype=='DUPC':
                        print("{: <10}{: <10}{: <20.5f}{: <20.5f}{: <10}{: <10}{: <10}".format(time,i,scd_host,float(n_neibs),nchol,ndppc,ndupc),file=outfile_dupc)
    
    def create_scd_histogram(self,scdfile):
        #grofile_output=self.temppath+'/calc_scd_for'+str(self.lipidmolecule)+'.gro'      
        time_resid_to_scd={}
        with open(scdfile,"r") as sfile:
            sfile.readline()
            for line in sfile:
                cols=line.split()
                time=cols[0]
                res=cols[1]
                scd=cols[3]
                keystring=time+'_'+res
                time_resid_to_scd.update({keystring:scd})
        for lipid in self.molecules:
            if scdfile[-4:]=='.dat':
                scdfile=scdfile[:-4]
            data_output=''.join([scdfile,'_',lipid,'.dat'])
            with open(data_output,"w") as outfile:
                print("{: <10} {: <10} {: <20}".format("Time","Lipid","Lipid_Scd"),file=outfile)
                endtime=int(float(time))
                for i in range(1,self.NUMBEROFPARTICLES+1):
                    restype=self.resid_to_lipid[i]
                    if restype!=lipid:
                        continue
                    print("Working on residue {} ".format(i),end="\r")
                    for t in range(self.t_start,endtime+1,self.dt):
                        time=str(t)+'.0'
                        scd_host=float(time_resid_to_scd[time+'_'+str(i)])
                        print("{: <10}{: <10}{: <20.5f}".format(time,i,scd_host),file=outfile)



    def get_orderparamconfigs(self):
        '''Extracting the mean configurations of a lipid with respect to the number of cholesterol Neighbors''' 
#         def calculate_distance(atomcoords1,atomcoords2):    #expect numpy arrays as np.array([x,y,z])
#             diffvector=atomcoords2-atomcoords1
#             distance=np.linalg.norm(diffvector)
#             return distance          
        getcoords,time=self.trajectory_to_gro()
        endtime=int(float(time))
        neiblist=self.find_all_neighbors()
        scds_of_atoms={}
        print(strftime("%H:%M:%S :", localtime()),"Write to file...")
            #endtime=self.determine_traj_length()
        for t in range(self.t_start,endtime+1,self.dt):                       #the "time" variable should be the last declared time, thus the last read frame!
            time=str(float(t))
            print("at time {} ...".format(t),end="\r") 
            sys.stdout.flush()
            for res in range(1,self.NUMBEROFPARTICLES+1):
                lipidmolecule=self.resid_to_lipid[res]
                if lipidmolecule!='DPPC':
                    continue
                neighbors=neiblist[res][float(time)]
                nchol=[self.resid_to_lipid[i] for i in neighbors].count('CHL1')
                for tail in self.scd_tail_atoms_of[lipidmolecule]:
                    tmpscd=[]
                    for atomindex in range(len(tail)-1): ### -1 because last res is not taken (taking index len() implies "range+1") 
                        atm1,atm2=tail[atomindex],tail[atomindex+1]    ### Attention: In the tail list only Scd-specific (every 2nd) atom is included!! Thus: atomindex-atomindex+1
                        coords_atm1,coords_atm2=np.array(getcoords[(lipidmolecule,time,str(res),atm1)],np.array(getcoords[(lipidmolecule,time,str(res),atm2)]))
                        diffvector=coords_atm1-coords_atm2
                        normdiffvector=np.linalg.norm(diffvector)
                        cos=np.dot(diffvector,[0,0,1])/normdiffvector
                        tmpscd.append(0.5 * (3 * cos**2 - 1))
                    tailscd=round(sum(tmpscd)/len(tmpscd),1)
                    if nchol not in scds_of_atoms.keys():
                        scds_of_atoms.update({nchol:{tailscd:[]}})
                    try:                        
                        scds_of_atoms[nchol][tailscd].append(tmpscd)
                    except KeyError:
                        scds_of_atoms[nchol].update({tailscd:[tmpscd]})
                            
        avgscd_pernchol={}
        for nchol in scds_of_atoms.keys():
            outputfname=''.join(["orderconfigurations_",str(nchol),'.dat'])
            with open(outputfname,"w") as outf:
                print("{: <15} {: <30} {: <20}".format('avgscd','vectors','nsamples'),file=outf)
                for tailscd in scds_of_atoms[nchol].keys():
                    scdlists=scds_of_atoms[nchol][tailscd]
                    tmplists=[]
                    
                    for tailscds in scdlists:
                        for atomscd in tailscds:
                            try:
                                tmplists[tailscds.index(atomscd)].append(atomscd)
                            except IndexError:
                                tmplists.append([])
                                tmplists[tailscds.index(atomscd)].append(atomscd)                                
                    
                    avgscds=[sum(i)/len(i) for i in tmplists]
                    #avgscds=[scdcarb for item in scdlists for scdcarb in scdlists[scdlists.index(item)]]
                    
                    nsamples=len(scds_of_atoms[nchol][tailscd])
                    if nchol not in avgscd_pernchol.keys():
                        avgscd_pernchol.update({nchol:{tailscd:[]}})
                    try:
                        avgscd_pernchol[nchol][tailscd].append(avgscds)
                    except KeyError:
                        avgscd_pernchol[nchol].update({tailscd:[avgscds]})
                        
                    scdstring=''.join([" {: <15.4} ".format(scdval) for scdval in avgscds])
                    print("{: <10} {} {: <10}".format(tailscd,scdstring,nsamples),file=outf)
        with open("configurational_difference.dat","w") as outf:
            maxnchol=max(scds_of_atoms.keys())
            minnchol=min(scds_of_atoms.keys())
            print(avgscd_pernchol)
            for i in range(minnchol+1,maxnchol+1):
                print("\n",i)
                for tailscd in avgscd_pernchol[minnchol].keys():
                    print(tailscd)
                    if tailscd not in avgscd_pernchol[i].keys():
                        continue
                    else:
                        print(avgscd_pernchol[i][tailscd])
                        difflist=np.array(avgscd_pernchol[minnchol][tailscd])-np.array(avgscd_pernchol[i][tailscd])
                        print(difflist)
                        scdstring=''.join([" {: <15.4} ".format(scdval) for scdval in difflist[0]])
                        print("{: <10} {: <10} {}".format(i,tailscd,scdstring),file=outf)          
            

    def create_NLip_by_Ntot(self,neighborfile):
        lipid_pairs={}
        print(self.molecules)
        mol=self.molecules.copy()
        mol.remove('CHOL')
        for lipid in mol:
            for index in range(len(mol)):
                lipid_pairs.update({str(lipid)+'_'+str(mol[index]):{}})
        with open(neighborfile,"r") as neibfile:
            neibfile.readline()
            for line in neibfile:
                cols=line.split()
                if int(cols[2])==0:
                    continue
                for resid in cols[3].split(','):
                    host=self.resid_to_lipid[int(cols[0])]
                    neib=self.resid_to_lipid[int(resid)]
                    pair1=host+'_'+neib         
                    #pair2=neib+host         
                    time=cols[1]
                    try:
                        lipid_pairs[pair1][time].append(pair1)
                    except KeyError:
                        lipid_pairs[pair1].update({time:[]})
                        lipid_pairs[pair1][time].append(pair1)
        with open("dppc_pairs.dat","w") as pdppc, open("chol_pairs.dat","w") as pchol, open("dppc_chol_pairs.dat","w") as pdppchol:
            #for pair in sorted(lipid_pairs.keys()):
            for time in sorted(lipid_pairs['DPPC_CHL1'].keys()):
                print(time,len(lipid_pairs['DPPC_CHL1'][time]),file=pdppchol)
            for time in sorted(lipid_pairs['CHL1_CHL1'].keys()):
                print(time,0.5*len(lipid_pairs['CHL1_CHL1'][time]),file=pchol)
            for time in sorted(lipid_pairs['DPPC_DPPC'].keys()):
                print(time,0.5*len(lipid_pairs['DPPC_DPPC'][time]),file=pdppc)

#     def create_z_vs_(self,var):
#         getcoords=self.trajectory_to_gro(overwrite)
#         pass
            
            
    def H_real_vs_H_fit(self):
        pass
            
            
    ############################################################################################################################
    ############################################################################################################################
    ############################################################################################################################





#''' Execute '''
#starttime=datetime.now()
#print("Job started at:",strftime("%H:%M:%S", localtime()))
##
#obj = AnalyseLipidSystem("inputfile")
#
#obj.determine_neighbors()
#obj.create_indexfile()
#obj.calculate_scds()
#obj.calculate_energies("resindex_all","energy_recalculation.mdp")
#obj.write_energyfile()
#obj.create_EofScd_input('all_energies', 'scd_distribution.dat')
#obj.calculate_distance()
#obj.radialdistribution()



