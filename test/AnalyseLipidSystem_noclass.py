# #!/usr/bin/env python3

''' 
    #######    Program to analyse a membrane system  ########

    This will create files in the current folder.
    
    An input file with system information and for energy calculation a raw-.mdp-file is required.
    
    Specify global variables:
     g_number_of_lipids: the last lipid molecule to calculate in for testing, if all lipids are to be calculated set g_number_of_lipids='all'.
     gmx_exec: to the name of Gromacs binary.
     If you don't want backup files: Set gromacs env variable to -1 
     g_denominator: for energy calculation, to split the neighbors of host in different groups.

'''

##################################################

from datetime import datetime
import os
import subprocess
import sys
from time import localtime, strftime
import numpy as np

import myAnalysisTools as mytools
import myFileoutputTools as filecreator

def display_runtime(func):
    starttime=datetime.now()
    print("Job started at:",strftime("%c:", localtime()))
    func
    print("Finished job at:",strftime("%c:", localtime()))
    finishtime=datetime.now()
    print("It took ",finishtime-starttime," to finish the job.")

class AnalyseLipidSystem:
    
    gmx_exec = 'gmx' #'gmx_5.0.4_mpi'
    os.environ["GMX_MAXBACKUP"] = "-1"
    g_denominator=40 #62
    g_number_of_lipids='all' #'all'
    
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


        
        '''_absolute_ paths to  md-files  '''
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
                    
        self.central_atom_of={\
                                'DPPC':'P',\
                                'CHL1':'O3',\
                                'DUPC':'P',\
                                'DOPC':'P'}

        #########################################################################################################
        
        self.index_to_resid,self.resid_to_lipid=self.index_conversion_dict()
        self.system_size,self.number_of_lipids=self.determine_systemsize_and_number_of_lipids()
        global g_number_of_lipids
        if g_number_of_lipids == 'all':
            g_number_of_lipids=self.number_of_lipids
        print('Total number of atoms: {}\nNumber of lipids: {}\n\n ___________________'.format(self.system_size,self.number_of_lipids))
        sys.stdout.flush()
        
        #########################################################################################################
        
    


    ############################################################################################################################
    ############################################################################################################################
    ############################################################################################################################

    ''' preparation tools '''

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


    def find_all_neighbors(self):
        ''' Returns a list of all neighbors being in the cutoff distance at least once in the trajectory. Neighborfile is !required! and is output of "determine_neighbors()" '''
        neighborfile="neighbor_info"
        neib=[0,[]]
        neiblist=[[]]
        with open(neighborfile,"r") as neibmap:
            fileheader=neibmap.readline() 
            for line in neibmap:
                cols=line.split()
                if neib[0]==int(cols[0]) and int(cols[2])!=0:
                    neib[1]+=[int(x) for x in cols[3].split(",")]
                    resindex=int(cols[0])
                    neiblist[resindex][1]+=neib[1]
                    neiblist[resindex][1]=list(set(neiblist[resindex][1]))
                    if neiblist[resindex][0] in neiblist[resindex][1]:
                        print("Residue {} is its own neighbor.".format(cols[0]))
                        break
                elif int(cols[2])!=0:
                    neib=[0,[]]
                    neib[0]=int(cols[0])
                    neib[1]+=[int(x) for x in cols[3].split(",")]
                    neiblist+=[neib]
                elif int(cols[2])==0:
                    print("No neighbors of residue {} at time {}.".format(cols[0],cols[1]))
                else:
                    print("Something went wrong on line: \n '{}'".format(line))
        sys.stdout.flush()
        return neiblist

    def create_indexfile(self):
        print("\n_____Creating index file____\n")
        ''' Creates an indexfile containing all indices of atom of each residue in system (resid_X) and all indices of all atoms in system.   '''
        #OUTPUT IS:    resindex_all.ndx     | in the cwd!
        for i in range(1,self.number_of_lipids+1):
            print("Working on resiue {}".format(i),end='\r')
            selectionfile=self.temppath+'/tmp_selectionfile'
            with open(selectionfile,"w") as sf:
                if self.resid_to_lipid[i] == 'DPPC':
                    selectionstring = "resid_"+str(i)+"=resid "+str(i)+" and resname DPPC;\nresid_"+str(i)+';'
                    sf.write(selectionstring)
                elif self.resid_to_lipid[i] == 'CHL1':
                    selectionstring="resid_"+str(i)+"=resid "+str(i)+" and resname CHL1;\nresid_"+str(i)+';'
                    sf.write(selectionstring)
                elif self.resid_to_lipid[i] == 'DUPC':
                    selectionstring="resid_"+str(i)+"=resid "+str(i)+" and resname DUPC;\nresid_"+str(i)+';'
                    sf.write(selectionstring)
            outputindex=self.indexpath+"/resid_"+str(i)+".ndx"
            gmx_select_arglist=[gmx_exec,'select','-s',self.gropath,'-sf',selectionfile,'-on',outputindex]
            out,err=self.exec_gromacs(gmx_select_arglist)
            with open("gmx_select.log","w") as logfile, open(outputindex,"r") as output_index, open("resindex_all.ndx","a") as resindex_all:
                logfile.write(err.decode()); logfile.write(150*'_'); logfile.write(out.decode()); logfile.write(150*'_')
                filecontent=output_index.readlines()
                resindex_all.write(''.join(filecontent)+'\n\n')
        make_ndx_output=self.temppath+'/make_ndx_system.ndx'
        gmx_make_ndx_arglist=[gmx_exec,'make_ndx','-f',self.gropath,'-o',make_ndx_output]
        inp_str=b'keep 0\nq\n'
        out,err=self.exec_gromacs(gmx_make_ndx_arglist,inp_str)
        with open("gmx_make_ndx.log","w") as logfile, open(make_ndx_output,"r") as output, open("resindex_all.ndx","a") as resindex_all:
            logfile.write(err.decode()); logfile.write(150*'_'); logfile.write(out.decode()); logfile.write(150*'_')
            filecontent=output.readlines()
            resindex_all.write(''.join(filecontent)+"\n\n")
   

    def determine_traj_length(self):
        trjpath=self.trjpath
        gmxarglist=[gmx_exec,'check','-f',trjpath]
        out,err=self.exec_gromacs(gmxarglist)
        err=err.split()
        endtime=err[err.index(b'Step')+1]
        return int(endtime.decode())
   
   

    def get_xy_angle(self,res,moltype,time,getcoords=None):
        res=str(res)
        tot_xyangle=0
        if getcoords==None:
            getcoords=self.trajectory_to_gro(time)[0]
        '''
        for tail in self.scd_tail_atoms_of[moltype]:
            xyangle_tail=0
            for atomindex in range(len(tail)-1): ### -1 because last res is not taken (taking index len() implies "range+1") 
                atm1=tail[atomindex]    ### Attention: In the tail list only Scd-specific (every 2nd) atom is included!! Thus: atomindex-atomindex+1
                atm2=tail[atomindex+1]
                resid=str(res)
                xyz_atom1=np.array(getcoords[''.join([moltype,str(time),resid,atm1])])
                xyz_atom2=np.array(getcoords[''.join([moltype,str(time),resid,atm2])])
                normed_diffvec=(xyz_atom1-xyz_atom2)/np.linalg.norm(xyz_atom1-xyz_atom2)
                xyangle=np.arccos(np.dot(normed_diffvec,[1,0,0]))*180/np.pi
            xyangle_tail=xyangle/(len(tail)-1)              
            tot_xyangle=tot_xyangle+xyangle_tail
            
            xyz_atom11=np.array(getcoords[''.join([moltype,str(time),resid,tail[0]])])
            xyz_atom22=np.array(getcoords[''.join([moltype,str(time),resid,tail[-1]])])
            normed_diffvec2=(xyz_atom11-xyz_atom22)/np.linalg.norm(xyz_atom11-xyz_atom22)
            xyangle1=np.arccos(np.dot(normed_diffvec2,[1,0,0]))*180/np.pi
            print(xyangle1)
        tot_xyangle=tot_xyangle/len(self.scd_tail_atoms_of[moltype])
        return tot_xyangle
        '''
        tailcoords=[]
        for i in range(len(self.scd_tail_atoms_of[moltype])):
            tailcoords+=([getcoords[''.join([moltype,str(time),str(res),x])] for x in self.scd_tail_atoms_of[moltype][i]])
        #tailcoords=[getcoords[''.join([lipidmolecule,time,str(res),x])] for x in self.scd_tail_atoms_of[lipidmolecule][0]]+[getcoords[''.join([lipidmolecule,time,str(res),x])] for x in self.scd_tail_atoms_of[lipidmolecule][1]]          #Saves coordinates of residue [tail1]+[tail2]
        headcoords=np.array(getcoords[''.join([moltype,str(time),str(res),self.central_atom_of[moltype]])])
        geocenter=self.calculate_geometriccenter(tailcoords)
        tiltangle=np.arccos(np.dot((geocenter-headcoords),[1,0,0])/np.linalg.norm(geocenter-headcoords))
        return (tiltangle*180/np.pi)
    
    
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

    ''' ____________________________________________ main calculations tools __________________________________________'''
    
    
    
    
        
           
           
    


    ############################################################################################################################
    ############################################################################################################################
    ############################################################################################################################

    ''' __________________________________________ Create Inputfiles ____________________________________________________ ''' 


    