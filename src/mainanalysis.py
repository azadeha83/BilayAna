''' Set of different tools for analysis
    Specify more here:
        -
 '''
import re
import os
import sys 
import numpy as np
from time import localtime, strftime

from src.systeminfo import mysystem
global mysystem
from src import lipidmolecules
from src import gromacstoolautomator as gmxauto

class Scd():
    ''' All about calculating the lipid Scd order parameter '''
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
            getcoords=gmxauto.trajectory_to_gro()[0]
            
        if calculation_scheme=='off':
            scds_of_atoms=[]
            scds_of_tails=[]
            #scds_of_tails_corrected=[]            
            for tail in lipidmolecules.scd_tail_atoms_of[lipidmolecule]:
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
        getcoords, time = gmxauto.trajectory_to_gro()
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
            for t in range(mysystem.t_start, endtime+1, mysystem.dt):                       #the "time" variable should be the last declared time, thus the last read frame!
                time=float(t)
                print("at time {} ...".format(t),end="\r") 
                sys.stdout.flush()
                #neibstraightness={}
                for res in range(1,mysystem.NUMBEROFMOLECULES+1):
                    lipidmolecule = mysystem.resid_to_lipid[res]
                    totalscd = self.scd_of_res(res, lipidmolecule, time, getcoords)[0]
                    print("{} \t {} \t {} \t {}".format(time,res,lipidmolecule,totalscd),file=scdfile)
        if separate == 'on':
            self.create_scd_histogram(scd_outputfile)

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
        for lipid in mysystem.molecules:
            if scdfile[-4:]=='.dat':
                scdfile=scdfile[:-4]
            data_output=''.join([scdfile,'_',lipid,'.dat'])
            with open(data_output,"w") as outfile:
                print("{: <10} {: <10} {: <20}".format("Time","Lipid","Lipid_Scd"),file=outfile)
                endtime=int(float(time))
                for i in range(1,mysystem.NUMBEROFMOLECULES+1):
                    restype = mysystem.resid_to_lipid[i]
                    if restype != lipid:
                        continue
                    print("Working on residue {} ".format(i),end="\r")
                    for t in range(mysystem.t_start, endtime+1, mysystem.dt):
                        time = str(t)+'.0'
                        scd_host=float(time_resid_to_scd[time+'_'+str(i)])
                        print("{: <10}{: <10}{: <20.5f}".format(time,i,scd_host),file=outfile)

def get_neighbor_of(self, hostres, time):
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


############################################################################################################################
############################################################################################################################
    ############################################################################################################################


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

def calc_gofr(self):
    print("\n_____Calculating radial distribution function ____\n")
    os.makedirs(self.datapath+'/rdf',exist_ok=True)
    makeselection=1
    return makeselection

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