import AnalyseLipidSystem as main

def trajectory_to_gro(self,frame='all',overwrite='off'):
        print('Converting trajectory-file to structure-file...\n')
        if "CHOL" in self.molecules:
            molecules_new=self.molecules.copy()
            molecules_new.remove("CHOL")
            print(molecules_new)
        else:
            molecules_new=self.molecules
        try:
            os.mkdir(self.datapath+'/grofiles')
        except OSError:
            pass
        getcoords={}
        #lipid={}
        list_of_needed_atoms=[] 
        for lipidmolecule in molecules_new:
            grofile_output=self.datapath+'/grofiles'+'/calc_scd_for'+str(lipidmolecule)+'.gro'
            #if os.path.isfile(grofile_output) and overwrite=='off':
            #        print(".gro file of residue {} already exists. Skipping.".format(lipidmolecule))
            #else:
            sys.stdout.flush()
            for item in self.scd_tail_atoms_of[lipidmolecule]:
                list_of_needed_atoms+=item
            print(strftime("%H:%M:%S :", localtime()),"Processing {} ...".format(lipidmolecule))
            inp_str=str(lipidmolecule).encode()
             # gmx_traj_arglist = [gmx_exec,'traj','-s',self.tprpath, '-f',self.trjpath,'-oxt',grofile_output]
            if os.path.isfile(grofile_output) and overwrite=='off':
                #print('File exists, will not overwrite',grofile_output)
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
                for line in grofile:
                    if 't=' in line:
                        time=str(float(line[line.index('t=')+2:].strip()))   #to get rid of obsolete decimals
                        print("...at time {}".format(time),end="\r")
                    if frame!='all' and float(time)<float(frame):
                        continue
                    elif frame!='all' and float(time)>float(frame):
                        break
                    elif float(self.t_end)>=float(time):
                        for atom in (list_of_needed_atoms+['P','O3']):
                            if str(atom) in line[9:]:
                                resid=line[:5].strip()
                                lipidtype=line[5:9]
                                #resid2lipid.update({int(resid):lipidtype})
                                coordinates=[float(x) for x in line[20:44].split()]
                                keystring=lipidtype+time+resid+atom
                                getcoords.update({keystring:coordinates})
                    else:
                        break
        return getcoords,float(time)-1000.0        

     
    def scd_of_res(self,calculation_scheme,res,lipidmolecule,time,getcoords=None, maxdiff=0.0):
        normtotalcounter=0
        def calculate_distance(atomcoords1,atomcoords2):    #expect numpy arrays as np.array([x,y,z])
            diffvector=atomcoords2-atomcoords1
            distance=np.linalg.norm(diffvector)
            return distance          
        totalscd=0
        if getcoords==None:
            getcoords=self.trajectory_to_gro(time)[0]
        if calculation_scheme=='off':
            for tail in self.scd_tail_atoms_of[lipidmolecule]:
                scd_of_tail=0
                for atomindex in range(len(tail)-1): ### -1 because last res is not taken (taking index len() implies "range+1") 
                    atm1=tail[atomindex]    ### Attention: In the tail list only Scd-specific (every 2nd) atom is included!! Thus: atomindex-atomindex+1
                    atm2=tail[atomindex+1]
                    resid=str(res)
                    coordinates_of_atom1=np.array(getcoords[''.join([lipidmolecule,time,resid,atm1])])
                    coordinates_of_atom2=np.array(getcoords[''.join([lipidmolecule,time,resid,atm2])])
                    diffvector=coordinates_of_atom1-coordinates_of_atom2
                    normdiffvector=np.linalg.norm(diffvector)
                    cos=np.dot(diffvector,[0,0,1])/normdiffvector
                    scd_of_atoms = 0.5 * (3 * cos**2 - 1)
                    scd_of_tail += scd_of_atoms
                scd_of_tail=scd_of_tail/(len(tail)-1)  
                totalscd=totalscd+scd_of_tail
            totalscd=totalscd/len(self.scd_tail_atoms_of[lipidmolecule])
            return (totalscd,)
        
        elif calculation_scheme=='Length_Based':
            for tail in self.scd_tail_atoms_of[lipidmolecule]:
                reallength=0
                startcoords=np.array([getcoords[''.join([lipidmolecule,time,str(res),tail[0]])]])
                endcoords=np.array([getcoords[''.join([lipidmolecule,time,str(res),tail[-1]])]])
                diffvector=endcoords-startcoords
                end_to_end_distance=calculate_distance(startcoords,endcoords)
                cos_tilt=np.dot(diffvector,[0,0,1])/np.linalg.norm(diffvector)
                tiltangle=np.arccos(cos_tilt)
                if tiltangle>np.pi/2:
                    tiltangle=abs(tiltangle-np.pi)
                scd_of_tail=0
                for atomindex in range(len(tail)-1): ### -1 because last res is not taken (taking index len() implies "range+1") 
                    atm1=tail[atomindex]    ### Attention: In the tail list only Scd-specific (every 2nd) atom is included!! Thus: atomindex-atomindex+1
                    atm2=tail[atomindex+1]
                    resid=str(res)
                    coordinates_of_atom1=np.array(getcoords[''.join([lipidmolecule,time,resid,atm1])])
                    coordinates_of_atom2=np.array(getcoords[''.join([lipidmolecule,time,resid,atm2])])
                    diffvector=coordinates_of_atom1-coordinates_of_atom2
                    normdiffvector=np.linalg.norm(diffvector)
                    cos=np.dot(diffvector,[0,0,1])/normdiffvector
                    reallength+=calculate_distance(coordinates_of_atom1,coordinates_of_atom2)
                    scd_of_atoms = 0.5 * (3 * cos**2 - 1)
                    scd_of_tail += scd_of_atoms
                scd_of_tail=scd_of_tail/(len(tail)-1)
                straightness=end_to_end_distance/reallength
                if straightness>=0.90:
                    #print('TILT!',res)
                    order_angle=np.arccos(((2*scd_of_tail+1)/3)**0.5)
                    #new_cos=(np.cos(order_angle-tiltangle)[0],np.cos(order_angle+tiltangle)[0])        
                    new_cos=np.cos(order_angle-tiltangle)[0]
                    scd_of_tail=0.5 * (3 * new_cos**2 - 1)    
                totalscd=totalscd+scd_of_tail
            totalscd=totalscd/len(self.scd_tail_atoms_of[lipidmolecule])
            return (totalscd, tiltangle*180/np.pi, straightness)
        
        elif calculation_scheme=='COM_Based':
            ### Calculates Center of mass of tail atoms and rescales ALL order angles by ang(COMxyz-Pxyz,[0,0,1])
            tailcoords=[]
            for i in range(len(self.scd_tail_atoms_of[lipidmolecule])):
                tailcoords+=([getcoords[''.join([lipidmolecule,time,str(res),x])] for x in self.scd_tail_atoms_of[lipidmolecule][i]])
            #tailcoords=[getcoords[''.join([lipidmolecule,time,str(res),x])] for x in self.scd_tail_atoms_of[lipidmolecule][0]]+[getcoords[''.join([lipidmolecule,time,str(res),x])] for x in self.scd_tail_atoms_of[lipidmolecule][1]]          #Saves coordinates of residue [tail1]+[tail2]
            headcoords=np.array(getcoords[''.join([lipidmolecule,time,str(res),self.central_atom_of[lipidmolecule]])])
            geocenter=self.calculate_geometriccenter(tailcoords)
            tiltangle=np.arccos(np.dot((geocenter-headcoords),[0,0,1])/np.linalg.norm(geocenter-headcoords))
            if tiltangle>np.pi/2:
                tiltangle=abs(tiltangle-np.pi)
            #print("Tilt angle for residue {} is {}".format(res,tiltangle*180/np.pi),end='\r')
            for tail in self.scd_tail_atoms_of[lipidmolecule]:
                scd_of_tail=0
                for atomindex in range(len(tail)-1): ### -1 because last res is not taken (taking index len() implies "range+1") 
                    atm1=tail[atomindex]    ### Attention: In the tail list only Scd-specific (every 2nd) atom is included!! Thus: atomindex-atomindex+1
                    atm2=tail[atomindex+1]
                    resid=str(res)
                    coordinates_of_atom1=np.array(getcoords[''.join([lipidmolecule,time,resid,atm1])])
                    coordinates_of_atom2=np.array(getcoords[''.join([lipidmolecule,time,resid,atm2])])
                    diffvector=coordinates_of_atom1-coordinates_of_atom2
                    normdiffvector=np.linalg.norm(diffvector)
                    cos=np.dot(diffvector,[0,0,1])/normdiffvector
                    chain_angle=np.arccos((cos**2)**0.5)
                    #if lipidmolecule=='CHL1':
                    #    print(res, chain_angle*180/np.pi,tiltangle*180/np.pi, chain_angle*180/np.pi-tiltangle*180/np.pi)
                    cos=np.cos(chain_angle-tiltangle)
                    #cos=np.cos(chain_angle-tiltangle)
                    scd_of_atoms = 0.5 * (3 * cos**2 - 1)
                    scd_of_tail += scd_of_atoms
                scd_of_tail=scd_of_tail/(len(tail)-1)
                totalscd=totalscd+scd_of_tail
            totalscd=totalscd/len(self.scd_tail_atoms_of[lipidmolecule])
            return (totalscd, tiltangle*180/np.pi)

        elif calculation_scheme=='New':
            for tail in self.scd_tail_atoms_of[lipidmolecule]:
                reallength=0
                startcoords=np.array([getcoords[''.join([lipidmolecule,time,str(res),tail[0]])]])
                endcoords=np.array([getcoords[''.join([lipidmolecule,time,str(res),tail[-1]])]])
                diffvector=endcoords-startcoords
                end_to_end_distance=calculate_distance(startcoords,endcoords)
                cos_tilt=np.dot(diffvector,[0,0,1])/np.linalg.norm(diffvector)
                tiltangle=np.arccos(cos_tilt)
                if tiltangle>np.pi/2:
                    tiltangle=abs(tiltangle-np.pi)
                scd_of_tail=0
                for atomindex in range(len(tail)-1): ### -1 because last res is not taken (taking index len() implies "range+1") 
                    atm1=tail[atomindex]    ### Attention: In the tail list only Scd-specific (every 2nd) atom is included!! Thus: atomindex-atomindex+1
                    atm2=tail[atomindex+1]
                    resid=str(res)
                    coordinates_of_atom1=np.array(getcoords[''.join([lipidmolecule,time,resid,atm1])])
                    coordinates_of_atom2=np.array(getcoords[''.join([lipidmolecule,time,resid,atm2])])
                    diffvector=coordinates_of_atom1-coordinates_of_atom2
                    normdiffvector=np.linalg.norm(diffvector)
                    cos=np.dot(diffvector,[0,0,1])/normdiffvector
                    reallength+=calculate_distance(coordinates_of_atom1,coordinates_of_atom2)
                    scd_of_atoms = 0.5 * (3 * cos**2 - 1)
                    scd_of_tail += scd_of_atoms
                scd_of_tail=scd_of_tail/(len(tail)-1)
                straightness=end_to_end_distance/reallength
                if straightness>=0.95:
                    neighbor_straightness=[]
                    neighbor_tilts=[]
                    neighbors=self.get_neighbor_of(res,time)
                    for neib in neighbors:
                        neibtype=self.resid_to_lipid[int(neib)]
                        scd,neibtilt,neib_straight=self.scd_of_res('Length_Based',neib,neibtype,time,getcoords)
                        if neib_straight>=0.95:
                            neighbor_straightness+=neib_straight
                            #xyangle=self.get_xy_angle(neib,neibtype,time,getcoords)
                        neighbor_tilts+=[neibtilt]
                    if len(neighbor_tilts)!=0: 
                        avg_tilt=sum(neighbor_tilts)/len(neighbor_tilts)
                        if 0.5*len(neighbor_straightness)<=len(neighbors) and abs(avg_tilt-tiltangle*180/np.pi)<=maxdiff:
                            #print('TILT!',res)
                            order_angle=np.arccos(((2*scd_of_tail+1)/3)**0.5)
                            #new_cos=(np.cos(order_angle-tiltangle)[0],np.cos(order_angle+tiltangle)[0])        
                            new_cos=np.cos(order_angle-avg_tilt*np.pi/180)[0]
                            scd_of_tail=0.5 * (3 * new_cos**2 - 1)    
                totalscd=totalscd+scd_of_tail
            totalscd=totalscd/len(self.scd_tail_atoms_of[lipidmolecule])
            return (totalscd, tiltangle*180/np.pi, straightness)        
        
        
              
    def create_scdfile(self,include_tilt='off',separate='on',maxdiff=0.0):           
        print("\n_____Extracting Scd values____\n\nTilt inclusion: {}\n".format(include_tilt)) 
        #resid2lipid=self.index_conversion_dict()[1] 
        getcoords,time=self.trajectory_to_gro()
        endtime=int(float(time))
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
            for t in range(0,endtime+1,self.dt):                       #the "time" variable should be the last declared time, thus the last read frame!
                time=''.join([str(t),'.0'])
                print("at time {} ...".format(t),end="\r") 
                sys.stdout.flush()
                for res in range(1,g_number_of_lipids+1):
                    lipidmolecule=self.resid_to_lipid[res]
                    totalscd=self.scd_of_res(include_tilt,res,lipidmolecule,time,getcoords,maxdiff)[0]
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
            for residue in range(1,g_number_of_lipids+1):
                print(". . . Working on residue: {} . . .".format(residue),end="\r")
                sys.stdout.flush()
                selectionfile=self.create_selectionfile_neighborsearch(residue)
                indexoutput=self.indexpath+'/neighbors_of_residue{}.ndx'.format(residue)
                datafileoutput=self.datapath+'/neighborfiles'+'/neighbors_of_residue{}.dat'.format(residue)
                if os.path.isfile(datafileoutput) and overwrite==False:
                    print("Neighbor file of residue {} already exists. Skipping.".format(residue))
                    pass
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
            firstline=ninfo.readline()
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
                firstline=infofile.readline()
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
                                
            
    ############################################################################################################################

    def calculate_energies(self,resindex_all,mdp_raw,overwrite=True,start=1,end=-1):
        neiblist=self.find_all_neighbors()
        if end==-1:
            end=self.number_of_lipids
        print('\n____Rerunning MD for energyfiles, yielding xvgtables with relevant energies.____\nCaution mdp-file must not have energy_grps indicated!\n')
        for i in range(start,end+1):
           print('\n',strftime("%H:%M:%S :", localtime()),'Working on lipid '+str(i)+'...')
           n_neibs=len(neiblist[i][1])
           if n_neibs % g_denominator == 0:
               number_of_groupfragments=(n_neibs//g_denominator)
           else:
               number_of_groupfragments=(n_neibs//g_denominator)+1

           print("Needing {} energy run(s)".format(number_of_groupfragments)) 
           
           for groupfragment in range(number_of_groupfragments):
                sys.stdout.flush()
                g_energy_output=self.energypath+'/xvgtables/energies_residue'+str(i)+'_'+str(groupfragment)+'.xvg'
                if os.path.isfile(g_energy_output) and overwrite==False:
                    print("Xvgtable for lipid {} part {} already exists. Will skip this calculation.".format(i,groupfragment))
                    continue
                    
                groupblockstart=groupfragment*g_denominator
                groupblockend=(groupfragment+1)*g_denominator
                energygroup_indeces=[neiblist[i][0]]+neiblist[i][1][groupblockstart:groupblockend]
                energygroup_strings=' '.join(["resid_"+str(x) for x in energygroup_indeces])
                select_energies_coulomb='\n'.join(["Coul-SR:resid_"+str(neiblist[i][0])+"-resid_"+str(x) for x in neiblist[i][1][groupblockstart:groupblockend]])
                select_energies_LJ='\n'.join(["LJ-SR:resid_"+str(neiblist[i][0])+"-resid_"+str(x) for x in neiblist[i][1][groupblockstart:groupblockend]])
                select_all_relevant_energies=select_energies_coulomb+"\n"+select_energies_LJ
                groupfragment=str(groupfragment) 


                #Create Mdpfile:
                try:
                    os.mkdir(self.energypath+'/mdpfiles')
                except OSError:
                    pass
                mdpfile=self.energypath+'/mdpfiles/energy_mdp_recalc_resid'+str(i)+'_'+groupfragment+'.mdp'
                with open(mdpfile,"w") as mdpfile_rerun, open(mdp_raw,"r") as mdpfile_raw:
                    mdp_raw_content=mdpfile_raw.readlines()
                    energygrpline='energygrps\t\t\t='+energygroup_strings+'\n'
                    mdp_raw_content.append(energygrpline)
                    mdpfile_rerun.write('\n'.join(mdp_raw_content)+'\n')
                
                
                #Create TPRFILE with GROMPP:
                print(strftime("%H:%M:%S :", localtime()),'...Creating .tpr-file...')
                try:
                    os.mkdir(self.energypath+'/tprfiles')
                except OSError:
                    pass
                tprfile_energyrerun=self.energypath+'/tprfiles/mdrerun_resid'+str(i)+'_'+groupfragment+'.tpr'
                mdpoutfile=self.energypath+'/mdpfiles'+'/mdrerun_resid'+str(i)+'_'+groupfragment+'.mdp'
                # grompp_arglist=[gmx_exec,'-nobackup','grompp','-f',mdpfile,'-p',self.toppath,'-c',self.gropath,'-o',tprfile_energyrerun,'-n',resindex_all,'-po',mdpoutfile]
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
                energyfile_output=self.energypath+'/edrfiles/energyfile_resid'+str(i)+'_'+groupfragment+'.edr'
                logoutput_file=self.energypath+'/logfiles'+'/mdrerun_resid'+str(i)+'.log'
                trajoutput="traj"+str(i)+'_'+groupfragment+'.trr'
                # mdrun_arglist=[gmx_exec,'-nobackup','mdrun','-s',tprfile_energyrerun,'-rerun',self.trjpath,'-e',energyfile_output,'-o',trajoutput,'-g',logoutput_file,'-nt','8']
                mdrun_arglist=[gmx_exec,'mdrun','-s',tprfile_energyrerun,'-rerun',self.trjpath,'-e',energyfile_output,'-o',trajoutput,'-g',logoutput_file,'-nt','8']
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
                g_energy_output=self.energypath+'/xvgtables/energies_residue'+str(i)+'_'+groupfragment+'.xvg'
                g_energy_arglist=[gmx_exec,'energy','-f',energyfile_output,'-s',tprfile_energyrerun,'-o',g_energy_output]
                inp_str=select_all_relevant_energies.encode()
                out,err=self.exec_gromacs(g_energy_arglist,inp_str)
                with open("gmx_energy.log","a") as logfile:
                    logfile.write(err.decode()); logfile.write(150*'_'); logfile.write(out.decode()); logfile.write(150*'_')
    

    def write_energyfile(self):
        print("______________Writing energy file____________\n")
        neiblist=self.find_all_neighbors()
        with open("all_energies","w") as energyoutput:
            energyoutput.write("{: <10}{: <10}{: <10}{: <20}{: <20}{: <20}\n".format("Time","Host","Neighbor","VdW","Coul","Etot"))
            for resid in range(1,g_number_of_lipids+1):
                n_neibs=len(neiblist[resid][1])
                if n_neibs % g_denominator == 0:
                    number_of_groupfragments=(n_neibs//g_denominator)
                else:
                    number_of_groupfragments=(n_neibs//g_denominator)+1
                for part in range(number_of_groupfragments):
                    print("Working on residue {}".format(resid),end="\r")
                    groupblockstart=part*g_denominator
                    groupblockend=(part+1)*g_denominator
                    neighbors_part_are=neiblist[resid][1][groupblockstart:groupblockend]
                    energypath=self.energypath+'/xvgtables/energies_residue'+str(resid)+'_'+str(part)+'.xvg'
                    with open(energypath,"r") as energyfile:
                        res_to_row={}
                        for energyline in energyfile: #folderlayout is: time Coul_resHost_resNeib LJ_resHost_resNeib ...
                            energyline_cols=energyline.split()
                            if '@ s' in energyline:                     #creating a dict to know which column(energies) belong to which residue
                                rownumber=int(energyline_cols[1][1:])+1                 #time is at row 0 !
                                resnumber=energyline_cols[3].split("resid_")[2][:-1]
                                energytype=energyline_cols[3].split("-")[0][1:]
                                res_to_row.update({energytype+resnumber:rownumber})
                            if '@' not in energyline and '#' not in energyline:     #pick correct energies from energyfile and print
                                time=float(energyline_cols[0])
                                for x in neighbors_part_are:
                                     vdw=energyline_cols[res_to_row['LJ'+str(x)]]
                                     coul=energyline_cols[res_to_row['Coul'+str(x)]]
                                     Etot=float(vdw)+float(coul)
                                     print("{: <10}{: <10}{: <10}{: <20}{: <20}{: <20.5f}".format(time,resid,x,vdw,coul,Etot), file=energyoutput)
       

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

    ############################################################################################################################
    
    
    
    def calculate_distance(self):
        print("\n____Calculating distances____\n")
        neiblist=self.find_all_neighbors()
        try:
            os.mkdir(self.datapath+'/distcalc')
        except OSError:
            pass
        for i in range(1,g_number_of_lipids+1):
            sys.stdout.flush()
            print("Working on lipid ",i,'...',end='\r')
            host=str(neiblist[i][0])
            neibors=[str(x) for x in neiblist[i][1]]
            selectstring='resid '+host+' and name P O3;\n'
            for x in neibors:
                selectstring+='resid '+x+' and name P O3;\n'
            selectionfile=self.temppath+'/seldist'
            outputfile=self.datapath+'/distcalc'+'/dist_to_hostresid'+str(i)+'.xvg'
            with open(selectionfile,"w") as sf:
                sf.write(selectstring)    
            g_pairdist_arglist=[gmx_exec,'-nobackup','pairdist','-f',self.trjpath,'-s',self.tprpath,'-sf',selectionfile,'-o',outputfile]
            out,err=self.exec_gromacs(g_pairdist_arglist)
            with open("gmx_pairdist.log","w") as logfile:
                logfile.write(err.decode()); logfile.write(150*'_'); logfile.write(out.decode()); logfile.write(150*'_')
        with open("all_distances.dat","w") as all_distances:
            # print("Time \t Host \t Neib \t distance/nm",file=all_distances)
            print("Time \t Host \t Neib \t distance/nm",file=all_distances)
            for resid in range(1,g_number_of_lipids+1):
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
                            neighbors_are=neiblist[resid][1]
                            for x in neighbors_are:
                                 distance=distanceline_cols[res_to_row[str(x)]]
                                 print("{} \t {} \t {} \t {}".format(time_distanceline,resid,x,distance),file=all_distances)

    def create_scdmap(self,scd_distribution,overwrite=False):
        grofile_output=self.temppath+'/calc_scd_for'+str(lipidmolecule)+'.gro'
        if os.path.isfile(datafileoutput) and overwrite==False:
            print("Writing trajectory as .gro")
            inp_str=str(lipidmolecule).encode()
            gmx_traj_arglist = [gmx_exec,'trjconv','-s',self.tprpath, '-f',self.trjpath,\
                                '-o',grofile_output,\
                                '-b', str(self.t_start), '-e', str(self.t_end),\
                                '-dt', str(self.dt),\
                                '-pbc', 'whole',\
                                ]
            out,err=self.exec_gromacs(gmx_traj_arglist,inp_str)      #Creates a gro file containing {timeframe:resid:atoms:xyz}
            with open("gmx_traj.log","w") as logfile:
                logfile.write(err.decode()); logfile.write(150*'_'); logfile.write(out.decode()); logfile.write(150*'_')
            with open(grofile_output,"r") as grofile:
                print(strftime("%H:%M:%S :", localtime()),"...read data output...")
                for line in grofile:
                    if 't=' in line:
                        time=str(float(line[line.index('t=')+2:].strip()))   #to get rid of obsolete decimals
                        print("...at time {}".format(time),end="\r")
                    for atom in list_of_needed_atoms:
                        if str(atom) in line:
                            resid=line[:5].strip()
                            lipidtype=line[5:9]
                            #resid2lipid.update({int(resid):lipidtype})
                            coordinates=[float(x) for x in line[20:44].split()]
                            keystring=lipidtype+time+resid+atom
            

