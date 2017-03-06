import subprocess





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
       

   


class Energy():
    ''' All about energy calculation '''
    
class Scd():
    ''' All about scds '''
