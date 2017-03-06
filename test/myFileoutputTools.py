




def create_Eofr_input(self,energyfile,distancefile):
        time_pair_to_E={}
        time_pair_to_r={}
        with open(energyfile,"r") as efile, open(distancefile,"r") as dfile:
            firstline_efile=efile.readline()
            firstline_dfile=dfile.readline()
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
            for i in range(1,g_number_of_lipids+1):
                for t in range(0,endtime+1,self.dt):
                    time=str(t)+'.0'
                    neighbors_are=neiblist[i][1]
                    for neib in neighbors_are:
                        respair=str(i)+'-'+str(neib)
                        Etot=time_pair_to_E[time+'_'+respair]
                        dist=time_pair_to_r[time+'_'+respair]
                        print("{: <10} {: <10} {: <10} {: < 10} {: <20.5f}".format(time,i,neib,dist,Etot),file=outfile)

    def create_EofScd_input(self,energyfile,scdfile,neighborfile):
        time_pair_to_E={}
        time_resid_to_scd={}
        neighbors_of_host={}
        hosts_without_neib=[]
        with open(energyfile,"r") as efile, open(scdfile,"r") as sfile, open(neighborfile,"r") as nfile:
            firstline_efile=efile.readline()
            firstline_sfile=sfile.readline()
            firstline_nfile=nfile.readline()
            for line in efile:
                cols=line.split()
                time=cols[0]
                respair=cols[1]+'-'+cols[2]
                Etot=cols[5]
                VDW=cols[4]
                COUL=cols[3]
                keystringTOT="TOT"+time+'_'+respair
                keystringVDW="VDW"+time+'_'+respair
                keystringCOUL="COUL"+time+'_'+respair
                time_pair_to_E.update({keystringTOT:Etot})
                time_pair_to_E.update({keystringVDW:VDW})
                time_pair_to_E.update({keystringCOUL:COUL})
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
        with open("Eofscd_dppc_dppc.dat","w") as outfile_dppc_dppc, open("Eofscd_dppc_chol.dat","w") as outfile_dppc_chol,open("Eofscd_chol_chol.dat","w") as outfile_chol_chol, open("Eofscd_dupc_chol.dat","w") as outfile_dupc_chol:
            if 'DPPC' in self.molecules:
                print("{: <10}{: <10}{: <20}{: <10}{: <20}{: <20}{: <20}{: <20}{: <20}{: <20}".format("Time","Host","Host_Scd","Neib","Neib_Scd","DeltaScd","AvgScd","Etot","Evdw","Ecoul"),file=outfile_dppc_dppc)
            if 'CHL1' in self.molecules and 'DPPC' in self.molecules: 
                print("{: <10}{: <10}{: <20}{: <10}{: <20}{: <20}{: <20}{: <20}{: <20}{: <20}".format("Time","Host","Host_Scd","Neib","Neib_Scd","DeltaScd","AvgScd","Etot","Evdw","Ecoul"),file=outfile_dppc_chol)
            if 'CHL1' in self.molecules and 'DUPC' in self.molecules: 
                print("{: <10}{: <10}{: <20}{: <10}{: <20}{: <20}{: <20}{: <20}{: <20}{: <20}".format("Time","Host","Host_Scd","Neib","Neib_Scd","DeltaScd","AvgScd","Etot","Evdw","Ecoul"),file=outfile_dupc_chol)
            if 'CHL1' in self.molecules: 
                print("{: <10}{: <10}{: <20}{: <10}{: <20}{: <20}{: <20}{: <20}{: <20}{: <20}".format("Time","Host","Host_Scd","Neib","Neib_Scd","DeltaScd","AvgScd","Etot","Evdw","Ecoul"),file=outfile_chol_chol)
            endtime=int(float(time))
            for i in range(1,g_number_of_lipids+1):
                print("Working on residue {} ".format(i),end="\r")
                type_host=self.resid_to_lipid[i]
                for t in range(0,endtime+1,self.dt):
                    time=str(t)+'.0'
                    if time+'_'+str(i) in hosts_without_neib:
                        continue
                    neighbors_are=neighbors_of_host[time+'_'+str(i)].split(',')
                    for neib in neighbors_are:
                        type_neib=self.resid_to_lipid[int(neib)]
                        type_pair=type_host+'_'+type_neib
                        scd_host=float(time_resid_to_scd[time+'_'+str(i)])
                        scd_neib=float(time_resid_to_scd[time+'_'+str(neib)])
                        delta_scd=abs(scd_host-scd_neib)
                        avg_scd=(scd_host+scd_neib)/2
                        respair=str(i)+'-'+str(neib)
                        keystrTOT="TOT"+time+'_'+respair
                        keystrVDW="VDW"+time+'_'+respair
                        keystrCOUL="COUL"+time+'_'+respair
                        Etot=time_pair_to_E[keystrTOT]
                        VDW=time_pair_to_E[keystrVDW]
                        COUL=time_pair_to_E[keystrCOUL]
                        if type_pair=='DPPC_DPPC':
                            print("{: <10}{: <10}{: <20.5f}{: <10}{: <20.5f}{: <20.5f}{: <20.5f}{: <20.5f}{: <20.5f}{: <20.5f}".format(time,i,scd_host,neib,scd_neib,delta_scd,avg_scd,float(Etot),float(VDW),float(COUL)),file=outfile_dppc_dppc)
                        if type_pair=='DPPC_CHL1' or type_pair=='CHL1_DPPC':
                            print("{: <10}{: <10}{: <20.5f}{: <10}{: <20.5f}{: <20.5f}{: <20.5f}{: <20.5f}{: <20.5f}{: <20.5f}".format(time,i,scd_host,neib,scd_neib,delta_scd,avg_scd,float(Etot),float(VDW),float(COUL)),file=outfile_dppc_chol)
                        if type_pair=='DUPC_CHL1':
                            print("{: <10}{: <10}{: <20.5f}{: <10}{: <20.5f}{: <20.5f}{: <20.5f}{: <20.5f}{: <20.5f}{: <20.5f}".format(time,i,scd_host,neib,scd_neib,delta_scd,avg_scd,float(Etot),float(VDW),float(COUL)),file=outfile_dupc_chol)
                        if type_pair=='CHL1_CHL1':
                            print("{: <10}{: <10}{: <20.5f}{: <10}{: <20.5f}{: <20.5f}{: <20.5f}{: <20.5f}{: <20.5f}{: <20.5f}".format(time,i,scd_host,neib,scd_neib,delta_scd,avg_scd,float(Etot),float(VDW),float(COUL)),file=outfile_chol_chol)

    def create_NofScd_input(self,scdfile,neighborfile):
        time_resid_to_scd={}
        neighbors_of_host={}
        hosts_without_neib=[]
        with open(scdfile,"r") as sfile, open(neighborfile,"r") as nfile:
            firstline_sfile=sfile.readline()
            firstline_nfile=nfile.readline()
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
            for i in range(1,g_number_of_lipids+1):
                restype=self.resid_to_lipid[i]
                print("Working on residue {} ".format(i),end="\r")
                for t in range(0,endtime+1,self.dt):
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
            firstline_sfile=sfile.readline()
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
                for i in range(1,g_number_of_lipids+1):
                    restype=self.resid_to_lipid[i]
                    if restype!=lipid:
                        continue
                    print("Working on residue {} ".format(i),end="\r")
                    for t in range(0,endtime+1,self.dt):
                        time=str(t)+'.0'
                        scd_host=float(time_resid_to_scd[time+'_'+str(i)])
                        print("{: <10}{: <10}{: <20.5f}".format(time,i,scd_host),file=outfile)



        



    def create_NLip_by_Ntot(self,neighborfile):
        lipid_pairs={}
        print(self.molecules)
        mol=self.molecules.copy()
        mol.remove('CHOL')
        for lipid in mol:
            for index in range(len(mol)):
                lipid_pairs.update({str(lipid)+'_'+str(mol[index]):{}})
        with open(neighborfile,"r") as neibfile:
            header=neibfile.readline()
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
                    except:
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

    def create_z_vs_(self,var):
        getcoords=self.trajectory_to_gro(overwrite)
        pass
            
            
    def H_real_vs_H_fit(self):
        pass


