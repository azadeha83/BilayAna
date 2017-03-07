''' All about creating energy files '''

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


