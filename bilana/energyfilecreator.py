''' All about creating energy files '''

#from bilana.systeminfo import mysystem
#global mysystem
from bilana import gromacstoolautomator, lipidmolecules

def create_Eofr_input(self,energyfile,distancefile):    # Ugly! Needs to be refactored.
    time_pair_to_E = {}
    time_pair_to_r = {}
    with open(energyfile, "r") as efile, open(distancefile, "r") as dfile:
        efile.readline()
        dfile.readline()
        for line in efile:
            cols = line.split()
            time = cols[0]
            respair = cols[1]+'-'+cols[2]
            Etot = cols[5]
            keystring = time+'_'+respair
            time_pair_to_E.update({keystring:Etot})
        for line in dfile:
            cols = line.split()
            time = cols[0]
            respair = cols[1]+'-'+cols[2]
            distance = cols[3]
            keystring = time+'_'+respair
            time_pair_to_r.update({keystring:distance})
    with open("Eofr_dat","w") as outfile:
        print(\
            '{: <10} {: <10} {: <10} {: < 10} {: <20}\n'\
            .format("Time", "Host", "Neib", "Distance", "Etot"), file=outfile)
        endtime = int(float(time))
        neiblist = self.get_neighbor_dict()
        for i in range(1, self.NUMBEROFPARTICLES+1):
            for t in range(self.t_start, endtime+1, self.dt):
                time = str(float(t))
                neighbors_are = neiblist[i][float(t)]
                for neib in neighbors_are:
                    respair = str(i)+'-'+str(neib)
                    Etot = time_pair_to_E[time+'_'+respair]
                    dist = time_pair_to_r[time+'_'+respair]
                    print("{: <10} {: <10} {: <10} {: < 10} {: <20.5f}".format(time, i, neib, dist, Etot), file=outfile)

def create_NofScd_input(self, scdfile, neighborfile):   # Ugly! Needs to be refactored.
    time_resid_to_scd = {}
    neighbors_of_host = {}
    hosts_without_neib = []
    with open(scdfile, "r") as sfile, open(neighborfile, "r") as nfile:
        sfile.readline()
        nfile.readline()
        for line in sfile:
            cols = line.split()
            time = cols[0]
            res = cols[1]
            scd = cols[3]
            keystring = time+'_'+res
            time_resid_to_scd.update({keystring:scd})
        for line in nfile:
            cols = line.split()
            time = str(float(cols[1]))
            host = cols[0]
            if int(cols[2]) == 0:
                #print(host, "has no neighbors at time", time)
                hosts_without_neib += [time+'_'+host]
                continue
            neighbors = cols[3]
            keystring = time+'_'+host
            neighbors_of_host.update({keystring:neighbors})
    with open("Nofscd_dppc.dat","w") as outfile_dppc,\
        open("Nofscd_chol.dat","w") as outfile_chol,\
        open("Nofscd_dupc.dat","w") as outfile_dupc:
        print(\
            '{: <10} {: <10} {: <20} {: <20}'\
            '{: <10} {: <10} {: <10}'\
            .format("Time", "Host", "Host_Scd", "N Neighbor",\
                       "N Chol", "N DPPC", "N DUPC"), file=outfile_dppc)
        print(\
            '{: <10} {: <10} {: <20} {: <20}'\
            '{: <10} {: <10} {: <10}'\
            .format("Time", "Host", "Host_Scd", "N Neighbor",\
                       "N Chol", "N DPPC", "N DUPC"), file=outfile_chol)
        print(\
            '{: <10} {: <10} {: <20} {: <20}'\
            '{: <10} {: <10} {: <10}'\
            .format("Time", "Host", "Host_Scd", "N Neighbor",\
                       "N Chol", "N DPPC", "N DUPC"), file=outfile_dupc)
        endtime = int(float(time))
        for i in range(1, self.NUMBEROFPARTICLES+1):
            restype = self.resid_to_lipid[i]
            #print("Working on residue {} ".format(i), end="\r")
            for t in range(self.t_start, endtime+1, self.dt):
                time = float(t)
                if time+'_'+str(i) in hosts_without_neib:
                    continue
                n_neibs = len(neighbors_of_host[time+'_'+str(i)].split(','))
                neibindexlist = neighbors_of_host[time+'_'+str(i)].split(',')
                neibtypelist = [self.resid_to_lipid[int(resid)] for resid in neibindexlist]
                nchol = neibtypelist.count('CHL1')
                ndppc = neibtypelist.count('DPPC')
                ndupc = neibtypelist.count('DUPC')
                scd_host = float(time_resid_to_scd[time+'_'+str(i)])
                if restype == 'DPPC':
                    print(\
                          '{: <10}{: <10}{: <20.5f}{: <20.5f}'\
                          '{: <10}{: <10}{: <10}'\
                          .format(time, i, scd_host, float(n_neibs),\
                                  nchol, ndppc, ndupc), file=outfile_dppc)
                elif restype == 'CHL1':
                    print(\
                          '{: <10}{: <10}{: <20.5f}{: <20.5f}'\
                          '{: <10}{: <10}{: <10}'\
                          .format(time, i, scd_host, float(n_neibs),\
                                  nchol, ndppc, ndupc), file=outfile_chol)
                elif restype == 'DUPC':
                    print(\
                          '{: <10}{: <10}{: <20.5f}{: <20.5f}'\
                          '{: <10}{: <10}{: <10}'\
                          .format(time, i, scd_host, float(n_neibs),\
                                  nchol, ndppc, ndupc), file=outfile_dupc)

class EofScd():
    ''' All about writing E(Scd)-files '''
    def __init__(self, systeminfo, parts, energyfilename, scdfilename):
        self.mysystem = systeminfo
        self.energyfilename = energyfilename
        self. scdfilename = scdfilename
        self.neiblist = gromacstoolautomator.Neighbors(systeminfo).get_neighbor_dict()
        if parts == 'complete':
            self.interactionskey = ['']
            self.parts = ''
        elif parts == 'head-tail':
            #interactions=['head-tail','head-head','tail-tail']
            self.interactionskey = ['h_h', 'h_t', 'h_w', 't_t', 't_w', 'w_w']
            self.parts = parts
        elif parts == 'head-tailhalfs':
            #interactions=['head-tail12','tail12-tail12','head-tail22','tail22-tail22']
            self.interactionskey = ['h_h', 'h_t12', 'h_t22', 'h_w',\
                                    't12_t12', 't12_t22', 't12_w',\
                                    't22_t22', 't22_w',\
                                    'w_w']
            self.parts = parts
        elif parts == 'carbons':
            self.interactionskey = []
            for i in range(lipidmolecules.shortestchain):
                for j in range(lipidmolecules.shortestchain):
                    if j > i:
                        continue
                    else:
                        self.interactionskey.append('C{0}_C{1}'.format(i,j))
            self.parts = parts
        self.lipidpairs = []
        for lipid1 in self.mysystem.molecules:   # Creates a list of pairs
            lipid1index = self.mysystem.molecules.index(lipid1) # Like 'DPPC_DPPC'
            for lipid2 in self.mysystem.molecules: 
                lipid2index = self.mysystem.molecules.index(lipid2)
                if lipid2index > lipid1index:  # DPPC_CHOL = CHOL_DPPC
                    break
                self.lipidpairs.append(''.join([lipid2, '_', lipid1]))

    def create_eofscdfile(self):
        print("______________Creating EofScd input file____________\n")
        timetoenergy, lasttime1 = self.read_energyinput(self.energyfilename)
        timetoscd, lasttime2 = self.read_scdinput(self.scdfilename)
        endtime = int(min(lasttime1, lasttime2))
        for pair in self.lipidpairs:
            self.write_outputfile(timetoenergy, timetoscd, endtime, pair)

    def write_outputfile(self, timetoenergy, timetoscd, endtime, wantedpair):
        outname = ''.join(['Eofscd', wantedpair, self.parts, '.dat'])
        with open(outname, "w") as outf:
            print(\
                  '{: ^10}{: ^10}{: ^15}{: ^8}{: ^10}'\
                  '{: ^15}{: ^15}{: ^15}'\
                  '{: ^15}{: ^15}'\
                  '{: ^15}{: ^5}'\
                  .format("Time", "Host", "Host_Scd", "Neib", "Neib_Scd",\
                            "Interaction", "DeltaScd", "AvgScd",\
                            "Etot", "Evdw",\
                            "Ecoul", "NChol"),\
                  file=outf)
            for host in range(1, self.mysystem.NUMBEROFMOLECULES+1):
                type_host = self.mysystem.resid_to_lipid[host]
                for t in range(self.mysystem.t_start, endtime+1, self.mysystem.dt):
                    #print("Working on residue {} at {}".format(host, t), end="\r")
                    t = float(t)
                    neighbors = self.neiblist[host][float(t)]
                    nchol = [self.mysystem.resid_to_lipid[neib] for neib in neighbors].count('CHL1')
                    for neib in neighbors:
                        type_neib = self.mysystem.resid_to_lipid[neib]
                        pair = (''.join([type_host, '_', type_neib]),\
                                ''.join([type_neib, '_', type_host]))
                        if neib < host or (pair[0] != wantedpair or pair[1] != wantedpair):
                            continue
                        respair = (host, neib)
                        #type_neib=self.resid_to_lipid[int(neib)]
                        #type_pair=type_host+'_'+type_neib
                        scd_host = timetoscd[(t, host)]
                        scd_neib = timetoscd[(t, neib)]
                        delta_scd = abs(scd_host-scd_neib)
                        avg_scd = (scd_host+scd_neib)/2
                        for inter in self.interactionskey:
                            if (t, respair, inter) in timetoenergy.keys():
                                Etot, VDW, COUL = timetoenergy[(t, respair, inter)]
                            else:
                                continue
                            print(\
                                  '{: <10}{: <10}{: <15.5f}{: <10}{: <15.5f}'\
                                  '{: <15}{: <15.5f}{: <15.5f}'\
                                  '{: <15.5f}{: <15.5f}'\
                                  '{: <15.5f}{: <5}'\
                                  .format(t, host, scd_host, neib, scd_neib,\
                                            inter, delta_scd, avg_scd,\
                                            float(Etot), float(VDW),\
                                            float(COUL), nchol),\
                                  file=outf)

    def read_energyinput(self, energyfile):
        timetoenergy = {}
        with open(energyfile, "r") as efile:
            efile.readline()
            for line in efile:
                #cols = [x.strip() for x in line.split(' ')]
                cols = line.split()
                if int(cols[1]) > int(cols[2]):
                    continue
                time = float(cols[0])
                respair = (int(cols[1]), int(cols[2]))
                #print(time, respair, end='\r')
                Etot = float(cols[6])
                VDW = float(cols[5])
                COUL = float(cols[4])
                inttype = cols[3]
                timetoenergy.update({(time, respair, inttype):(Etot, VDW, COUL)})
        return timetoenergy, time

    def read_scdinput(self, scdfile):
        timetoscd = {}
        with open(scdfile, "r") as sfile:
            sfile.readline()
            for line in sfile:
                cols = line.split()
                time = float(cols[0])
                res = int(cols[1])
                scd = float(cols[3])
                timetoscd.update({(time, res):scd})
        return timetoscd, time
###################################

