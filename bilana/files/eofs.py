'''
    Module that control the creation of files containing E(S) data
'''
from .io import *
from .. import log
from ..analysis import neighbors
from ..analysis.neighbors import Neighbors
from ..definitions import lipidmolecules

LOGGER = log.LOGGER

class EofScd(Neighbors):
    '''
        Class to assemble necessary data to get the lipid interaction energy as a function of the lipids order parameter
            Files that are read:
                - File that contains data of every lipids order at each frame (scd_distribution.dat)
                - File that contains the interaction energy of all lipids at each frame
                - A mapping with neighbors (defined by our cutoff) of each lipid at each frame

        Optional a part of the lipid can be defined. Parts are:
            complete:       Interaction energy of two complete neighboring lipid molecules
            head-tail:      Interaction energy between the groups of part of lipids
            head-tailhalfs: Interaction energy where a distinction of the upper and lower part of the lipid chains is made
            carbons:        Interaction energy between the different carbon atoms of lipids. (A lot of data is produced here)
    '''
    def __init__(self, part, outputfile_tag='', inputfilename="inputfile", energyfilename="all_energies.dat", scdfilename="scd_distribution.dat", orientationfilename="orientation.dat"):
        super().__init__(inputfilename)
        self.energyfilename = energyfilename
        self.scdfilename = scdfilename
        self.orientationfilename = orientationfilename
        self.neiblist = neighbors.get_neighbor_dict()
        self.components = self.molecules
        self.outputfile_tag = outputfile_tag
        self.part = part

        # Define names for different interaction groups
        if part == 'complete':
            self.interactionskey = ['']
            self.part = ''
        elif part == 'head-tail':
            self.interactionskey = ['h_h', 'h_t', 'h_w', 't_t', 't_w', 'w_w']
        elif part == 'head-tailhalfs':
            self.interactionskey = ['h_h', 'h_t12', 'h_t22', 'h_w',\
                                    't12_t12', 't12_t22', 't12_w',\
                                    't22_t22', 't22_w',\
                                    'w_w']
        elif part == 'carbons':
            self.interactionskey = []
            for i in range(lipidmolecules.SHORTESTCHAIN):
                for j in range(lipidmolecules.SHORTESTCHAIN):
                    if j > i:
                        continue
                    else:
                        self.interactionskey.append('C{0}_C{1}'.format(i,j))

        self.lipidpairs = []
        for lipid1 in self.molecules:   # Creates a list of pairs
            lipid1index = self.molecules.index(lipid1) # Like 'DPPC_DPPC'
            for lipid2 in self.molecules:
                lipid2index = self.molecules.index(lipid2)
                if lipid2index > lipid1index:  # Ignore double entries: Lipid1_Lipid2 = Lipid2_Lipid1
                    break
                self.lipidpairs.append(''.join([lipid2, '_', lipid1]))        

    def create_eofscdfile(self):
        ''' Creates files of data of E(S) for each lipid part'''
        timetoscd, endtime = read_scdinput(self.scdfilename)
        timetoorientation = read_orientationinput(self.orientationfilename)
        for pair in self.lipidpairs:
            self.assemble_eofs(self.energyfilename, timetoscd, timetoorientation, pair, endtime)

    def create_selfinteractionfile(self):
        ''' Creates files of data of E(S) for each lipid '''
        timetoscd, endtime = read_scdinput(self.scdfilename)
        for lipid in self.molecules:
            self.assemble_eofs_selfinteraction(self.energyfilename, lipid, timetoscd, endtime)

    def assemble_eofs(self, energyfile, timetoscd, timetoorientation, lipidpair, endtime):
        ''' Read energyfile (all_energies.dat) in line and assemble from this data E(S) for a specific lipid pair '''
        components = self.components.copy() # dont change self.components

        if 'CHL1' in self.molecules:
            components += ["CHL1_host"] # Chols that are neighbor to host
            components += ["CHL1_neib"] # Chols ... to neib
            components += ['CHL1_both'] # Chols ... to both
            components += ['CHL1_partboth'] # Chols ... not to both

        if 'CHOL' in self.molecules:
            components += ["CHOL_host"] # Chols that are neighbor to host
            components += ["CHOL_neib"] # Chols ... to neib
            components += ['CHOL_both'] # Chols ... to both
            components += ['CHOL_partboth'] # Chols ... not to both
        
        if 'ERG' in self.molecules:
            components += ["ERG_host"] # Ergs that are neighbor to host
            components += ["ERG_neib"] # Ergs ... to neib
            components += ['ERG_both'] # Ergs ... to both
            components += ['ERG_partboth'] # Ergs ... not to both
        
        # Create Eofs output name
        if energyfile == 'all_energies.dat':
            outname = ''.join(['Eofscd', lipidpair, self.part, self.outputfile_tag, '.dat'])
        else:
            outname = ''.join(['Eofscd', lipidpair, self.part, self.outputfile_tag, '.dat'])

        LOGGER.debug("Will write to: %s", outname)
        LOGGER.debug("Opening energyfile: %s", energyfile )
        with open(energyfile, "r") as efile, open(outname, "w") as outf:
            # Print header
            print(\
                  '{: ^8}{: ^8}{: ^15}{: ^8}{: ^15}'\
                  '{: ^15}{: ^18}{: ^15}'\
                  '{: ^15}{: ^15}'\
                  '{: ^15}{: ^15}'\
                  .format("Time", "Host", "Host_Scd", "Neib", "Neib_Scd",
                            "Interaction", "DeltaScd", "AvgScd",
                            "Etot", "Evdw",
                            "Ecoul", "Ntot")\
                            + (len(components)*'{: <15}').format(*components) + '{: <15}{: <15}{: <15}'.format("Facing", "Notfacing", "Orient_flag"), file=outf)
            
            lipid1 = lipidpair.split('_')[0]
            lipid2 = lipidpair.split('_')[1]
            
            # Read all_energies file
            efile.readline()
            for line in efile:
                cols = line.split()
                
                if lipid1 == lipid2:

                    if int(cols[1]) > int(cols[2]):
                        continue

                time = float(cols[0])
                if time > self.t_end:
                    LOGGER.warning("Maybe did not use whole trajectory ( %s (real) vs %s (used) )",
                        endtime, self.t_end
                        )
                    continue
                if time < float(self.t_start) or time % self.dt != 0:
                    continue
                elif time > float(endtime):
                    continue
                
                # Get neighbor dictionary
                host, neib = int(cols[1]), int(cols[2])
                try:
                    neighbors = self.neiblist[host][float(time)]
                    neighbors_neib = self.neiblist[neib][float(time)]
                    

                except KeyError:
                    raise KeyError("Failed at host/neib: {}/{}\ttime: {}".format(host, neib, time))
               
                if neib not in neighbors: # in all_energies all neighbors that occur during simulation are stored
                    continue # skip if neib id not neighbor at this time
                
                type_host = self.resid_to_lipid[host]
                type_neib = self.resid_to_lipid[neib]
                type_pair = ('{0}_{1}'.format(type_host, type_neib), '{1}_{0}'.format(type_host, type_neib))
                
                if lipid1 == lipid2:
                    if neib < host or (type_pair[0] != lipidpair and type_pair[1] != lipidpair):  # Ignore second pair
                        # energies same for res1-res2 and res2-res1 so take just 1st version
                        continue
                else:
                    if type_host == type_neib or type_host != lipid1:
                        continue

                # Get energy values
                Etot = float(cols[6])
                COUL = float(cols[5])
                VDW  = float(cols[4])
                interactiontype = cols[3]

                # Get orientation data                
                # This part looks for all the sterol molecules in the neighbors of the host and neib (pair) and check if they are facing to them
                n_sterol_facing_host = 0
                n_sterol_notfacing_host = 0

                for n in list(set(neighbors)-set([neib])):
                    if self.resid_to_lipid[n] in lipidmolecules.STEROLS:                        
                        orientation_neibtohost = timetoorientation[(time, n, host)]
                        
                        if orientation_neibtohost == 1:
                            n_sterol_facing_host += 1
                        else:
                            n_sterol_notfacing_host += 1
                
                n_sterol_facing_neib = 0
                n_sterol_notfacing_neib = 0
                
                for n in list(set(neighbors_neib)-set([host])):
                    if self.resid_to_lipid[n] in lipidmolecules.STEROLS:

                        orientation_neibtoneib = timetoorientation[(time, n, neib)]

                        if orientation_neibtoneib == 1:
                            n_sterol_facing_neib += 1
                        else:
                            n_sterol_notfacing_neib += 1

                # This part assigns a flag to the orientation of the sterol-sterol or lipid-sterol with respect to each other:
                # 0: not facing, 1: one is facing to the other, 2: two are facing to each other
                orient_flag = 0

                if (lipidpair.split('_')[0] not in lipidmolecules.STEROLS) and (lipidpair.split('_')[1] in lipidmolecules.STEROLS):
                    
                    if type_host in lipidmolecules.STEROLS:

                        orientation_steroltohost = timetoorientation[(time, host, neib)]
                    else:
                        orientation_steroltohost = timetoorientation[(time, neib, host)]

                    if orientation_steroltohost == 0:
                        orient_flag = 0
                    else:
                        orient_flag = 1
                
                elif (lipidpair.split('_')[0] in lipidmolecules.STEROLS) and (lipidpair.split('_')[1] in lipidmolecules.STEROLS):
                    
                    orientation_steroltosterol1 = timetoorientation[(time, host, neib)]
                    orientation_steroltosterol2 = timetoorientation[(time, neib, host)]
                    
                    if orientation_steroltosterol1 == 0 and orientation_steroltosterol2 == 0:
                        orient_flag = 0
                    elif orientation_steroltosterol1 == 1 and orientation_steroltosterol2 == 1:
                        orient_flag = 2
                    else:
                        orient_flag = 1
                
                sterol_comp_list = []
                sterol_comp_list.append(orient_flag)

                # Get number of neighbors of type
                pair_neibs = list(set(neighbors+neighbors_neib)-set([host])-set([neib]))
        
                ntot = len(pair_neibs)
                neib_comp_list = []
                for lip in self.components:
                    ncomp = [self.resid_to_lipid[N] for N in pair_neibs].count(lip)
                    neib_comp_list.append(ncomp)
                
                # Get number of cholesterol neighbors, that are neighbor to both PLs
                if 'CHL1' in self.molecules:
                    host_chol = [self.resid_to_lipid[N] for N in neighbors if N != neib].count('CHL1')
                    neib_comp_list.append(host_chol)

                    hostneib_chol = [self.resid_to_lipid[N] for N in neighbors_neib].count('CHL1')
                    neib_comp_list.append(hostneib_chol)

                    shared_neighbors = list(set(neighbors) & set(neighbors_neib))
                    shared_chol = [self.resid_to_lipid[N]\
                        for N in shared_neighbors].count('CHL1')
                    
                    not_shared_neighbors = list(set(neighbors) ^ set(neighbors_neib))
                    not_shared_chol = [self.resid_to_lipid[N]\
                        for N in not_shared_neighbors].count('CHL1')
                    
                    neib_comp_list.append(shared_chol)
                    neib_comp_list.append(shared_chol + not_shared_chol/2)
                
                if 'CHOL' in self.molecules:
                    host_chol = [self.resid_to_lipid[N] for N in neighbors if N != neib].count('CHOL')
                    neib_comp_list.append(host_chol)

                    hostneib_chol = [self.resid_to_lipid[N] for N in neighbors_neib].count('CHOL')
                    neib_comp_list.append(hostneib_chol)

                    shared_neighbors = list(set(neighbors) & set(neighbors_neib))
                    shared_chol = [self.resid_to_lipid[N]\
                        for N in shared_neighbors].count('CHOL')
                    
                    not_shared_neighbors = list(set(neighbors) ^ set(neighbors_neib))
                    not_shared_chol = [self.resid_to_lipid[N]\
                        for N in not_shared_neighbors].count('CHOL')
                    
                if 'ERG' in self.molecules:
                    host_erg = [self.resid_to_lipid[N] for N in neighbors if N != neib].count('ERG')
                    neib_comp_list.append(host_erg)

                    hostneib_erg = [self.resid_to_lipid[N] for N in neighbors_neib].count('ERG')
                    neib_comp_list.append(hostneib_erg)

                    shared_neighbors = list(set(neighbors) & set(neighbors_neib))
                    shared_erg = [self.resid_to_lipid[N]\
                        for N in shared_neighbors].count('ERG')
                    
                    not_shared_neighbors = list(set(neighbors) ^ set(neighbors_neib))
                    not_shared_erg = [self.resid_to_lipid[N]\
                        for N in not_shared_neighbors].count('ERG')

                    neib_comp_list.append(shared_erg)
                    neib_comp_list.append(shared_erg + not_shared_erg/2)
                
                # Get order values
                scd_host = timetoscd[(time, host)]
                scd_neib = timetoscd[(time, neib)]
                delta_scd = abs(scd_host-scd_neib)
                avg_scd = (scd_host+scd_neib)/2

                # Write output line
                print(
                      '{: <10}{: <10}{: <15.5f}{: <10}{: <15.5f}'
                      '{: <15}{: <15.5f}{: <15.5f}'
                      '{: <15.5f}{: <15.5f}'
                      '{: <15.5f}{: <15}'\
                      .format(time, host, scd_host, neib, scd_neib,
                              interactiontype, delta_scd, avg_scd,
                              float(Etot), float(VDW),
                              float(COUL), ntot)\
                      + (len(neib_comp_list)*'{: <15}').format(*neib_comp_list)+'{: <15}{: <15}{: <15}'.format(n_sterol_facing_host+n_sterol_facing_neib,
                      n_sterol_notfacing_host+n_sterol_notfacing_neib,orient_flag),file=outf)

    def assemble_eofs_selfinteraction(self, energyfile, lipid, timetoscd, endtime):
        ''' Same as assemble eofs  '''
        outname = ''.join(['Eofscd_self', lipid, self.part, '.dat'])
        with open(energyfile, "r") as efile, open(outname, "w") as outf:
            print(\
                  '{: ^8}{: ^8}{: ^15}'
                  '{: ^15}{: ^15}{: ^15}'
                  '{: ^15}{: ^15}{: ^15}{: ^15}{: ^10}{: ^10}{: ^10}'\
                  .format("Time", "Host", "Host_Scd",
                            "Etot",  "Evdw_tot", "Ecoul_tot",
                            "Evdw_SR", "Evdw_14", "Ecoul_SR", "Ecoul_14", "Ntot", "NChol", "NDPPC",
                            ),
                  file=outf)
            efile.readline()
            for line in efile:
                #cols = [x.strip() for x in line.split(' ')]
                cols = line.split()
                host = int(cols[1])
                type_host = self.resid_to_lipid[host]
                time = float(cols[0])
                Etot = float(cols[2])
                coul_sr = float(cols[4])
                vdw_sr = float(cols[3])
                coul_14 = float(cols[6])
                vdw_14 = float(cols[5])
                coul_tot = float(cols[8])
                vdw_tot = float(cols[7])
                if type_host != lipid:
                    continue
                if time < float(self.t_start) or time % self.dt != 0:
                    continue
                elif time > float(endtime):
                    continue
                neibs = self.neiblist[host][time]
                nchol = [self.resid_to_lipid[N] for N in neibs].count('CHL1')
                ndppc = [self.resid_to_lipid[N] for N in neibs].count('DPPC')
                ntot = len(neibs)
                scd_host = timetoscd[(time, host)]
                print(\
                      '{: <10}{: <10}{: <15.5f}'
                      '{: <15.5f}{: <15.5f}{: <15.5f}'
                      '{: <15.5f}{: <15.5f}{: <15.5f}{: <15.5f}{: <10}{: <10}{: <10}'
                      .format(time, host, scd_host,
                                Etot, vdw_tot, coul_tot,
                                vdw_sr, vdw_14, coul_sr, coul_14, ntot, nchol, ndppc),
                      file=outf)
