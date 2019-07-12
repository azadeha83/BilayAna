'''
    Module that control the creation of files containing N(S) data
'''
from .io import read_energyinput, read_neighborinput, read_scdinput
from .. import log
from ..analysis.neighbors import Neighbors
from ..definitions import lipidmolecules

LOGGER = log.LOGGER

class NofScd(Neighbors):
    '''
        This class handles the assembly of data from neighbor mapping (neighbor_info) and order parameter data (scd_distribution.dat)

    '''
    def __init__(self, inputfilename="inputfile"):
        super().__init__(inputfilename)
        self.components = self.molecules
    def create_NofScd_input(self,
        scdfilename="scd_distribution.dat",
        neighborfilename="neighbor_info"):
        time_to_scd, endtime1 = read_scdinput(scdfilename)
        time_to_neiblist, hosts_without_neib = read_neighborinput(neighborfilename)
        if neighborfilename != 'neighbor_info':
            nofscdname =  'Nofscd.dat_{}'.format(neighborfilename)
        else:
            nofscdname = "Nofscd.dat"
        with open(nofscdname, "w") as outfile:

            # Print header
            print(
                '{: <10}{: <10}{: <15}{: <20}{: <20}'\
                .format("Time", "Host", "Lipid_type", "Host_Scd", "Ntot")\
                +('{: ^10}'*len(self.components)).format(*self.components),
                file=outfile)

            for res in self.MOLRANGE:
                host_type = self.resid_to_lipid[res]
                if self.t_end < int(endtime1):
                    LOGGER.warning("Trajectory is longer (%s ps) than given end time (%s ps)",
                        endtime1, self.t_end)
                for t in range(self.t_start, self.t_end+1, self.dt):
                    time = float(t)
                    if t > int(endtime1):
                        break
                    if (time, res) in hosts_without_neib:
                        continue
                    neibindexlist = list(set(time_to_neiblist[(time, res)].split(',')))
                    n_neibs = len(neibindexlist)
                    neibtypelist = [self.resid_to_lipid[int(resid)] for resid in neibindexlist]
                    neib_comp_list = []
                    for lip in self.components:
                        number_neib  = neibtypelist.count(lip)
                        neib_comp_list.append(number_neib)
                    scd_host = float(time_to_scd[(time, res)])
                    print(
                          '{: <10}{: <10}{: <10}{: <20.5f}{: <20.5f}'\
                          .format(time, res, host_type, scd_host, float(n_neibs))\
                          +('{: ^10}'*len(neib_comp_list)).format(*neib_comp_list),
                          file=outfile)
