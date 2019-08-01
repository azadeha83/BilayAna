
import numpy as np
import MDAnalysis as mda
from . import neighbors
from ..definitions import lipidmolecules
from .. import log

LOGGER = log.LOGGER

def get_neighbor_of(hostres, time):
    'returns list of resids of neighbors of host resid; Inputs hostres: Resid host, time: Time as integer'
    time=float(time)
    with open("neighbor_info","r") as ninfo:
        ninfo.readline()
        for line in ninfo:
            cols = line.split()
            res = str(cols[0])
            t = cols[1]
            if res == str(hostres) and t == ''.join([str(time),'00']):
                try:
                    return cols[3].split(',')
                except IndexError:
                    return []
    print(time, "I should never get here...")

def calc_averagedistance(self, distancefile):
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


def write_neighbortype_distr(systeminfo, fname="neighbortype_distribution.dat"):
    '''
        Creates datafile < fname > with columns:
        < time >  < residue > < resname > < N comp 1 > < N comp 2 > ... < N comp i >
    '''
    neiblist = neighbors.get_neighbor_dict()
    components = systeminfo.molecules
    with open(fname, "w") as outf:
        print("{: <12}{: <10}{: <7}".format("time", "resid", "lipidtype")\
            + (len(components)*'{: ^7}').format(*components),
            file=outf)
        for resid in systeminfo.MOLRANGE:
            lipidtype = systeminfo.resid_to_lipid[resid]
            for time in range(systeminfo.t_start, systeminfo.t_end, systeminfo.dt):
                neibs = neiblist[resid][float(time)]
                neib_comp_list = []
                for lip in components:
                    ncomp = [systeminfo.resid_to_lipid[N] for N in neibs].count(lip)
                    neib_comp_list.append(ncomp)
                print("{: <12}{: <10}{: <7}".format(time, resid, lipidtype)\
                    + (len(neib_comp_list)*'{: ^7}').format(*neib_comp_list),
                    file=outf)

def get_distance(systeminfo, refstr, selstr, outputname=None, dim=2):
    ''' '''
    outformat_str = "{: <20}{: <15}{: <15}{: >20}{: >15}{: >15}{: >15}"
    outformat = "{: <20}{: <15}{: <15}{: >20}{: >15.5f}{: >15.5f}{: >15.5f}"
    u = systeminfo.universe
    len_traj = len(u.trajectory)
    if outputname is None:
        outputname = "distance_ref{}_sel{}_dim{}.dat".format(refstr, selstr, dim).replace(" ", "_")
    with open(outputname, "w") as outf:
        print(outformat_str.format("time", "resid", "resname", "distance", "x", "y", "z"), file=outf)
        for t in range(len_traj):
            time = u.trajectory[t].time
            if systeminfo.t_end < time or systeminfo.t_start > time:
                continue
            LOGGER.info("Time %s", time)
            sel     = u.select_atoms(selstr)
            pos_ref = u.select_atoms(refstr).center_of_mass()
            pos_sel = sel.positions
            box     = u.dimensions
            if dim == 2:
                pos_ref[2] = 0
                pos_sel[:,2] = 0
                box[2] = 1
            elif dim == 3:
                pass
            else:
                raise ValueError("Dimension not possible. Use either 2 or 3")
            LOGGER.debug("Sel %s", sel,)
            LOGGER.debug("Resids %s", sel.resids)
            pos_array = mda.lib.distances.distance_array(pos_ref, pos_sel, box=box)[0]
            LOGGER.debug("Array: %s", pos_array)
            for ind, res in enumerate(sel.resids):
                resname = sel[ind].resname
                pos = sel[ind].position
                dist = pos_array[ind]
                print(outformat.format(time, res, resname, dist, *pos), file=outf)
