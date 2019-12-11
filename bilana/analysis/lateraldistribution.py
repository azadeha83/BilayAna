import numpy as np
import MDAnalysis as mda
from . import neighbors
from . import protein
from ..import log
from ..import common as cm

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

#def calc_averagedistance(self, distancefile):
#    distdata={}
#    neiblist=self.find_all_neighbors()
#    with open(distancefile,"r") as df:
#        df.readline()
#        for line in df:
#            col=line.split()
#            time=float(col[0])
#            if float(time)<float(self.t_start):
#                continue
#            host=int(col[1])
#            print(host,time,end='\r')
#            neib=int(col[2])
#            if self.resid_to_lipid[host]!='DPPC' and self.resid_to_lipid[neib]!='DPPC':
#                continue
#            distance=float(col[3])
#            try:
#                neibs=neiblist[host][time]
#                nchol=[self.resid_to_lipid[int(i)] for i in neibs].count('CHL1')
#            except KeyError:
#                print('Attention: Neighbor file seems not to be complete time "{}" is missing'.format(time))
#            try:
#                distdata[nchol]+=[distance]
#            except KeyError:
#                distdata.update({nchol:[distance]})
#    with open("average_distance_DPPC.dat",'w') as outputf:
#        for key in distdata:
#            avg=sum(distdata[key])/len(distdata[key])
#            print("{0: <5} {1: <20}".format(key,avg),file=outputf)

def create_prot_pl_distance(neighbor_inst, outputfilename="protein_distances.dat"):
    ''' Creates file with columns
            < Time > < Residue > < distance > 
    '''
    u = neighbor_inst.universe
    refatoms = neighbor_inst.reference_atom_selection
    traj_len = len(neighbor_inst.universe.trajectory)
    leaflets = [[[],], [[],]]

    with open(outputfilename, "w") as outf:
        print("{: <20}{: <20}{: <20}{: <20}".format("Time", "Residue", "distance", "leaflet"), file=outf)
        for t in range(traj_len):
            time = u.trajectory[t].time
            if neighbor_inst.t_end < time or neighbor_inst.t_start > time:
                continue

            refatomgrp = u.select_atoms(refatoms)
            LOGGER.debug("Found %s atoms: %s ", len(refatomgrp.atoms), refatomgrp.atoms)
            leaflets_tmp = leaflets
            leaflets = neighbor_inst.get_ref_positions("atom", refatomgrp) # leaflets=( [(resid1, pos1), ...], [(residN, posN), ...] )

            # Get correct protein position
            ref_res_prot = protein.get_reference_resids() # [pos_leaf1, pos,leaf2]
            sel1 = "resid {}".format(' '.join([str(i) for i in ref_res_prot[0]]))
            sel2 = "resid {}".format(' '.join([str(i) for i in ref_res_prot[1]]))
            prot_pos_leaf1 = u.select_atoms( sel1 ).center_of_mass()
            prot_pos_leaf2 = u.select_atoms( sel2 ).center_of_mass()
            leaflets = [ [(0, prot_pos_leaf1), *leaflets[0]], [(0, prot_pos_leaf2), *leaflets[1]] ]

            if len(leaflets[0]) != len(leaflets[1]):
                #raise ValueError("Number of found ref positions differs: {} vs {}".format(len(leaflets[0]), len(leaflets[1])))
                LOGGER.warning("frame %s: Number of found ref positions differs: %s vs %s", t, len(leaflets[0]), len(leaflets[1]))
                LOGGER.warning("Difference1: %s", set([i[0] for i in leaflets[0]]).symmetric_difference(set([i[0] for i in leaflets_tmp[0]])))

            boxdim = u.dimensions.copy()
            cutoff = neighbor_inst.cutoff

            for leafndx, all_coords_per_leaflet in enumerate(leaflets):
                if leafndx == 0:
                    leafndx = 1
                elif leafndx == 1:
                    leafndx = 0


                for i in range(len(all_coords_per_leaflet)): # Make it "2D"
                    all_coords_per_leaflet[i][1][2] = 0      # By setting z values to 0

                prot_res, prot_pos = all_coords_per_leaflet[0]
                position_array =  np.array([pos for resid, pos in all_coords_per_leaflet]) # Get all positions in leaflet selection
                dist_array = mda.lib.distances.distance_array(prot_pos, position_array, box=boxdim)[0] # output is [ [[dist1], [dist2], ...] ]

                neiblist = []
                for resndx, distance in enumerate(dist_array):
                    LOGGER.debug("prot %s, resndx %s, neib %s", prot_res, resndx, all_coords_per_leaflet[resndx][0])
                    LOGGER.debug("protpos vs neibpos: %s vs %s", prot_pos, all_coords_per_leaflet[resndx][1] )
                    LOGGER.debug("box %s", boxdim)
                    #LOGGER.debug("dist %s and dist_self %s", distance, dist_helper(prot_pos,  all_coords_per_leaflet[resndx][1], boxdim))
                    resid_neib = all_coords_per_leaflet[resndx][0]
                    if resid_neib == 0:
                        continue
                    neiblist.append(resid_neib)
                    print("{: <20}{: <20}{: <20}{: <20}".format(time, resid_neib, distance/10.0, leafndx), file=outf)
                if len(neiblist) != len(list(set(neiblist))):
                    print( "ATATATATATTENTION", neiblist )


def write_neighbortype_distr(neighbor_instance, fname="neighbortype_distribution.dat"):
    '''
        Creates datafile < fname > with columns:
        < time >  < residue > < resname > < N comp 1 > < N comp 2 > ... < N comp i >
    '''
    neiblist = neighbors.get_neighbor_dict()
    components = neighbor_instance.molecules
    with open(fname, "w") as outf:
        print("{: <12}{: <10}{: <7}".format("time", "resid", "lipidtype")\
            + (len(components)*'{: ^7}').format(*components),
            file=outf)
        for resid in neighbor_instance.MOLRANGE:
            LOGGER.debug("At res %s", resid)
            lipidtype = neighbor_instance.resid_to_lipid[resid]
            for time in range(neighbor_instance.t_start, neighbor_instance.t_end + neighbor_instance.dt, neighbor_instance.dt):
                LOGGER.debug("At time %s", time)
                neibs = neiblist[resid][float(time)]
                neib_comp_list = []
                for lip in components:
                    ncomp = [neighbor_instance.resid_to_lipid[N] for N in neibs].count(lip)
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

def separate_quadrants(coords, vec, ang_to_x):
    '''
        Separates coords into 4 quadrants like
        vec needs to have the same root!!!
        _______
        |  vec|
        |\ |0/|
        |1\|/2|
        | /|\ |
        |/ |3\|
        -------        
    '''
    def is_within(p, vecleft, vecright):
        ''' Checks in which area spanned by AxB point p lies 
            Both vectors must face upwards!!!
        '''
        #         Upper    Left    Right  Lower
        areas = [(1, -1), (-1, -1), (1, 1), (-1, 1)]
        det1 = np.linalg.det(np.matrix([p, vecleft]))
        det2 = np.linalg.det(np.matrix([p, vecright]))
        areandx = areas.index( (np.sign(det1), np.sign(det2)) )
        return areandx
    areas = []
    vecleft  = cm.rotate_2d(vec, 2*np.pi - ang_to_x)
    vecright  = cm.rotate_2d(vec, ang_to_x)
    for p in coords:
        area = is_within(p, vecleft, vecright)
        areas.append(area)
    return np.array(areas)



def create_scd_map(sysinfo, orderfilename="scd_distribution.dat"):
    ''' 
        Creates map with order parameter entries and lipid + protein relative positions
        protein_pos phi r S
                      /r
                     /phi
                    /______
                  TMD
                
    '''

