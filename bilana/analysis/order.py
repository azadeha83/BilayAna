'''
    This module focuses on the analysis of structural features of lipids in a bilayer
        - The order parameter S is calculate using class Order
        - The tilt of lipids is calculated in calc_avg_tilt

'''
import numpy as np
import pandas as pd
import MDAnalysis as mda

from . import neighbors
from .neighbors import Neighbors
from .. import log
from ..systeminfo import SysInfo
from ..definitions import lipidmolecules
from ..definitions.structure_formats import REGEXP_GRO

LOGGER = log.LOGGER

class Order(Neighbors):
    '''
        This class handles the calculation of the order of lipid chains
            - Create order distribution file using create_orderfile
    '''
    LOGGER = LOGGER

    def __init__(self, inputfilename="inputfile"):
        super().__init__(inputfilename)
        self.atomlist = lipidmolecules.scd_tail_atoms_of
        self.neiblist = neighbors.get_neighbor_dict()
        self.components = self.molecules

    def create_orderfile(self, mode="CC", outputfile='scd_distribution.dat', with_tilt_correction="tilt.csv"):
        '''
            Calculate order parameter S = 0.5 (3cos^2(a)-1) of lipids.
            a is the angle between the bilayer normal and (set with mode):
                CC -- Vector between Cn-Cn+2 carbon atoms of chain. Averaged over all the angles

            if with tilt_correction avg tilt angle per time is read from name given in variable
        '''
        if with_tilt_correction:
            try:
                dat = pd.read_csv(with_tilt_correction)
            except FileNotFoundError:
                calc_tilt(self)
                dat = pd.read_csv(with_tilt_correction)
        if mode == "CC":
            calc_s = self.scc_of_res
        else:
            raise ValueError("Mode not yet implement use mode=CC (default)")
        with open(outputfile, "w") as scdfile:
            components = self.molecules
            print("{: <12}{: <10}{: <7}{: <15}".format("Time", "Residue", "Type", "Scd")\
                + (len(components)*'{: ^7}').format(*components),
                file=scdfile)
            len_traj = len(self.universe.trajectory)
            for t in range(len_traj):
                time = self.universe.trajectory[t].time
                if time % self.dt != 0 or self.t_start > time or self.t_end < time:
                    continue
                LOGGER.info("At time %s", time)
                for res in self.MOLRANGE:
                    LOGGER.debug("At time %s and residue %s", time, res)
                    resname = self.resid_to_lipid[res]
                    LOGGER.debug("resname %s", resname)
                    if resname[:-2] not in lipidmolecules.TAIL_ATOMS_OF.keys()\
                        and resname not in lipidmolecules.STEROLS\
                        and resname not in lipidmolecules.PROTEINS:
                        LOGGER.debug("Skipping")
                        continue
                    leaflet = self.res_to_leaflet[res]
                    if with_tilt_correction:
                        ang_corr = np.array(dat.loc[(dat.time == time) & (dat.leaflet == leaflet)][["x", "y", "z"]])[0]
                        LOGGER.debug("Corrected angle %s", ang_corr)
                    else:
                        ang_corr = None
                    scd_value = calc_s(self.universe, res, tilt_correction=ang_corr)
                    neibs = self.neiblist[res][float(time)]
                    neib_comp_list = []
                    LOGGER.debug("Calculated order: %s", scd_value)
                    for lip in self.components:
                        ncomp = [self.resid_to_lipid[N] for N in neibs].count(lip)
                        neib_comp_list.append(ncomp)
                    print("{: <12.2f}{: <10}{: <7}{: <15.8}".format(
                        time, res, resname, scd_value)\
                        + (len(neib_comp_list)*'{: ^7}').format(*neib_comp_list),
                        file=scdfile)

    def create_orientationfile(self, outputfile="orientation.dat"):
        ''' Assigning orientation of the sterol molecule with respect to the neighboring molecules '''

        u = self.universe
        len_traj = len(u.trajectory)

        with open(outputfile, "w") as orifile:
            print("{: <12}{: <10}{: <10}{: <15}{: <15}{: <15}".format("Time", "Resid", "Resname", "Neib_resid", "Neib_resn", "Orientation")\
                ,file=orifile)
            
            for t in range(len_traj):
                time = u.trajectory[t].time 
                if time % self.dt != 0 or time > self.t_end or time < self.t_start:
                    continue
                LOGGER.info("At time %s", time)

                for res in self.MOLRANGE:
                    resn = self.resid_to_lipid[res]
                    if resn not in lipidmolecules.STEROLS:
                        continue

                    neibs = self.neiblist[res][float(time)]
                    if len(neibs) != 0:
                        for neib in neibs:
                            orientation = self.calc_orientation(self.universe, res, neib)
                            resname_neib = self.resid_to_lipid[neib]
                            print("{: <12.2f}{: <10}{: <10}{: <15}{: <15}{: <15}".format(
                                time, res, resn, neib, resname_neib, orientation),
                                file=orifile)
    
    @staticmethod
    def scc_of_res(mda_uni, resid, tilt_correction=None):
        ''' Calculate the order parameter  '''
        resinfo = mda_uni.atoms.select_atoms("resid {}".format(resid))
        resname = list(set(resinfo.resnames))

        if tilt_correction is None:
            tilt_correction = [0, 0, 1] # Use z-axis

        if len(resname) > 1:
            raise ValueError("Atomselection resulted in selection of multiple residues")
        else:
            resname = resname[0]

        tailatms = lipidmolecules.scd_tail_atoms_of(resname)
        scds_of_tails = []
        for tail in tailatms:
            scds_of_atoms = []
            for atomindex in range(len(tail)-1):
                atm1, atm2 = tail[atomindex], tail[atomindex+1]
                coords12 = resinfo.atoms.select_atoms("name {} {}".format(atm1, atm2)).positions
                diffvector = np.subtract(*coords12)
                diffvector /= np.linalg.norm(diffvector)

                cos_angle = np.dot(diffvector, tilt_correction)
                scds_of_atoms += [0.5 * (3 * cos_angle**2 - 1)]
            scds_of_tails += [np.array(scds_of_atoms).mean()]

        LOGGER.debug("Res:  %s -- Scc of atoms: %s and of tails %s", resid, scds_of_atoms, scds_of_tails)
        return np.array(scds_of_tails).mean()

    def calc_orientation(self, mda_uni, resid, neib_resid):
        ''' Calculate the orientation angle of sterol molecule with its neighbors '''
                
        resn = self.resid_to_lipid[resid]

        C8_sterol = mda_uni.select_atoms("resid {} and name {}".format(resid, lipidmolecules.head_atoms_of(resn)[8])).positions
        C10_sterol = mda_uni.select_atoms("resid {} and name {}".format(resid, lipidmolecules.head_atoms_of(resn)[10])).positions
        C13_sterol = mda_uni.select_atoms("resid {} and name {}".format(resid, lipidmolecules.head_atoms_of(resn)[13])).positions
        #C19_sterol = mda_uni.select_atoms("resid {} and name {}".format(resid, lipidmolecules.head_atoms_of(resn)[19])).positions # just for checking
        
        v1 = np.subtract(C8_sterol,C10_sterol)
        v2 = np.subtract(C13_sterol,C10_sterol)
        v3 = np.cross(v1,v2) # vector peripendicular to the plane of sterol molecule
        #v4 = np.subtract(C19_sterol,C10_sterol) # just for checking
        
        neib_resn = self.resid_to_lipid[neib_resid]
        
        if neib_resn in lipidmolecules.STEROLS:

            C10_neib_sterol = mda_uni.select_atoms("resid {} and name {}".format(neib_resid, lipidmolecules.head_atoms_of(neib_resn)[6])).positions
            a = np.subtract(C10_neib_sterol,C10_sterol)
            
        elif neib_resn[:-2] in lipidmolecules.TAIL_ATOMS_OF.keys():

            P_lip = mda_uni.select_atoms("resid {} and name P".format(neib_resid)).positions
            a = np.subtract(P_lip,C10_sterol)

        orient_angle = np.arccos(np.dot(a[0], v3[0])/(np.linalg.norm(a)*np.linalg.norm(v3))) * (180/np.pi)
        #orient_angle1 = np.arccos(np.dot(v3[0], v4[0])/(np.linalg.norm(v3)*np.linalg.norm(v4))) * (180/np.pi) # just for checking
        orient_flag = 0
        
        if orient_angle < 90:
            orient_flag += 1
        
        LOGGER.debug("Res:  %s -- neib_resid: %s with orientation of %s", resid, neib_resid, orient_angle)
        return orient_flag

def calc_tilt(sysinfo, include_neighbors="global", filename="tilt.csv"):
    ''' End to end vector of each tail, average vector is substracted depending on include_neighbors variable
        include_neighbors can be "global", "local"(not implemented)
    '''
    u = sysinfo.universe
    len_traj = len(u.trajectory)

    # Lists as input to create pandas.DataFrame
    leaflet = []
    times = []
    angles = []
    vectors = []

    for t in range(len_traj):

        time = u.trajectory[t].time
        if time % sysinfo.dt != 0 or time > sysinfo.t_end or time < sysinfo.t_start:
            continue
        LOGGER.info("At time %s", time)

        leaflist = [[], []]
        for res in sysinfo.MOLRANGE:

            mda_res = u.residues[ res - 1 ]
            resn = sysinfo.resid_to_lipid[res]
            leaf = sysinfo.res_to_leaflet[res]

            if resn in lipidmolecules.STEROLS:
                continue

            t1_xyz = mda_res.atoms.select_atoms("name {}".format(lipidmolecules.tailcarbons_of(resn)[0][0])).positions
            t2_xyz = mda_res.atoms.select_atoms("name {}".format(lipidmolecules.tailcarbons_of(resn)[1][0])).positions

            tail1_xyz = mda_res.atoms.select_atoms("name {}".format(lipidmolecules.tailcarbons_of(resn)[0][-1])).positions - t1_xyz
            tail2_xyz = mda_res.atoms.select_atoms("name {}".format(lipidmolecules.tailcarbons_of(resn)[1][-1])).positions - t2_xyz
            tail1_xyz = tail1_xyz/np.linalg.norm(tail1_xyz)
            tail2_xyz = tail2_xyz/np.linalg.norm(tail2_xyz)

            leaflist[leaf] += [tail1_xyz, tail2_xyz]

        avg_vec_leaf1 = np.array(leaflist[0]).mean(axis=0)
        avg_vec_leaf2 = np.array(leaflist[1]).mean(axis=0)
        avg_vec_leaf1 = (avg_vec_leaf1/np.linalg.norm(avg_vec_leaf1))[0]
        avg_vec_leaf2 = (avg_vec_leaf2/np.linalg.norm(avg_vec_leaf2))[0]
        ang_l1 = angle_to_axis(avg_vec_leaf1)
        ang_l2 = angle_to_axis(avg_vec_leaf2)

        if ang_l1 > 90:
            ang_l1 = np.abs(ang_l1 - 180)
        if ang_l2 > 90:
            ang_l2 = np.abs(ang_l2 - 180)

        angles += [ang_l1, ang_l2]
        vectors +=  [avg_vec_leaf1, avg_vec_leaf2]
        times += [time, time]
        leaflet += [0, 1]

    vectors = np.array(vectors)
    dat = pd.DataFrame({"time":times, "leaflet":leaflet, "angle":angles, "x":vectors[:,0], "y":vectors[:,1], "z":vectors[:,2]})
    dat.to_csv(filename, index=False)


def angle_to_axis(vec: np.array, axis=np.array([0,0,1])) -> float:
    ''' Calculates angle of vector to axis (default: z axis)
        Returns angle in degree
    '''
    return np.arccos(np.dot(vec, axis)/np.linalg.norm(vec)) * (180/np.pi)
