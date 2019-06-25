'''
    This module focuses on the analysis of structural features of lipids in a bilayer
        - The order parameter S is calculate using class Order
        - The tilt of lipids is calculated in calc_avg_tilt

'''
import re
import numpy as np
import pandas as pd
import MDAnalysis as mda

from . import neighbors
from .neighbors import Neighbors
from .. import log
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

    def create_orderfile(self, mode="CC", outputfile='scd_distribution.dat'):
        '''
            Calculate order parameter S = 0.5 (3cos^2(a)-1) of lipids.
            a is the angle between the bilayer normal and (set with mode):
                CC -- Vector between Cn-Cn+2 carbon atoms of chain. Averaged over all the angles
        '''
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
                if self.t_start > time or self.t_end < time:
                    continue
                for res in self.MOLRANGE:
                    LOGGER.debug("At time %s and residue %s", time, res)
                    resname = self.resid_to_lipid[res]
                    LOGGER.debug("resname %s", resname)
                    if resname[:-2] not in lipidmolecules.TAIL_ATOMS_OF.keys()\
                        and resname not in lipidmolecules.STEROLS\
                        and resname not in lipidmolecules.PROTEINS:
                        LOGGER.debug("Skipping")
                        continue
                    scd_value = calc_s(self.universe, res)
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

    @staticmethod
    def scc_of_res(mda_uni, resid):
        ''' Calculate the order parameter  '''
        scds_of_atoms = []
        scds_of_tails = []
        resinfo = mda_uni.atoms.select_atoms("resid {}".format(resid))
        resname = list(set(resinfo.resnames))
        if len(resname) > 1:
            raise ValueError("Atomselection resulted in selection of multiple residues")
        else:
            resname = resname[0]
        # tailatms is like [Sn1atomlist, Sn2atomlist]
        tailatms = lipidmolecules.scd_tail_atoms_of(resname)
        for tail in tailatms:
            for atomindex in range(len(tail)-1):
                atm1, atm2 = tail[atomindex], tail[atomindex+1]
                coords12 = resinfo.atoms.select_atoms("name {} {}".format(atm1, atm2)).positions
                #diffvector = np.subtract(*[np.array(ar) for ar in coords12])
                diffvector = np.subtract(*coords12)
                normdiffvector = np.linalg.norm(diffvector)
                cos_angle = np.dot(diffvector, [0,0,1])/normdiffvector
                scds_of_atoms += [0.5 * (3 * cos_angle**2 - 1)]
            scds_of_tails += [np.array(scds_of_atoms).mean()]
        LOGGER.debug("Res:  %s -- Scc of atoms: %s and of tails %s", resid, scds_of_atoms, scds_of_tails)
        return np.array(scds_of_tails).mean()

def calc_tilt_princ_axis(sysinfo, neibdict=None, filename="tilt.dat"):
    ''' Calculate tilt angle (angle between principal axis of lipid and z-axis [0, 0, 1]) '''
    LOGGER.info("Calculating tilt with neibdict=%s to file %s", neibdict, filename)
    u = sysinfo.universe
    len_traj = len(u.trajectory)
    frames = []
    for t in range(len_traj):
        time = u.trajectory[t].time
        if time % sysinfo.dt != 0 or time > sysinfo.t_end or time < sysinfo.t_start:
            continue
        LOGGER.info("At time %s", time)
        resnames = []
        leaflet = []
        angles1 = []
        angles2 = []
        angles3 = []
        residues = []
        for res in sysinfo.MOLRANGE:
            mda_res = u.residues[ res - 1 ]
            residues.append(res)
            resnames.append(sysinfo.resid_to_lipid[res])
            leaflet.append(sysinfo.res_to_leaflet[res])
            ang1, ang2, ang3 = [np.arccos(np.dot(i, np.array([0,0,1]))) * (180/np.pi) for i in mda_res.atoms.principal_axes()]
            if ang1 >= 90:
                ang1 = np.abs(ang1 - 180)
            if ang2 >= 90:
                ang2 = np.abs(ang2 - 180)
            if ang3 >= 90:
                ang3 = np.abs(ang3 - 180)
            ang1 = [ang1]
            ang2 = [ang2]
            ang3 = [ang3]
            #if neibdict is not None:
            #    try:
            #        neibs = neibdict[res][time]
            #    except KeyError:
            #        continue
            #    for n in neibs:
            #        ang_n = (np.arccos(np.dot(res.atoms.principal_axes()[2], np.array([0,0,1]))) * (180/np.pi))
            #        if ang_n >= 90:
            #            ang_n = np.abs(ang_n - 180)
            #        ang.append(ang_n)
            ang1 = np.array(ang1).mean()
            ang2 = np.array(ang2).mean()
            ang3 = np.array(ang3).mean()
            angles1.append(ang1)
            angles2.append(ang2)
            angles3.append(ang3)
            #LOGGER.debug("%s %s", res, ang)
        dat = pd.DataFrame({"resid":residues, "angle1":angles1, "angle2":angles2, "angle3":angles3, "resname":resnames, "leaflet":leaflet})
        dat["time"] = time
        frames.append(dat)
    final = pd.concat(frames)
    final.to_csv(filename, index=False)

def calc_tilt_vectorwise(sysinfo, include_neighbors="global", filename="tilt.csv"):
    ''' End to end vector of each tail, average vector is substracted depending on include_neighbors variable
        include_neighbors can be "global", "local"(not implemented)
    '''
    #LOGGER.info("Calculating tilt with neibdict=%s to file %s", neibdict, filename)
    u = sysinfo.universe
    #zaxis = np.array([0.0, 0.0, 1.0])
    len_traj = len(u.trajectory)
    leaflet = []
    times = []
    angles = []
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
            if resn != "DPPC":
                continue
            t1_xyz = mda_res.atoms.select_atoms("name {}".format(lipidmolecules.tailcarbons_of(resn)[0][0])).positions
            t2_xyz = mda_res.atoms.select_atoms("name {}".format(lipidmolecules.tailcarbons_of(resn)[1][0])).positions

            tail1_atmstr = ' '.join([i for i in lipidmolecules.tailcarbons_of(resn)[0][-2:]])
            tail2_atmstr = ' '.join([i for i in lipidmolecules.tailcarbons_of(resn)[1][-2:]])

            tail1_xyz = (mda_res.atoms.select_atoms("name {}".format(tail1_atmstr)).positions - t1_xyz).mean(axis=0)
            tail1_xyz = tail1_xyz/np.linalg.norm(tail1_xyz)

            tail2_xyz = (mda_res.atoms.select_atoms("name {}".format(tail2_atmstr)).positions - t2_xyz).mean(axis=0)
            tail2_xyz = tail2_xyz/np.linalg.norm(tail2_xyz)

            leaflist[leaf] += [tail1_xyz, tail2_xyz]

        avg_vec_leaf1 = np.array(leaflist[0]).mean(axis=0)
        avg_vec_leaf2 = np.array(leaflist[1]).mean(axis=0)
        ang_l1 = angle_to_axis(avg_vec_leaf1)
        ang_l2 = angle_to_axis(avg_vec_leaf2)

        if ang_l1 > 90:
            ang_l1 = np.abs(ang_l1 - 180)
        if ang_l2 > 90:
            ang_l2 = np.abs(ang_l2 - 180)

        angles += [ang_l1, ang_l2]
        times += [time, time]
        leaflet += [0, 1]

    dat = pd.DataFrame({"time":times, "leaflet":leaflet, "angle":angles})
    dat.to_csv(filename, index=False)

#def calc_avg_tilt(inputfilename="inputfile", outputfname="bilayer_tilt.dat", average_over_neibs=True):
#    ''' Calculate the tilt of blocks of molecules (host + its neighbors) '''
#    systeminfo = Neighbors(inputfilename)
#    neiblist = neighbors.get_neighbor_dict()
#    uni = mda.Universe(systeminfo.gropath, systeminfo.trjpath)
#    len_traj = len(uni.trajectory)
#    with open(outputfname, "w") as outf:
#        print("{: <15}{: <5}{: <15}{: <15}{: <15}{: <9}{: <15}{: <15}{: <15}"\
#            .format('time', 'res', 'scd', 'tiltangle', 'cos_tiltangle', 'leaflet', 'x' , 'y', 'z'), file=outf)
#        for res in systeminfo.MOLRANGE:
#            LOGGER.info("At res %s", res)
#            for t in range(len_traj):
#                time = uni.trajectory[t].time
#                resname = systeminfo.resid_to_lipid[res]
#                xyz = uni.select_atoms("resid {} and name {}".format(res , lipidmolecules.central_atom_of(resname))).positions[0]
#                leaflet = systeminfo.res_to_leaflet[res]
#                try:
#                    neibs = neiblist[res][time]
#                except KeyError:
#                    continue
#
#                if average_over_neibs:
#                    lipid_block = [res] + neibs
#                else:
#                    lipid_block = [res]
#
#                LOGGER.debug("Resid %s at time %s and neibs %s", res, time, neibs)
#
#                cos_tiltangle = tilt_of_molecules(uni, lipid_block, time)
#                tiltangle = np.arccos(cos_tiltangle)*(180/np.pi)
#                if abs(tiltangle) > 90:
#                    tiltangle = tiltangle - 180
#                scd = 0.5 * (3 * cos_tiltangle**2 - 1)
#                print("{: <15}{: <5}{: <15.5f}{: <15.5f}{: <15.5f}{: <9}{: <15.5f}{: <15.5f}{: <15.5f}"\
#                    .format(time, res, scd, tiltangle, cos_tiltangle, leaflet, *xyz), file=outf)

def tilt_of_molecules(mda_universe, list_of_resids, time):
    ''' Calculates the tilt angle of a group of molecules in list_of_resids
        at a time <time / ps>
        Tilt is calculated as follows:
        angle between vector spanned by
            COM of HEAD
            COM of TAILHALF1
            COM of TAILHALF2
        and z axis
    '''
    def get_indices(wanted_atms, atoms_list):
        ''' Returns the mask to choose correct indices in mda.positions '''
        return np.isin(atoms_list, wanted_atms)

    # Set frame
    frame_n = int(time / mda_universe.trajectory.dt)
    mda_universe.trajectory[frame_n]

    # Get the tilt per resid
    tilts = []
    for resid in list_of_resids:
        #resid = resid - 1 # Because mda is weird
        selection = mda_universe.select_atoms("resid {}".format(resid))
        atms_res = selection.names
        resname = selection.resnames[0]

        LOGGER.debug("resid %s\natms_res %s\nresname %s", resid, atms_res, resname)
        #if lipidmolecules.is_sterol(resname):
        #    continue


        # Choose atoms along which tilt is to be calculated
        head_atmnames = lipidmolecules.head_atoms_of(resname)
        tail_atmnames = lipidmolecules.tailcarbons_of(resname)
        LOGGER.debug("Atomnames of tail: %s", tail_atmnames)
        tail_atmnames = [i for lt in tail_atmnames for i in lt]
        LOGGER.debug("Atomnames of head: %s", head_atmnames)
        tail1 = tail_atmnames[int(np.round(len(tail_atmnames)/2)):]
        tail2 = tail_atmnames[:int(np.round(len(tail_atmnames)/2))]
        LOGGER.debug("Names of tail1: %s\nand tail2: %s: ", tail1, tail2)

        # Get coords from mda_universe
        LOGGER.debug("Selection matrix %s", get_indices(head_atmnames, atms_res))
        head_com  = selection[get_indices(head_atmnames, atms_res)].center_of_mass(pbc=True)
        tail1_com = selection[get_indices(tail1, atms_res)].center_of_mass(pbc=True)
        tail2_com = selection[get_indices(tail2, atms_res)].center_of_mass(pbc=True)
        LOGGER.debug("COM of:\nhead:%s\ntail1:%s\ntail2:%s", head_com, tail1_com, tail2_com)

        # Calculate average tilt between vectors head -> tail1 and head -> tail2
        zvec = np.array([0, 0, 1])
        diff_ht1 = head_com-tail1_com
        diff_ht2 = head_com-tail2_com
        LOGGER.debug("Diff1: %s, diff2: %s", diff_ht1, diff_ht2)
        cos1 = np.dot(diff_ht1, zvec) / (np.linalg.norm(diff_ht1))
        cos2 = np.dot(diff_ht2, zvec) / (np.linalg.norm(diff_ht2))
        LOGGER.debug("cos1: %s, cos2 %s", cos1, cos2)
        avg_cos = (cos1 + cos2) / 2.0
        tilts.append(avg_cos)
    tiltangle = np.mean(tilts)
    LOGGER.debug("Average tilt of group %s", tiltangle)
    return tiltangle

def write_tilt_distr(systeminfo, filename="tilt.dat"):
    neiblist = neighbors.get_neighbor_dict()
    uni = mda.Universe(systeminfo.gropath, "{}/xtc/{}_{}.xtc".format(systeminfo.mdfilepath, systeminfo.system, systeminfo.temperature))
    len_traj = len(uni.trajectory)
    frames = []
    #with open(outputfname, "w") as outf:
    #    print("{: <15}{: <5}{: <15}{: <15}{: <15}{: <9}{: <15}{: <15}{: <15}"\
    #        .format('time', 'res', 'scd', 'tiltangle', 'cos_tiltangle', 'leaflet', 'x' , 'y', 'z'), file=outf)
    for t in range(len_traj):
        time = uni.trajectory[t].time
        if time % systeminfo.dt != 0 or time > systeminfo.t_end or time < systeminfo.t_start:
            continue
        LOGGER.info("At time %s", time)
        resnames = []
        leaflet = []
        angles = []
        residues = []
        for res in systeminfo.MOLRANGE:
            mda_res = uni.residues[ res - 1 ]
            residues.append(res)
            resnames.append(mda_res.resname)
            leaflet.append(systeminfo.res_to_leaflet[res])
            #try:
            #    neibs = neiblist[res][time]
            #except KeyError:
            #    continue

            tilt_angle = calc_tilt(mda_res)
            angles.append(tilt_angle)

            #if average_over_neibs:
            #    lipid_block = [res] + neibs
            #else:
            #    lipid_block = [res]

            LOGGER.debug("Resid %s at time %s ", res, time)

        dat = pd.DataFrame({"resid":residues, "tilt":angles, "resname":resnames, "leaflet":leaflet})
        dat["time"] = time
        frames.append(dat)
            #cos_tiltangle = tilt_of_molecules(uni, lipid_block, time)
            #tiltangle = np.arccos(cos_tiltangle)*(180/np.pi)
            #if abs(tiltangle) > 90:
            #    tiltangle = tiltangle - 180
            #scd = 0.5 * (3 * cos_tiltangle**2 - 1)
            #print("{: <15}{: <5}{: <15.5f}{: <15.5f}{: <15.5f}{: <9}{: <15.5f}{: <15.5f}{: <15.5f}"\
            #    .format(time, res, scd, tiltangle, cos_tiltangle, leaflet, *xyz), file=outf)
    final = pd.concat(frames)
    final.to_csv(filename, index=False)

def calc_tilt(mda_res, mode="end_to_end"):
    ''' Calculate tilt angle with respect to z axis '''
    resn = mda_res.resname
    head_xyz = mda_res.atoms.select_atoms("name {}".format(' '.join(lipidmolecules.head_atoms_of(resn)))).center_of_mass()
    if mode == "end_to_end":
        atmstr = ' '.join([i for t in np.array(lipidmolecules.tailcarbons_of(resn))[:,-3:] for i in t])
        tail_end = mda_res.atoms.select_atoms("name {}".format(atmstr)).center_of_mass()
        angle = np.arccos( np.dot( head_xyz - tail_end, [0, 0, 1] ) / np.linalg.norm(head_xyz - tail_end) ) * (180 / np.pi)
        if angle >= 90.0:
            angle = np.abs(angle - 180)
        return angle

def angle_to_axis(vec: np.array, axis=np.array([0,0,1])) -> float:
    ''' Calculates angle of vector to axis (default: z axis)
        Returns angle in degree
    '''
    return np.arccos(np.dot(vec, axis)/np.linalg.norm(vec)) * (180/np.pi)
