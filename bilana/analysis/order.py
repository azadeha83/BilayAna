'''
    This module focuses on the analysis of structural features of lipids in a bilayer
        - The order parameter S is calculate using class Order
        - The tilt of lipids is calculated in calc_avg_tilt

'''
import re
import numpy as np
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
                if float(self.t_start) > time:
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

def calc_avg_tilt(inputfilename="inputfile", outputfname="bilayer_tilt.dat"):
    ''' Calculate the tilt of blocks of molecules (host + its neighbors) '''
    systeminfo = Neighbors(inputfilename)
    neiblist = systeminfo.get_neighbor_dict()
    uni = mda.Universe(systeminfo.gropath, systeminfo.trjpath)
    len_traj = len(uni.trajectory)
    with open(outputfname, "w") as outf:
        print("{: <15}{: <5}{: <15}{: <15}{: <15}{: <9}{: <15}{: <15}{: <15}"\
            .format('time', 'res', 'scd', 'tiltangle', 'cos_tiltangle', 'leaflet', 'x' , 'y', 'z'), file=outf)
        for res in systeminfo.MOLRANGE:
            for t in range(len_traj):
                time = uni.trajectory[t].time
                resname = systeminfo.resid_to_lipid[res]
                xyz = uni.select_atoms("resid {} and name {}".format(res -1, lipidmolecules.central_atom_of(resname))).positions[0]
                leaflet = systeminfo.res_to_leaflet[res]
                try:
                    neibs = neiblist[res][time]
                except KeyError:
                    continue

                lipid_block = [res] + neibs

                LOGGER.debug("Resid %s at time %s and neibs %s", res, time, neibs)

                cos_tiltangle = tilt_of_molecules(uni, lipid_block, time)
                tiltangle = np.arccos(cos_tiltangle)*(180/np.pi)
                if abs(tiltangle) > 90:
                    tiltangle = tiltangle - 180
                scd = 0.5 * (3 * cos_tiltangle**2 - 1)
                print("{: <15}{: <5}{: <15.5f}{: <15.5f}{: <15.5f}{: <9}{: <15.5f}{: <15.5f}{: <15.5f}"\
                    .format(time, res, scd, tiltangle, cos_tiltangle, leaflet, *xyz), file=outf)

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
        resid = resid - 1 # Because mda is weird
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
