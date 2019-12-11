'''
    This module focuses on the analysis of structural features of lipids in a bilayer
        - The order parameter S is calculate using class Order
        - The tilt of lipids is calculated in calc_avg_tilt

'''
import os
import numpy as np
import pandas as pd
import multiprocessing as mp
import MDAnalysis as mda

from . import neighbors
from .neighbors import Neighbors
from .. import log
from ..systeminfo import SysInfo
from ..common import loop_to_pool
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
        self.neiblist = neighbors.get_neighbor_dict()
        # Change dict entries from neiblist[host][time] to neiblist[time][host]

    def create_orderfile(self, mode="CC", outputfile='scd_distribution.dat', with_tilt_correction="tilt.csv", parallel=True):
        '''
            Calculate order parameter S = 0.5 (3cos^2(a)-1) of lipids.
            a is the angle between the bilayer normal and (set with mode):
                CC -- Vector between Cn-Cn+2 carbon atoms of chain. Averaged over all the angles

            if with tilt_correction avg tilt angle per time is read from name given in variable
        '''
        if not self.trjpath_whole:
            LOGGER.error("Is input trajectory in mda.universe made whole? If not, this may result in spurious results")
            raise FileNotFoundError("No whole trajectory found")

        def catch_callback(result):
            ''' Catches callback value of pool.apply_async workers and puts into outputlist'''
            outputlist.append(result)
        neiblist_t = pd.DataFrame(self.neiblist).transpose().to_dict()
        outputlist = [] # Don't change this name as function catch_callback uses it!

        if parallel:
            pool = mp.Pool(len(os.sched_getaffinity(0)), maxtasksperchild=1)

        # If tilt correction is activated read file with tilt information or create it
        if with_tilt_correction:
            try:
                dat = pd.read_csv(with_tilt_correction)
            except FileNotFoundError:
                LOGGER.info("Could not find tilt file. Creating new one.")
                calc_tilt(self)
                dat = pd.read_csv(with_tilt_correction)

        # Define how the order parameter has to be calculated
        if mode == "CC":
            calc_s = self.scc_of_res
        elif mode == "profile":
            calc_s = self.get_order_profile
        else:
            raise ValueError("Mode not yet implement use mode=CC (default)")

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
    def scc_of_res(res, tailatms, all_atmnames, positions, tilt_correction=None):
        ''' Calculate the order parameter  '''

        if tilt_correction is None:
            tilt_correction = [0, 0, 1] # Use z-axis


        scds_of_tails = []
        for tail in tailatms:
            scds_of_atoms = []
            for atomindex in range(len(tail)-1):

                atm1, atm2 = tail[atomindex], tail[atomindex+1]

                mask1 =  ( all_atmnames.names == tail[ atomindex ], all_atmnames.resids == res )
                mask2 =  ( all_atmnames.names == tail[ atomindex + 1 ], all_atmnames.resids == res )
                atm1, atm2 = positions[ mask1[0] & mask1[1] ][0], positions[ mask2[0] & mask2[1] ][0]

                diffvector = atm2 - atm1
                diffvector /= np.linalg.norm(diffvector)
                LOGGER.debug("Diffvector %s", diffvector)

                cos_angle = np.dot(diffvector, tilt_correction)
                LOGGER.debug("Resulting cos %s", cos_angle)
                scds_of_atoms.append( 0.5 * ( ( 3 * (cos_angle**2)) - 1 )  )

            scds_of_tails.append( np.array(scds_of_atoms).mean() )
        LOGGER.debug("Scds of res %s: %s", res, scds_of_tails)

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
    DEBUG = False
    if DEBUG:
        LOGGER.setLevel("DEBUG")
    u = sysinfo.universe
    len_traj = len(u.trajectory)
    inpargs = []
    #angles  = []
    #vectors = []
    #times   = []
    #leaflet = []


    LOGGER.info("gather input arguments")
    for t in range(len_traj):

        time = u.trajectory[t].time
        if time % sysinfo.dt != 0 or time > sysinfo.t_end or time < sysinfo.t_start:
            continue
        #LOGGER.info("At time %s", time)

        inpargs.append( (time, sysinfo.MOLRANGE, sysinfo.resid_to_lipid, sysinfo.res_to_leaflet, u.atoms, u.atoms.positions) )

        #Leaflist = [[], []]
        #for res in sysinfo.MOLRANGE:

        #    mda_res = u.atoms.select_atoms("resid {}".format(res))
        #    resn = sysinfo.resid_to_lipid[res]
        #    leaf = sysinfo.res_to_leaflet[res]

        #    if resn in lipidmolecules.STEROLS:
        #        continue

        #    # first carbon of tail
        #    t1_xyz = mda_res.atoms.select_atoms("name {}".format(lipidmolecules.tailcarbons_of(resn)[0][0])).positions
        #    t2_xyz = mda_res.atoms.select_atoms("name {}".format(lipidmolecules.tailcarbons_of(resn)[1][0])).positions

        #    # get vector c_last - c_first
        #    tail1_xyz = mda_res.atoms.select_atoms("name {}".format(lipidmolecules.tailcarbons_of(resn)[0][-1])).positions - t1_xyz
        #    tail2_xyz = mda_res.atoms.select_atoms("name {}".format(lipidmolecules.tailcarbons_of(resn)[1][-1])).positions - t2_xyz

        #    # normalize vectors
        #    tail1_xyz = tail1_xyz / np.linalg.norm(tail1_xyz)
        #    tail2_xyz = tail2_xyz / np.linalg.norm(tail2_xyz)

        #    # leaf is either 0 or 1: append avg vector to correct leaflet list
        #    leaflist[leaf] += [tail1_xyz, tail2_xyz]

        ## get average tilt vector per leaflet
        #avg_vec_leaf1 = np.array(leaflist[0]).mean(axis=0)
        #avg_vec_leaf2 = np.array(leaflist[1]).mean(axis=0)
        #avg_vec_leaf1 = (avg_vec_leaf1/np.linalg.norm(avg_vec_leaf1))[0]
        #avg_vec_leaf2 = (avg_vec_leaf2/np.linalg.norm(avg_vec_leaf2))[0]

        ## calculate the angle of average vector to z axis
        #ang_l1 = angle_to_axis(avg_vec_leaf1)
        #ang_l2 = angle_to_axis(avg_vec_leaf2)
        #if ang_l1 > 90:
        #    ang_l1 = np.abs(ang_l1 - 180)
        #if ang_l2 > 90:
        #    ang_l2 = np.abs(ang_l2 - 180)

        ## fill lists this will become line in the final dataframe
        #angles  += [ang_l1, ang_l2]
        #vectors +=  [avg_vec_leaf1, avg_vec_leaf2]
        #times   += [time, time]
        #leaflet += [0, 1]

    LOGGER.info("running tasks")
    if parallel and not DEBUG:
        output = loop_to_pool(_calc_tilt_output, inpargs)
    else:
        output = []
        for inp in inpargs:
            print("at time:", inp[0])
            output.append( _calc_tilt_output(*inp) )
    LOGGER.info("sort and write dataframe")
    dat = pd.concat(output)
    dat = dat.sort_values(by=["time"])

    #vectors = np.array(vectors)
    #dat = pd.DataFrame({"time":times, "leaflet":leaflet, "angle":angles, "x":vectors[:,0], "y":vectors[:,1], "z":vectors[:,2]})
    dat.to_csv(filename, index=False)

def _calc_tilt_output(time, molrange, resid_to_lipid, res_to_leaflet, all_atms, all_positions,):
        leaflist = [[], []]

        for res in molrange:

            #mda_res = u.atoms.select_atoms("resid {}".format(res))
            #mda_res = all_atms[all_atms.resids == res]
            resn = resid_to_lipid[res]
            leaf = res_to_leaflet[res]

            if resn in lipidmolecules.STEROLS:
                continue

            # first carbon of tails
            #t1_xyz = mda_res.atoms.select_atoms("name {}".format(lipidmolecules.tailcarbons_of(resn)[0][0])).positions
            #t2_xyz = mda_res.atoms.select_atoms("name {}".format(lipidmolecules.tailcarbons_of(resn)[1][0])).positions
            masks1 = ( all_atms.names == lipidmolecules.tailcarbons_of(resn)[0][0], all_atms.resids == res )
            masks2 = ( all_atms.names == lipidmolecules.tailcarbons_of(resn)[1][0], all_atms.resids == res )
            t1_xyz = all_positions[ masks1[0] & masks1[1] ][0]
            t2_xyz = all_positions[ masks2[0] & masks2[1] ][0]

            # get vector c_last - c_first
            #tail1_xyz = mda_res.atoms.select_atoms("name {}".format(lipidmolecules.tailcarbons_of(resn)[0][-1])).positions - t1_xyz
            #tail2_xyz = mda_res.atoms.select_atoms("name {}".format(lipidmolecules.tailcarbons_of(resn)[1][-1])).positions - t2_xyz
            masks1 = ( all_atms.names == lipidmolecules.tailcarbons_of(resn)[0][-1], all_atms.resids == res )
            masks2 = ( all_atms.names == lipidmolecules.tailcarbons_of(resn)[1][-1], all_atms.resids == res )
            tail1_xyz = all_positions[ masks1[0] & masks1[1] ][0] - t1_xyz
            tail2_xyz = all_positions[ masks2[0] & masks2[1] ][0] - t2_xyz

            # normalize vectors
            tail1_xyz = tail1_xyz / np.linalg.norm(tail1_xyz)
            tail2_xyz = tail2_xyz / np.linalg.norm(tail2_xyz)

            LOGGER.debug("tail1_xyz %s", tail1_xyz)

            # leaf is either 0 or 1: append avg vector to correct leaflet list
            leaflist[leaf] += [tail1_xyz, tail2_xyz]
        LOGGER.debug("leaflist: %s", leaflist)

        # get average tilt vector per leaflet
        avg_vec_leaf1 = np.array(leaflist[0]).mean(axis=0)
        avg_vec_leaf2 = np.array(leaflist[1]).mean(axis=0)
        avg_vec_leaf1 = (avg_vec_leaf1/np.linalg.norm(avg_vec_leaf1))
        avg_vec_leaf2 = (avg_vec_leaf2/np.linalg.norm(avg_vec_leaf2))

        LOGGER.debug("avg_vec_leaf1 %s", avg_vec_leaf1)

        # calculate the angle of average vector to z axis
        ang_l1 = angle_to_axis(avg_vec_leaf1)
        LOGGER.debug("ang1 %s", ang_l1)
        ang_l2 = angle_to_axis(avg_vec_leaf2)
        if ang_l1 > 90:
            ang_l1 = np.abs(ang_l1 - 180)
        if ang_l2 > 90:
            ang_l2 = np.abs(ang_l2 - 180)

        # fill lists this will become line in the final dataframe
        angles  = [ang_l1, ang_l2]
        vectors =  [avg_vec_leaf1, avg_vec_leaf2]
        times   = [time, time]
        leaflet = [0, 1]
        vectors = np.array(vectors)
        dat = pd.DataFrame({"time":times, "leaflet":leaflet, "angle":angles, "x":vectors[:,0], "y":vectors[:,1], "z":vectors[:,2]})
        return dat



def angle_to_axis(vec: np.array, axis=np.array([0,0,1])) -> float:
    ''' Calculates angle of vector to axis (default: z axis)
        Returns angle in degree
    '''
    return np.arccos(np.dot(vec, axis)/np.linalg.norm(vec)) * (180/np.pi)
