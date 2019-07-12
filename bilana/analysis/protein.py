'''
    Protein related analysis
'''
import numpy as np
import MDAnalysis as mda
from .. import log
LOGGER = log.LOGGER


def calc_tilt(universe: mda.Universe, resid_up, resid_down, fname="TMD_tilt.dat"):
    ''' Calculate tilt related to angle between zaxis and resid_down --> resid_up
        Takes COM of resids
    '''
    fstr2 = '{: <15}{: <20}'
    fstr  = '{: <15}{: <20.5f}'
    with open(fname, "w") as outf:
        print(fstr2.format("time", "tilt"), file=outf)
        for t in range(universe.trajectory.n_frames):
            time = universe.trajectory[t].time
            LOGGER.info("At %s", time)
            zaxis = np.array([0, 0, 1])
            sel_u = universe.select_atoms("resid {}".format(resid_up))
            sel_d = universe.select_atoms("resid {}".format(resid_down))
            pos_u = sel_u.center_of_mass()
            pos_d = sel_d.center_of_mass()
            costilt = np.dot((pos_d - pos_u), zaxis)/np.linalg.norm(pos_d - pos_u)
            angle = np.arccos(costilt) * (180/np.pi)
            if angle > 90:
                angle -= 180
            print(fstr.format(time, abs(angle)), file=outf)
