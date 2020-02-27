import os
import re
import numpy as np
import pandas as pd
from ..definitions import lipidmolecules
from ..common import exec_gromacs
from .. import log
import MDAnalysis as mda

LOGGER = log.LOGGER

def is_neighbor_in_leaflet(systeminfo_inst, neiblist):
    ''' Searches for interleaflet neighborhood '''
    leaflet_assign = systeminfo_inst.res_to_leaflet
    host_has_interleafletneib = []
    for host in systeminfo_inst.MOLRANGE:
        neibs_times = neiblist[host]
        host_leaflet = leaflet_assign[host]
        for t in neibs_times:
            for neib in neibs_times[t]:
                neib_leaflet = leaflet_assign[neib]
                if neib_leaflet != host_leaflet:
                    host_has_interleafletneib.append([host, neib])
    LOGGER.info(host_has_interleafletneib)


def molecule_leaflet_orientation(atompos1: np.array, atompos2: np.array, axis=np.array([0.0, 0.0, 1.0])) -> int:
    ''' Takes two positions and returns 0 or 1 depending if molecule is oriented upside down or up '''
    print(atompos1,atompos2)
    new_coords = atompos1 - atompos2
    cos = np.dot(new_coords, axis) / np.linalg.norm(new_coords)
    return ( 0 if cos <= 0 else 1 )

def create_leaflet_assignment_file_1(sysinfo_obj, verbosity="INFO"):
    ''' Creates a file with that assigns all lipids to upper or lower leaflet
                        !Attention!
            !Flip flops of Cholesterol are not considered! Though should it?
    '''
    LOGGER.setLevel(verbosity)
    outputfilename = 'leaflet_assignment.dat'
    grofile_path = sysinfo_obj.gropath
    coord_head, coord_base = None, None

    with open(grofile_path, "r") as gfile, open(outputfilename, "w") as outf:
        print("{: <7} {: <5}".format('resid', 'leaflet'), file=outf) # HEADER

                #========================RESID=======RESNAME======ATOMNAME=======INDEX===============X============Y============Z=======
        regexp = re.compile(r'^([\s,\d]{5})([\w,\s]{5})([\d,\w,\s]{5})([\s*\d+]{5})\s*(-?\d+\.\d+\s*-?\d+\.\d+\s*-?\d+\.\d+).*')
        old_resid = 0
        sum_upper = 0
        sum_lower = 0

        for line in gfile:
            match = regexp.match(line)

            if match:
                resid = int(match.group(1).split()[0])

                if not old_resid:
                    old_resid = resid

                resname = match.group(2).split()[0].upper()
                atomname = match.group(3).split()[0]
                coords = [float(i) for i in match.group(5).split()]

                if resname not in sysinfo_obj.molecules:
                    continue

                last_tail_atm = lipidmolecules.scd_tail_atoms_of(resname)[0][-1]
                head_atm = lipidmolecules.central_atom_of(resname)

                if old_resid != resid:
                    LOGGER.debug("Resid/resname/atmname %s/%s/%s", resid, resname, atomname)
                    LOGGER.debug("Coords head/base: %s/%s", coord_head, coord_base)

                    leaflet = molecule_leaflet_orientation(coord_head, coord_base)

                    if leaflet:
                        sum_lower += 1
                    else:
                        sum_upper += 1

                    print("{: <7} {: <5}".format(old_resid, leaflet), file=outf)

                    old_resid = resid
                    coord_head = coord_base = None

                if atomname == head_atm:
                    coord_head = np.array(coords)
                    LOGGER.debug("Added head %s", atomname)
                if atomname == last_tail_atm:
                    coord_base = np.array(coords)
                    LOGGER.debug("Added tail %s", atomname)

        if coord_base is not None and coord_head is not None:
            new_coords = coord_head - coord_base
            cos = np.dot(new_coords, np.array([0.0,0.0,1.0]))/np.linalg.norm(new_coords)
            if cos <= 0:
                sum_upper += 1
                leaflet = 0
            else:
                sum_lower += 1
                leaflet = 1
            print("{: <7} {: <5}".format(old_resid, leaflet), file=outf)
        LOGGER.info("UP: %s LOW: %s", sum_upper, sum_lower)


def create_leaflet_assignment_file(sysinfo_obj, verbosity="INFO"):
    ''' Creates a file with that assigns all lipids to upper or lower leaflet
                        !Attention!
            !Flip flops of Cholesterol are not considered! Though should it?
    '''
    LOGGER.setLevel(verbosity)
    outputfilename = 'leaflet_assignment.dat'
    
    try:
        grofile_path = sysinfo_obj.initgropath
    except:
        grofile_path = sysinfo_obj.gropath

    w = mda.Universe(grofile_path)
    lipids = w.select_atoms('resname {}'.format(' '.join(sysinfo_obj.molecules)))
    com_lipids = lipids.center_of_geometry()
    leaflet_assign_flag = 0

    with  open(outputfilename, "w") as outf:
        print("{: <7} {: <5}".format('resid', 'leaflet'), file=outf)

        for resid in sysinfo_obj.MOLRANGE:
            resn = sysinfo_obj.resid_to_lipid[resid]
            head = w.select_atoms('resid {} and name {}'.format(resid, lipidmolecules.head_atoms_of(resn)[0]))
            com_head = head.center_of_geometry()

            if com_head[2] > com_lipids[2]:
                leaflet_assign_flag = 1
            else:
                leaflet_assign_flag = 0

            print("{: <7} {: <5}".format(resid, leaflet_assign_flag), file=outf)

def calc_density(systeminfo, selstr, outname="density.xvg", overwrite=False, **kw_den):
    '''
        Uses density calculation of gromacs
        1. Get index file using gmx select with
            gmx select -f ... -select <selstr>
        2. Run
            gmx density -f ... -d Z
            NOTE: Additional flags can be set adding with kw_den like b=3 converted to -b 3
    '''
    os.makedirs(systeminfo.datapath + "densities", exist_ok=True)
    TRJ = systeminfo.trjpath_whole
    TPR = systeminfo.tprpath
    NDX = systeminfo.temppath + "temp.ndx"
    OUT = systeminfo.datapath + "densities/" + outname
    if not overwrite and os.path.exists(OUT):
        LOGGER.warning("Density file already exists")
        return OUT

    # Get index file containing respective residue indices
    LOGGER.info("Creating index file...")
    commandstring = 'gmx select -f {} -s {} -on {} -select'.format(TRJ, TPR, NDX) ## dont forget the selstr
    cmd = commandstring.split() + [selstr]
    out, err = exec_gromacs(cmd)
    with open("gmx_select.log", "w") as f:
        print(err, file=f)
        print(out, file=f)

    # Run gmx density
    LOGGER.info("Run gmx density...")
    additional_input = []
    if kw_den:
        keys = kw_den.keys()
        keys = ["-"+i for i in keys]
        vals = kw_den.values()
        for z in zip(keys, vals):
            additional_input += list(z)
    commandstring = "gmx density -f {} -s {} -n {} -o {}  -d Z".format(TRJ, TPR, NDX, OUT)
    cmd = commandstring.split() + additional_input
    out, err = exec_gromacs(cmd)
    with open("gmx_density.log", "w") as f:
        print(err, file=f)
        print(out, file=f)

    return OUT


def calc_thickness(universe, ref_atomname, fname="bilayer_thickness"):
    ''' Calculate thickness of a bilayer structure using reference atoms ref_atomname '''
    fstr = "{: <15}{: <15}"
    fstr2 = "{: <15}{: <15.3f}"
    with open(fname, "w") as outf:
        print(fstr.format("time", "thickness"), file=outf)
        for t in range(universe.trajectory.n_frames):
            time= universe.trajectory[t].time
            LOGGER.info("At %s", time)
            sel = universe.select_atoms("name {}".format(ref_atomname))
            pos = sel.positions
            upper = pd.Series(pos[len(pos)//2:, 2])
            downer= pd.Series(pos[:len(pos)//2, 2])
            thickness =  downer.mean() - upper.mean()
            print(fstr2.format(time, thickness), file=outf)
