
import numpy as np
import re
from ..definitions import lipidmolecules
from .. import log
LOGGER = log.LOGGER

def is_neighbor_in_leaflet(systeminfo_inst, neiblist):
    ''' Searches for interleaflet neighborhood '''
    leaflet_assign = systeminfo_inst.res_to_leaflet
    #nlipids = systeminfo_inst.NUMBEROFMOLECULES
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


def create_leaflet_assignment_file(sysinfo_obj, verbosity="INFO"):
    ''' Creates a file with that assigns all lipids to upper or lower leaflet
                        !Attention!
            !Flip flops of Cholesterol are not considered! Though should it?
    '''
    LOGGER.setLevel(verbosity)
    outputfilename = 'leaflet_assignment.dat'
    grofile_path = sysinfo_obj.gropath
    coord_head, coord_base = None, None
    with open(grofile_path, "r") as gfile, open(outputfilename, "w") as outf:
        print("{: <7} {: <5}".format('resid', 'leaflet'), file=outf)
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
                #logger.debug("Linepars: %s %s %s %s", resid, resname, atomname, coords)
                if resname not in sysinfo_obj.molecules:
                    continue
                last_tail_atm = lipidmolecules.scd_tail_atoms_of(resname)[0][-1]
                head_atm = lipidmolecules.central_atom_of(resname)
                #logger.debug("Atmn/last_tail %s/%s", atomname, last_tail_atm)
                if old_resid != resid:
                    LOGGER.debug("Resid/resname/atmname %s/%s/%s", resid, resname, atomname)
                    LOGGER.debug("Coords head/base: %s/%s", coord_head, coord_base)
                    new_coords = coord_head - coord_base
                    cos = np.dot(new_coords, np.array([0.0,0.0,1.0]))/np.linalg.norm(new_coords)
                    if cos <= 0:
                        sum_upper += 1
                        leaflet = 0
                    else:
                        sum_lower += 1
                        leaflet = 1
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
