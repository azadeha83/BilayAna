'''
    Protein related analysis
'''
import os
import re
import bisect
import numpy as np
import pandas as pd
import MDAnalysis as mda
import MDAnalysis.analysis.helanal as hl
from . import leaflets
from .. import log
from .. import systeminfo
from ..definitions import structure_formats
from ..definitions import lipidmolecules

LOGGER = log.LOGGER

def run_helanal(systeminfo,
        selection="name CA",
        origin_pdbfile="origin.pdb",
        matrix_filename="bending_matrix.dat",
        summary_filename="summary.txt",
        screw_filename="screw.xvg",
        tilt_filename="local_tilt.xvg",
        bend_filename="local_bend.xvg",
        twist_filename="unit_twist.xvg",
        ref_axis=None,
        ):
    path_to_analysis_files = "{}tmdanalysis/".format(systeminfo.datapath)
    os.makedirs(path_to_analysis_files, exist_ok=True)
    prefix = path_to_analysis_files
    hl.helanal_trajectory(systeminfo.universe,
        origin_pdbfile=origin_pdbfile,
        matrix_filename=matrix_filename,
        summary_filename=summary_filename,
        screw_filename=screw_filename,
        tilt_filename=tilt_filename,
        bend_filename=bend_filename,
        twist_filename=twist_filename,
        ref_axis=ref_axis,
        prefix=prefix,
        )

def define_bilayer_head_region(sysinfo, overwrite=False):
    '''
        Divides bilayer z range in head and tail regions
        e.g.:
            <region> <zlow> <zhigh>
            head1 -20 0
            tail1 0 20
    '''
    from scipy.interpolate import UnivariateSpline
    head_regions = []
    
    sel = "P=resname {} and name P;P;".format(' '.join(sysinfo.PL_molecules))

    output = leaflets.calc_density(sysinfo, sel, outname="densities_P.xvg", xvg="none", overwrite=True)

    dat = pd.read_table(output, header=None, delim_whitespace=True, names=["z", "dens"])
    LOGGER.debug("\n%s", dat)
    dat["zerogrps"] = (dat.dens.diff(1) != 0).astype(int).cumsum()
    dat = dat.groupby("zerogrps").mean() 
    dat["mask"] = dat.z < dat[dat.dens == 0].iloc[1].z
    for i, fr in  dat.groupby("mask"):

        ### Just use min max values of distributions as head region borders
        LOGGER.debug("Frame %s", i)
        LOGGER.debug("\n%s", fr)
        fr = fr[fr.dens > 0]
        r1 = fr.z.min()
        r2 = fr.z.max()
        head_regions.append( (r1, r2) )

        #### Use FWHM as head region borders
        #spline = UnivariateSpline(fr.z, fr.dens-np.max(fr.dens)/2, s=0)
        #r1, r2 = spline.roots() # find the roots
        #head_regions.append( (r1, r2) )

    LOGGER.debug("head_regions: %s", head_regions)
    return head_regions

def define_bilayer_center(sysinfo, overwrite=False):
    '''
        Defines the middlepoint of the bilayer
    '''
    atmnames = []
    for PL in sysinfo.PL_molecules:
        tail = lipidmolecules.tail_atoms_of(PL)
        tail = [i[-1] for i in tail]
        atmnames += tail
    atmnames = list( set( atmnames ) )
    sel = "tailends=resname {} and name {};tailends;".format( ' '.join(sysinfo.PL_molecules), ' '.join(atmnames) )
    output = leaflets.calc_density(sysinfo, sel, outname="densities_P.xvg", xvg="none", overwrite=overwrite)
    dat = pd.read_table(output, header=None, delim_whitespace=True, names=["z", "dens"])
    return dat.mean().z

def get_bilayer_regions(sysinfo, overwrite=False):
    ''' 
       Create file with
       region llower lupper
       head1   0.1    0.4
        .       .      .
        .       .      .
    '''
    headregions    = define_bilayer_head_region(sysinfo, overwrite=overwrite)
    bilayer_center = define_bilayer_center(sysinfo, overwrite=overwrite)
    return sorted([ *headregions[0], bilayer_center, *headregions[1], ])
    

def create_leaflet_assignment_prot(sysinfo, outputfilename="leaflet_assignment_prot.csv", overwrite=False):
    '''
        Append leaflet assignment for each protein residue to 
        leaflet_assignment.dat created in leaflet.create_leaflet_assignment_file
        < resid > < leaflet > < region >
    '''

    resids       = []
    leaflet_list = []
    regions      = []

    regions_name = {0:"solvent1", 1:"head1", 2:"tail1", 3:"tail2", 4:"head2", 5:"solvent2"}
    # list of region borders:  [ h1/t1, t1/t2, t2/h2, ]
    regions_range = get_bilayer_regions(sysinfo)
    for res in sysinfo.MOLRANGE_PROT:
        resname = sysinfo.resid_to_lipid[res]
        sel = "res=resname {} and name CA and resid {};res;".format( resname, res )
        output = leaflets.calc_density(sysinfo, sel, outname="densities_res{}.xvg".format(res), xvg="none", overwrite=overwrite)
        dat = pd.read_table(output, header=None, delim_whitespace=True, names=["z", "dens"])
        dat = dat[dat.dens > 0]
        dat["wt"] = dat.dens / dat.dens.sum()
        zmean = (dat.z  * dat.wt).sum()
        region_ind = bisect.bisect(regions_range, zmean)
        regname = regions_name[region_ind]
        if region_ind < 2: # In leaflet 1 
            leaflet = 0
        else:
            leaflet = 1
        resids.append(res)
        leaflet_list.append(leaflet)
        regions.append(regname)
        LOGGER.debug("Resid %s in %s with zmean %s and range %s", res, regname, zmean, regions_range)

    dat = pd.DataFrame({"resid":resids, "leaflet":leaflet_list, "region":regions})
    dat.to_csv(outputfilename, index=False)
    return dat 

def get_reference_resids(inputfilename="leaflet_assignment_prot.csv"):
    '''
        Returns list of two lists containing resids for both head regions
        List indices correspond to leaflet value
    '''
    try:
        dat = pd.read_csv(inputfilename)
    except FileNotFoundError:
        LOGGER.info("File %s not found", inputfilename)
        raise
        #create_leaflet_assignment_prot(systeminfo.SysInfo(), outputfilename=inputfilename)
    head1_res = list(dat[dat.region == "head1"].resid)
    head2_res = list(dat[dat.region == "head2"].resid)
    return [head1_res, head2_res] # Like dict with index == leaflet 

#def get_reference_positions(sysinfo, inputfilename="leaflet_assignment_prot.csv"):
#    '''
#        Returns a list with protein positions
#            [pos_leaf1, pos_leaf2]
#    '''
#    reference_resids = get_reference_resids(inputfilename=inputfilename)
#    pos_leaf1 = sysinfo.universe.sel
    
#def get_level_dict(filename="origin.pdb", path="datafiles/tmdanalysis"):
#    '''
#        Reads output from helanal <origin.pdb> and 
#        creates a dictionary with origin coordinates of all amino acids per frame
#    '''
#    outdict = {0:{}}
#    frame_counter = 0
#    fname = "{}/{}".format(path, filename)
#    regex = structure_formats.REGEXP_PDB
#
#    with open(path+filename, "r") as f:
#        for line in f:
#            match = regex.match(line)
#
#            if match is not None:
#                ld = match.groupdict()
#                outdict[frame_counter][ ld["resid"] ] = np.array([ld["X"], ld["Y"], ld["Z"]])
#                
#            if "ENDMDL" in line:
#                frame_counter += 1
#                outdict[frame_counter] = {}


def calc_tilt_end_to_end(universe: mda.Universe, resid_up, resid_down, fname="TMD_tilt.dat"):
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

