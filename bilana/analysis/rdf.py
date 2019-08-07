'''
    This module contains a class that defines own implementation to calculate the radial distribution function
    as well as a wrapper for the gmx rdf tool function radialdistribution
'''
import os
import re
import numpy as np
import bisect
from .. import log
from ..common import exec_gromacs
from ..systeminfo import SysInfo

LOGGER = log.LOGGER
GMXNAME = 'gmx'

def radialdistribution(systeminfo, ref, sel, seltype='atom', selrpos='atom', binsize=0.002, refprot=False):
    ''' Calculates 2D RDF of sel relative to ref only for one specific leaflet
        leaflet_assignment.dat (created in leaflets.py) is needed
        Using gmx rdf
    '''

    LOGGER.info("Calculating radial distribution function")
    LOGGER.info("Ref: %s\nSel: %s\n", ref, sel)
    os.makedirs(systeminfo.datapath+'/rdf', exist_ok=True)
    selectdict = {}

    if os.path.isfile("leaflet_assignment.dat"):
        res_in_leaf = []
        with open("leaflet_assignment.dat", "r") as leas:
            leas.readline()
            for line in leas:
                cols = line.split()
                if not int(cols[1]):
                    res_in_leaf.append(cols[0])
        select_one_leaflet = 'and resid '+' '.join(res_in_leaf)
    else:
        raise FileNotFoundError("Need leaflet_assignment.dat")
        #select_one_leaflet = 'and z<4'

    for selection in (('ref', ref), ('sel', sel)):

        if selection[0] == "ref" and refprot:
            selectstring = '{}'.format(selection[1])
        else:
            selectstring = '{} {}'.format(selection[1], select_one_leaflet)

        select_fname = '{}/select{}_{}'.format(systeminfo.temppath,
                                               selection[0], selection[1]).replace(" ", "_")
        selectdict.update({selection[1]:select_fname})

        with open(select_fname, "w") as selfile:
            print("Selection for {} is:\n{}\n".format(selection[0], selectstring))
            print(selectstring, file=selfile)

    outputfile = '{}/rdf/rdf_{}-{}.xvg'.format(systeminfo.datapath,
                                               ref, sel).replace(" ", "")
    outputfile_cn = '{}/rdf/nr_{}-{}.xvg'.format(systeminfo.datapath,
                                                 ref, sel).replace(" ", "")

    g_rdf_arglist = [
        GMXNAME, 'rdf', '-xy', '-xvg', 'none',
        '-f', systeminfo.trjpath, '-s', systeminfo.tprpath,
        '-o', outputfile,  '-cn', outputfile_cn,
        '-ref',  '-sf', selectdict[ref],
        '-sel',  '-sf', selectdict[sel],
        '-selrpos', selrpos, '-seltype', seltype, '-bin', str(binsize),
        ]
    out, err = exec_gromacs(g_rdf_arglist)
    rdf_log = 'rdf_{}-{}.log'.format(ref, sel).replace(" ", "")
    with open(rdf_log, "w") as logfile:
        print('STDERR:\n{}\n\nSTDOUT:\n{}'\
                .format(err, out), file=logfile)


############################# DEPRECATED ###########################################

#class calc_rdf_selfimplementation(SysInfo):
#    ''' This class holds an old implementation to calculate the radial distribution function
#        NOTE: DEPRECATED
#    '''
#    def __init__(self, refstring, selstring, neibdict,
#                 refstring2=None, selstring2=None,
#                 mode='SINGLE', nchol='off', outfname='rdf.dat',
#                 binwidth=0.002, overwriterdf=False,
#                 overwritegro=False, grofile=None, dim=2, inputfile="inputfile"):
#        '''
#            NOTE: This function is not really optimal
#            mode='COM' Calculate center of mass for selection
#            mode='ATM' Use plain selection
#            refstring --> the choice for references;
#                          string with similar syntax as in vmd
#            selstring --> choice for atoms to be considered around reference;
#                        #valid descriptions:
#                        #    -resname (only in combination with mode='COM')
#                        #    -atomname (can be list of names or single atom)
#         '''
#        super().__init__(inputfile)
#        self.binwidth = binwidth
#        self.overwrite = overwritegro
#        self.outname = outfname
#        self.grofile = grofile
#
#        self.dim = dim
#        self.mode = mode
#        self.ref = self.translate_selectionline(refstring)
#        self.sel = self.translate_selectionline(selstring)
#
#        if refstring2 is not None:
#            self.ref2 = self.translate_selectionline(refstring2)
#        else:
#            self.ref2 = self.ref
#        if selstring2 is not None:
#            self.sel2 = self.translate_selectionline(selstring2)
#        else:
#            self.sel2 = self.sel
#
#        print(self.ref, self.ref2)
#        if mode != 'COM' and self.ref[0] == 'resname'\
#                or self.sel[0] == 'resname':
#            raise ValueError('Combination of resname without mode=COM invalid')
#        self.nchol = nchol
#        if nchol != 'off':
#            self.neiblist = neibdict
#        if os.path.isfile(outfname) and not overwriterdf:
#            print("Data exists and I shall not overwrite.")
#        else:
#            self.run()
#
#    def run(self):
#        rdf = self.calc_rdf()
#        self.write_to_file(rdf, self.outname)
#
#    def write_to_file(self, rdf, outputfilename):
#        ''' Takes rdf dictionary and writes to file with outputfilename '''
#        print("Writing data to file")
#        with open(outputfilename, "w") as outf:
#            for r in sorted(rdf.keys()):
#                print('{: <10.3f}{: <10}'.format(r, rdf[r]), file=outf)
#
#    @staticmethod
#    def calc_distance(coord1, coord2, box_vectors, dim):
#        def create_pbc_copies(coord, boxvectors, dim):
#            pbc_copies = []
#            if dim == 2:
#                mirror_image_vectors = [
#                    (1, 0, 0),
#                    (0, 1, 0), (0, -1, 0),
#                    (-1, 0, 0),
#                    ]
#            elif dim == 3:
#                mirror_image_vectors = [
#                    (1, 0, 0),
#                    (0, 1, 0), (0, -1, 0),
#                    (-1, 0, 0),
#                    (0, 0, 1 ),
#                    (0, 0, -1),
#                    ]
#            for vec in mirror_image_vectors:
#                vec = np.array(vec[:dim])
#                pbc_translation = vec*boxvectors[:dim]
#                pbc_copy = coord[:dim] + pbc_translation
#                pbc_copies.append(pbc_copy)
#            return pbc_copies
#        distances = []
#        coords = [np.array(coord2[:dim])]
#        rmax = min(box_vectors[:dim])/2
#        #print("RMAX", rmax)
#        for xyz in create_pbc_copies(coord2, box_vectors, dim):
#            coords.append(xyz)
#        for coord_partner in coords:
#            distance = np.round(np.linalg.norm(coord_partner-coord1[:dim]), 5)
#            if distance <= rmax:
#                distances.append(distance)
#            #distances.append(distance)
#        #print("DISTANCES", distances)
#        if len(distances) > 1:
#            return min(distances)
#            #raise ValueError("More than one distance smaller rmax")
#        elif not distances:
#            return None
#        return distances[0]
#
#    def translate_selectionline(self, line):
#        ''' translates line from vmd syntax to
#            syntax usable for python and grofiles '''
#        print("Translating:", line)
#        choices = line.split(' and ')
#        conditions = []
#        for item in choices:
#            choicelist = item.split(' ')
#            description = choicelist[0]
#            choice_items = choicelist[1:]
#            conditions.append((description, choice_items))
#            print(description, choice_items)
#        return conditions
#
#    def gather_coords(self, grofile, time):
#        ''' !!! BEWARE: BAD CODE !!!
#            Reads grofile and creates histogram from choices
#            Center of mass is actually calculated as center of geometry!
#            Center of mass is always respective to resid
#        '''
#        def calculate_geometriccenter(coordinateinputlist, *args):
#            geocenter = np.array([0., 0., 0.])
#            for atomcoords in coordinateinputlist:
#                #print("ATOMCOORDS", atomcoords)
#                for ind, coord in enumerate(atomcoords):
#                    geocenter[ind]+=coord
#            geocenter= geocenter/len(coordinateinputlist)
#            #print(geocenter)
#            #sys.exit()
#            return np.array(geocenter)
#        name2colnr = {
#            'resid':1,
#            'resname':2,
#            'atom':3,
#            'index':4,
#            }
#        coordlist_ref = [[], []]
#        coordlist_sel = [[], []]
#        #========================RESID=======RESNAME======ATOMNAME=======INDEX===============X============Y============Z=======
#        regexp = re.compile(r'^([\s,\d]{5})([\w,\s]{5})([\d,\w,\s]{5})([\s*\d+]{5})\s*(-?\d+\.\d+\s*-?\d+\.\d+\s*-?\d+\.\d+).*')
#        regexp_boxline = re.compile(r'^\s+\d+\.\d+\s+\d+\.\d+\s+\d+\.\d+\s*$')
#        if self.mode == 'SINGLE':
#            for line in grofile:
#                match = regexp.match(line)
#                if match:
#                    #atomname = match.group(3).strip()
#                    resid = int(match.group(1).strip())
#                    xyz = np.array([float(i) for i in match.group(5).split()])
#                    #if xyz[2] < 4:
#                        #continue
#                    matched = True
#                    matched_cond = match.group(name2colnr[cond[0]]).strip()
#                    for cond in self.ref:
#                        if not matched_cond in cond[1]:
#                            matched = False
#                            break
#                    if matched:
#                        if self.nchol != 'off':
#                            neibs = self.neiblist[int(resid)][float(time)]
#                            neibs_chol = [self.resid_to_lipid[N]\
#                                          for N in neibs].count('CHL1')
#                            if neibs_chol == self.nchol:
#                                if self.res_to_leaflet[resid]:
#                                    coordlist_ref[0].append(xyz)
#                                else:
#                                    coordlist_ref[1].append(xyz)
#                        else:
#                            if self.dim == 2:
#                                if self.res_to_leaflet[resid]:
#                                    coordlist_ref[0].append(xyz)
#                                else:
#                                    coordlist_ref[1].append(xyz)
#                            else:
#                                coordlist_ref[0].append(xyz)
#                    matched = True
#                    for cond in self.sel:
#                        if not matched_cond in cond[1]:
#                            matched = False
#                            break
#                    if matched:
#                        if self.dim == 2:
#                            if self.res_to_leaflet[resid]:
#                                coordlist_sel[0].append(xyz)
#                            else:
#                                coordlist_sel[1].append(xyz)
#                        else:
#                            coordlist_sel[0].append(xyz)
#                elif regexp_boxline.match(line):
#                    boxvector = [float(i) for i in line.split()]
#                    break
#                else:
#                    continue
#        elif self.mode == 'COM':
#            COMlist_ref = []
#            COMlist_sel = []
#            resid_old = 1
#            #coordfilename = "z_coords_{}.xyz".format(time)
#            #with open(coordfilename, "w") as cfile:
#            for line in grofile:
#                match = regexp.match(line)
#                if match:
#                    #atomname = match.group(3).strip()
#                    resid = int(match.group(1).strip())
#                    if resid != resid_old:
#                        #print("NEW")
#                        #print("COMLIST_REF", COMlist_ref)
#                        if COMlist_ref:
#                            geo_center_ref = calculate_geometriccenter(
#                                COMlist_ref, self.dim)
#                            COMlist_ref = []
#                            if self.nchol != 'off':
#                                neibs = self.neiblist[int(resid_old)][float(time)]
#                                neibs_chol = [self.resid_to_lipid[N]\
#                                                for N in neibs].count('CHL1')
#                                if neibs_chol == self.nchol:
#                                    coordlist_ref[0].append(geo_center_ref)
#                            else:
#                                coordlist_ref[0].append(geo_center_ref)
#                        if COMlist_sel:
#                            geo_center_sel = calculate_geometriccenter(COMlist_sel, self.dim)
#                            #print("GEO_SEL_Z", geo_center_sel[-1])
#                            coordlist_sel[0].append(geo_center_sel)
#                            COMlist_sel = []
#                        #print("COORDLIST_REF", len(coordlist_ref[0]))
#                        #print("COORDDLIST_SEL", len(coordlist_sel))
#                    resid_old = resid
#                    xyz = np.array([float(i) for i in match.group(5).split()])
#                    if self.res_to_leaflet[resid]:
#                        continue
#                    matched = True
#                    for cond, cond2 in zip(self.ref, self.ref2):
#                        match_cond =  match.group(name2colnr[cond[0]]).strip()
#                        match_cond2 = match.group(name2colnr[cond2[0]]).strip()
#                        if not match_cond in cond[1]\
#                                and not match_cond2 in cond2[1]:
#                            matched = False
#                            break
#                    if matched:
#                        #print(resid)
#                        COMlist_ref.append(xyz)
#                    matched = True
#                    for cond in self.sel:
#                        if not match_cond in cond[1]:
#                            matched = False
#                            break
#                    if matched:
#                        COMlist_sel.append(xyz)
#                elif regexp_boxline.match(line):
#                    boxvector = [float(i) for i in line.split()]
#                    break
#                else:
#                    continue
#        return coordlist_sel, coordlist_ref, boxvector
#
#    def process_frame(self, grofile, time):
#        '''
#            Calculates the RDF for SEL around REF
#            for ONE TIMEFRAME in a GROFILE
#        '''
#        binwidth = self.binwidth
#        gofr = {}
#        coordlist_sel, coordlist_ref, boxvector = self.gather_coords(grofile,
#                                                                     time)
#        #print("COORDLISTS", len(coordlist_ref), len(coordlist_sel))
#        rmax = min(boxvector[:2])/2
#        histog = {}
#        histog_list = np.arange(0, rmax+2*binwidth, binwidth)
#        for i in histog_list: histog[i] = 0 # Initialize
#        for leaflet in range(1):
#            for coord1 in coordlist_ref[leaflet]:
#                for coord2 in coordlist_sel[leaflet]:
#                    if (coord1 == coord2).all():
#                        continue
#                    distance = self.calc_distance(coord1, coord2,
#                                                  boxvector, self.dim)
#                    if distance is None:
#                        continue
#                    interval = histog_list[bisect.bisect_right(histog_list,
#                                                               distance)
#                                          ]
#                    histog[interval] += 1
#        NumRefAtm = len(coordlist_ref[0])
#        NumSelAtm = len(coordlist_sel[0])
#        if self.dim == 2:
#            Box = (2*rmax)**2                   # Box area
#        elif self.dim == 3:
#            Box = (2*rmax)**3
#        if self.ref == self.sel:                # Box volume
#            Rho_zero = (NumSelAtm-1)/Box        # Density of ideal gas
#        else:
#            Rho_zero = NumSelAtm/Box            # Density of ideal gas
#        for r in histog_list:
#            if self.dim == 2:
#                bin_norm = np.pi*(((r+binwidth)**2)-(r**2))      # Area of bin
#            elif self.dim == 3:
#                bin_norm = (4/3)*(np.pi)*(((r+binwidth)**3)-(r**3)) # V of bin
#            a = histog[r]/(NumRefAtm*bin_norm*Rho_zero)
#            gofr[r] = a
#        return gofr
#    def calc_rdf(self):
#        ''' NOTE: Rewrite to use MDAnalysis '''
#        if self.grofile is None:
#            print("NOTE: This function is old and depends on old structure files.")
#            raise ValueError("Need trajectory file in gro format")
#        else:
#            grofile_output = self.grofile
#        gofr_all_frames = {}
#        with open(grofile_output,"r") as grofile:
#            for line in grofile:
#                if 't=' in line:
#                    time = float(line[line.index('t=')+2:].strip())
#                    print(time)
#                    if float(self.t_end) < time:
#                    #if time > 16000.:
#                        break
#                    else:
#                        gofr_frame = self.process_frame(grofile, time)
#                        for r in gofr_frame.keys():
#                            try:
#                                gofr_all_frames[r].append(gofr_frame[r])
#                            except KeyError:
#                                gofr_all_frames[r] = []
#                                gofr_all_frames[r].append(gofr_frame[r])
#        for r in sorted(gofr_all_frames.keys()):
#            gval = np.mean(gofr_all_frames[r])
#            gofr_all_frames[r] = gval
#            #print(r, gval)
#        return gofr_all_frames
#