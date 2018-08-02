''' Module for automating analysis using Gromacs '''

import subprocess
import re
import os
#import sys 
from time import localtime, strftime
import time as Time
import numpy as np
from bilana import lipidmolecules #,systeminfo
import bisect

#from src.systeminfo import mysystem
#global mysystem

#inputfile = sys.argv[1]
#mysystem = systeminfo.SysInfo(inputfile)


gmx_exec = 'gmx' #'gmx_5.0.4_mpi'
os.environ["GMX_MAXBACKUP"] = "-1"

def exec_gromacs(cmd,inp_str=None): 
    ''' Execute Gromacs commands.
        arglist (cmd) is list of arguments like "['gmx cmd', '-f', 'tprfile', '-e', 'en.edr']"
    '''
    if inp_str is None:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = proc.communicate()
    else:
        try:
            inp_str = inp_str.encode()
        except AttributeError:
            pass
        proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = proc.communicate(inp_str)
    proc.wait()
    proc.stdout.close()
    proc.stderr.close()
    if proc.returncode == 1:
        #print("Failed to execute command:", ' '.join(cmd))
        try:
            err = err.decode()
            print(err)
        except UnboundLocalError:
            pass
        raise ChildProcessError('Failed to execute command "{}" and input "{}"'.format(' '.join(cmd), inp_str))
    return out, err

def produce_gro(mysystem, grofilename='/traj_complete.gro'):
    ''' Converts trajectory .trr or .xtc (binary) file to .gro (ASCII) file '''
    print(strftime("%H:%M:%S :", localtime()), 'Converting trajectory-file to structure-file ...\n')
    generate_index_file(mysystem.gropath, mysystem.molecules)
    os.makedirs(mysystem.datapath+'/grofiles', exist_ok=True)
    grofile_output = ''.join([mysystem.datapath, '/grofiles/', grofilename])
    print(strftime("%H:%M:%S :", localtime()),"... Start conversion from .trj to .gro ...")
    print(mysystem.molecules)
    sorted_mols = []
    ### REFACTOR HERE!!! 
    if 'DPPC' in mysystem.molecules:
        sorted_mols.append('DPPC')
    if 'DUPC' in mysystem.molecules:
        sorted_mols.append('DUPC')
    if 'CHL1' in mysystem.molecules:
        sorted_mols.append('CHL1')
    if 'CHIM' in mysystem.molecules:
        sorted_mols.append('CHIM')
    if 'ERG' in mysystem.molecules:
        sorted_mols.append('ERG')
    inp_str = str('_'.join(sorted_mols)+'\n').encode()
    gmx_traj_arglist = [
        gmx_exec, 'trjconv', '-s', mysystem.tprpath, '-f', mysystem.trjpath,
        '-o', grofile_output, '-n', 'index.ndx',
        '-b', str(mysystem.t_start), #'-e', str(self.sysinfo.t_end),
        '-dt', str(mysystem.dt),
        '-pbc', 'whole',
        ]
    out, err = exec_gromacs(gmx_traj_arglist, inp_str)
    with open("gmx_traj_compl.log","w") as logfile:
        logfile.write(err.decode())
        logfile.write(150*'_')
        logfile.write(out.decode())
        logfile.write(150*'_')
    return grofile_output

def trajectory_to_gro(systeminfo, overwrite='off', atomlist=None, lipids='all'):
    os.makedirs(systeminfo.datapath+'/grofiles', exist_ok=True)
    getcoords = {}
    generate_index_file(systeminfo.gropath, systeminfo.molecules)
    print('Converting trajectory-file to structure-file ...\n')
    if lipids != 'all':
        molecules = [lipids]
    else:
        if "CHOL" in systeminfo.molecules:
            molecules_new = systeminfo.molecules.copy()
            molecules_new.remove("CHOL")
            #print(molecules_new)
        else:
            molecules_new = systeminfo.molecules
        molecules = molecules_new
    for lipidmolecule in molecules:
        grofile_output = ''.join([systeminfo.datapath, '/grofiles', '/calc_scd_for', str(lipidmolecule), '.gro'])
        if atomlist == None:
            atomlist = lipidmolecules.scd_tail_atoms_of[lipidmolecule]+['P', 'O3']
        print(strftime("%H:%M:%S :", localtime()),"Processing {} ...".format(lipidmolecule))
        inp_str = str(lipidmolecule).encode()
        if os.path.isfile(grofile_output) and overwrite=='off':
            pass
        else:
            print(strftime("%H:%M:%S :", localtime()),"... Start conversion from .trj to .gro ...")
            gmx_traj_arglist = [\
                gmx_exec,'trjconv','-s',systeminfo.tprpath, '-f',systeminfo.trjpath,\
                '-o',grofile_output,\
                '-b', str(systeminfo.t_start), '-e', str(systeminfo.t_end),\
                '-dt', str(systeminfo.dt),\
                '-pbc', 'whole',\
                ]
            out,err = exec_gromacs(gmx_traj_arglist,inp_str)      #Creates a gro file containing {timeframe:resid:atoms:xyz}
            with open("gmx_traj.log","w") as logfile:
                logfile.write(err.decode())
                logfile.write(150*'_')
                logfile.write(out.decode())
                logfile.write(150*'_')
        with open(grofile_output,"r") as grofile:
            print(strftime("%H:%M:%S :", localtime()),"... read data from .gro-file ...")
            regexp = re.compile(r'[\s]*\d+'+lipidmolecule)
            for line in grofile:
                if 't=' in line:
                    re_time = re.compile(r'.*t=\s+(\d+\.\d+).*')
                    time = float(re_time.match(line).group(1))
                    #print("...at time {}".format(time), end="\r")
                    if float(systeminfo.t_end) < time:
                        #print("breaking at", time)
                        break
                #print("Match", regexp.match(line), line[:15], end='\n')
                #print("Time match is", float(systeminfo.t_start)<=time, end='\n')
                if float(systeminfo.t_start) <= time and regexp.match(line) != None:
                    #print("Reading data at", time, end='\n')
                    atom = line[9:15].strip()
                    lipidtype = line[5:9]
                    if atom not in atomlist and lipidtype not in molecules: continue
                    resid = line[:5].strip()
                    coordinates = [float(x) for x in line[20:44].split()]
                    keytuple = (str(lipidtype),time,int(resid),str(atom))
                    getcoords.update({keytuple:coordinates})
        print(strftime("%H:%M:%S :", localtime()),"Finished reading.")
    return getcoords, time-1000.0

def generate_index_file(gropath, molecules):
    ''' Creates an indexfile with all relevant entries + one that contains all lipids '''
    molstrings = ['r {}'.format(i) for i in molecules]
    inp_str = ' | '.join(molstrings)
    inp_str += '\n q \n'
    #print(inp_str)
    inp_str = inp_str.encode()
    gmx_arglist = [
        gmx_exec, 'make_ndx', '-f', gropath,
        ]
    out, err = exec_gromacs(gmx_arglist, inp_str)
    with open("gmx_index.log","w") as logfile:
        logfile.write('__________ STDERR ___________')
        logfile.write(err.decode())
        logfile.write('__________ STDOUT ___________')
        logfile.write(out.decode())

class calc_rdf_selfimplementation():
    def __init__(self, _systeminfo, refstring, selstring, refstring2=None, selstring2=None, 
                 mode='SINGLE', nchol='off', outfname='rdf.dat', binwidth=0.002, overwriterdf=False,
                 overwritegro=False, grofile=None, dim=2):
        '''
            mode='COM' Calculate center of mass for selection
            mode='ATM' Use plain selection
            refstring --> the choice for references; string with same syntax as in vmd (NOT COMPLETE AT ALL!)
            selstring --> choice for atoms to be considered around reference; string with same syntax as in vmd (NOT COMPLETE AT ALL!)

                        #valid descriptions:
                        #    -resname (only in combination with mode='COM')
                        #    -atomname (can be list of names or single atom)
         '''
        self.sysinfo = _systeminfo
        self.binwidth = binwidth
        self.overwrite = overwritegro
        self.outname = outfname
        self.grofile = grofile

        self.dim = dim
        self.mode = mode
        self.ref = self.translate_selectionline(refstring)
        self.sel = self.translate_selectionline(selstring)
        
        if refstring2 is not None:
            self.ref2 = self.translate_selectionline(refstring2)
        else:
            self.ref2 = self.ref
        if selstring2 is not None:
            self.sel2 = self.translate_selectionline(selstring2)
        else:
            self.sel2 = self.sel

        print(self.ref, self.ref2)
        if mode != 'COM' and self.ref[0] == 'resname' or self.sel[0] == 'resname':
            raise ValueError('Combination of resname without mode=COM invalid')
        self.nchol = nchol
        if nchol != 'off':
            self.neiblist = Neighbors(_systeminfo).get_neighbor_dict()
        if os.path.isfile(outfname) and not overwriterdf:
            print("Data exists and I shall not overwrite.")
        else:
            self.run()

    def run(self):
        rdf = self.calc_rdf()
        self.write_to_file(rdf, self.outname)
    def write_to_file(self, rdf, outputfilename):
        ''' Takes rdf dictionary and writes to file with outputfilename '''
        print("Writing data to file")
        with open(outputfilename, "w") as outf:
            # for i in np.arange(min(rdf.keys), max(rdf.keys)+self.binwidth, self.binwidth):
            for r in sorted(rdf.keys()):
                print('{: <10.3f}{: <10}'.format(r, rdf[r]), file=outf)
    @staticmethod
    def calc_distance(coord1, coord2, box_vectors, dim):
        def create_pbc_copies(coord, boxvectors, dim):
            pbc_copies = []
            if dim == 2:
                mirror_image_vectors = [
                                        (1, 0, 0), 
                                        (0, 1, 0), (0, -1, 0),
                                        (-1, 0, 0),
                                        ]
            elif dim == 3:
                mirror_image_vectors = [
                                        (1, 0, 0),
                                        (0, 1, 0), (0, -1, 0),
                                        (-1, 0, 0),
                                        (0, 0, 1 ),
                                         (0, 0, -1),
                                        ]
                #===============================================================
                # mirror_image_vectors = [
                #                         (1, 1, 0), (1, 0, 0), (1, -1, 0),
                #                         (0, 1, 0), (0, -1, 0),
                #                         (-1, 1, 0), (-1, 0, 0), (-1, -1, 0),
                #                         (1, 1, 1), (1, 0, 1), (1, -1, 1),
                #                         (0, 1, 1), (0, -1, 1), (0, 0, 1 ),
                #                         (-1, 1, 1), (-1, 0, 1), (-1, -1, 1),
                #                         (1, 1, -1), (1, 0, -1), (1, -1, -1),
                #                         (0, 1, -1), (0, -1, -1), (0, 0, -1),
                #                         (-1, 1, -1), (-1, 0, -1), (-1, -1, -1),
                #                         ]
                #===============================================================
            for vec in mirror_image_vectors:
                vec = np.array(vec[:dim])
                pbc_translation = vec*boxvectors[:dim]
                pbc_copy = coord[:dim] + pbc_translation
                pbc_copies.append(pbc_copy)
            return pbc_copies
        distances = []
        coords = [np.array(coord2[:dim])]
        rmax = min(box_vectors[:dim])/2
        #print("RMAX", rmax)
        for xyz in create_pbc_copies(coord2, box_vectors, dim): coords.append(xyz)
        for coord_partner in coords:
            distance = np.round(np.linalg.norm(coord_partner-coord1[:dim]), 5)
            if distance <= rmax:
                distances.append(distance)
            #distances.append(distance)
        #print("DISTANCES", distances)
        if len(distances) > 1:
            return min(distances)
            #raise ValueError("More than one distance smaller rmax")
        elif len(distances) == 0:
            return None
        return distances[0]

    def translate_selectionline(self, line):
        ''' translates line from vmd syntax to syntax usable for python and grofiles '''
        print("Translating:", line)
        choices = line.split(' and ')
        conditions = []
        for item in choices:
            choicelist = item.split(' ')
            description = choicelist[0]
            choice_items = choicelist[1:]
            conditions.append((description, choice_items))
            print(description, choice_items)
        return conditions

    def gather_coords(self, grofile, time):
        ''' !!! BEWARE: BAD CODE !!!
            Reads grofile and creates histogram from choices 
            Center of mass is actually calculated as center of geometry!
            Center of mass is always respective to resid
        '''
        def calculate_geometriccenter(coordinateinputlist, dim):
            #print("INPUTLIST", coordinateinputlist)
            geocenter = np.array([0., 0., 0.])
            #===================================================================
            # if self.dim == 3:
            #     geocenter = np.array([0, 0, 0])
            # elif self.dim == 2:
            #     geocenter = np.array([0, 0])
            #===================================================================
            for atomcoords in coordinateinputlist:
                #print("ATOMCOORDS", atomcoords)
                for ind, coord in enumerate(atomcoords):
                    #print("IND, COORD", ind, coord)
                    geocenter[ind]+=coord
                    #print("GEOCENTER", geocenter[ind])
            geocenter= geocenter/len(coordinateinputlist)
            #print(geocenter)
            #sys.exit()
            return np.array(geocenter)
        #description_ref, choice_ref = self.ref
        #description_sel, choice_sel = self.sel
        name2colnr = {
                    'resid':1,
                    'resname':2,
                    'atom':3,
                    'index':4,
                    }
        coordlist_ref = [[], []]
        coordlist_sel = [[], []]
        #========================RESID=======RESNAME======ATOMNAME=======INDEX===============X============Y============Z=======
        regexp = re.compile(r'^([\s,\d]{5})([\w,\s]{5})([\d,\w,\s]{5})([\s*\d+]{5})\s*(-?\d+\.\d+\s*-?\d+\.\d+\s*-?\d+\.\d+).*')
        regexp_boxline = re.compile(r'^\s+\d+\.\d+\s+\d+\.\d+\s+\d+\.\d+\s*$')
        if self.mode == 'SINGLE':
            for line in grofile:
                match = regexp.match(line)
                if match:
                    #atomname = match.group(3).strip()
                    resid = int(match.group(1).strip())
                    xyz = np.array([float(i) for i in match.group(5).split()])
                    #if xyz[2] < 4:
                        #continue
                    matched = True
                    for cond in self.ref:
                        if not match.group(name2colnr[cond[0]]).strip() in cond[1]:
                            matched = False
                            break
                    if matched:
                        if self.nchol != 'off':
                            neibs = self.neiblist[int(resid)][float(time)]
                            neibs_chol = [self.sysinfo.resid_to_lipid[N] for N in neibs].count('CHL1')
                            if neibs_chol == self.nchol:
                                if self.sysinfo.res_to_leaflet[resid]:
                                    coordlist_ref[0].append(xyz)
                                else:
                                    coordlist_ref[1].append(xyz)
                        else:
                            if self.dim == 2:
                                if self.sysinfo.res_to_leaflet[resid]:
                                    coordlist_ref[0].append(xyz)
                                else:
                                    coordlist_ref[1].append(xyz)
                            else:
                                coordlist_ref[0].append(xyz)
                    matched = True
                    for cond in self.sel:
                        if not match.group(name2colnr[cond[0]]).strip() in cond[1]:
                            matched = False
                            break
                    if matched:
                        if self.dim == 2:
                            if self.sysinfo.res_to_leaflet[resid]:
                                coordlist_sel[0].append(xyz)
                            else:
                                coordlist_sel[1].append(xyz)
                        else:
                            coordlist_sel[0].append(xyz)
                elif regexp_boxline.match(line):
                    boxvector = [float(i) for i in line.split()]
                    break
                else:
                    continue
        elif self.mode == 'COM':
            COMlist_ref = []
            COMlist_sel = []
            resid_old = 1
            #coordfilename = "z_coords_{}.xyz".format(time)
            #with open(coordfilename, "w") as cfile:
            for line in grofile:
                match = regexp.match(line)
                if match:
                    #atomname = match.group(3).strip()
                    resid = int(match.group(1).strip())
                    if resid != resid_old:
                        #print("NEW")
                        #print("COMLIST_REF", COMlist_ref)
                        if COMlist_ref:
                            geo_center_ref = calculate_geometriccenter(COMlist_ref, self.dim)
                            #print("GEOCENTER_REF", geo_center_ref)
                            #print("M{} {} {} {}".format(resid, *geo_center_ref), file=cfile)
                            COMlist_ref = []
                            if self.nchol != 'off':
                                neibs = self.neiblist[int(resid_old)][float(time)]
                                neibs_chol = [self.sysinfo.resid_to_lipid[N] for N in neibs].count('CHL1')
                                if neibs_chol == self.nchol:
                                    coordlist_ref[0].append(geo_center_ref)
                            else:
                                coordlist_ref[0].append(geo_center_ref)
                        if COMlist_sel:
                            geo_center_sel = calculate_geometriccenter(COMlist_sel, self.dim)
                            #print("GEO_SEL_Z", geo_center_sel[-1])
                            coordlist_sel[0].append(geo_center_sel)
                            COMlist_sel = []
                        #print("COORDLIST_REF", len(coordlist_ref[0]))
                        #print("COORDDLIST_SEL", len(coordlist_sel))
                    resid_old = resid
                    xyz = np.array([float(i) for i in match.group(5).split()])
                    if self.sysinfo.res_to_leaflet[resid]:
                        continue
                    matched = True
                    for cond, cond2 in zip(self.ref, self.ref2):
                        if not match.group(name2colnr[cond[0]]).strip() in cond[1] and not match.group(name2colnr[cond2[0]]).strip() in cond2[1]:
                            matched = False
                            break
                    if matched:
                        #print(resid)
                        COMlist_ref.append(xyz)
                    matched = True
                    for cond in self.sel:
                        if not match.group(name2colnr[cond[0]]).strip() in cond[1]:
                            matched = False
                            break
                    if matched:
                        COMlist_sel.append(xyz)
                elif regexp_boxline.match(line):
                    boxvector = [float(i) for i in line.split()]
                    break
                else:
                    continue
        return coordlist_sel, coordlist_ref, boxvector

    def process_frame(self, grofile, time):
        ''' Calculates the RDF for SEL around REF for ONE TIMEFRAME in a GROFILE '''
        binwidth = self.binwidth
        gofr = {}
        print(Time.ctime(), "Searching atoms.")
        coordlist_sel, coordlist_ref, boxvector = self.gather_coords(grofile, time)
        #print("COORDLISTS", len(coordlist_ref), len(coordlist_sel))
        rmax = min(boxvector[:2])/2
        print(Time.ctime(), "Creating histogram.")
        histog = {}
        histog_list = np.arange(0, rmax+2*binwidth, binwidth)
        for i in histog_list: histog[i] = 0 # Initialize
        for leaflet in range(1):
            for coord1 in coordlist_ref[leaflet]:
                for coord2 in coordlist_sel[leaflet]:
                    if (coord1 == coord2).all():
                        continue
                    distance = self.calc_distance(coord1, coord2, boxvector, self.dim)
                    if distance is None:
                        continue
                    interval = histog_list[bisect.bisect_right(histog_list, distance)]
                    histog[interval] += 1
        print(Time.ctime(), "Normalisation.")
        NumRefAtm = len(coordlist_ref[0])
        NumSelAtm = len(coordlist_sel[0])
        if self.dim == 2:
            Box = (2*rmax)**2                   # Box area
        elif self.dim == 3:
            Box = (2*rmax)**3
        if self.ref == self.sel:                # Box volume
            Rho_zero = (NumSelAtm-1)/Box        # Density of ideal gas
        else:
            Rho_zero = NumSelAtm/Box            # Density of ideal gas
        for r in histog_list:
            if self.dim == 2:
                bin_norm = np.pi*(((r+binwidth)**2)-(r**2))             # Area of bin
            elif self.dim == 3:
                bin_norm = (4/3)*(np.pi)*(((r+binwidth)**3)-(r**3))     # Volume of bin
            #AvgRefPerVol = NumRefAtm / Volume
            a = histog[r]/(NumRefAtm*bin_norm*Rho_zero)
            #c = histog[r]/(NumRefAtm*Volume_bin*(512/(3.24**3)))
            #b = histog[r]/(NumRefAtm*binwidth*Rho_zero*(4*np.pi*r**2))
            #a = histog[r]/((r**2)*Volume_bin*
            #b = histog[r]/(NumRefAtm*Volume_bin*Rho_zero)
            #c = histog[r]/(NumRefAtm*Rho_zero)
            #d = histog[r]/(NumRefAtm*Volume*AvgNumPart)
            #print("r: {}\tN(r): {}\nVolume: {}\nNumRefAtm: {}\tNumSelAtm: {}\nRho: {}".format(r, histog[r], Volume_bin, NumRefAtm, NumSelAtm, Rho_zero))
            #print("{: <10.3f}{: <10.3f}{: <10.3f}{: <10.3f}".format(a))
            #print(r, a)
            gofr[r] = a
        return gofr
    def calc_rdf(self):
        if self.grofile is None:
            grofile_output = produce_gro(self.sysinfo)
        else:
            grofile_output = self.grofile
        gofr_all_frames = {}
        with open(grofile_output,"r") as grofile:
            print(strftime("%H:%M:%S :", localtime()),"... read data from .gro-file ...")
            for line in grofile:
                if 't=' in line:
                    time = float(line[line.index('t=')+2:].strip())
                    print(time)
                    if float(self.sysinfo.t_end) < time:
                    #if time > 16000.:
                        break
                    else:
                        gofr_frame = self.process_frame(grofile, time)
                        for r in gofr_frame.keys():
                            try:
                                gofr_all_frames[r].append(gofr_frame[r])
                            except KeyError:
                                gofr_all_frames[r] = []
                                gofr_all_frames[r].append(gofr_frame[r])
        for r in sorted(gofr_all_frames.keys()):
            gval = np.mean(gofr_all_frames[r])
            gofr_all_frames[r] = gval
            #print(r, gval)
        return gofr_all_frames

def radialdistribution(systeminfo, ref, sel,):
    ''' Calculates the RDF of sel relative to ref. '''
    print("\n_____Calculating radial distribution function ____\n")
    print("Ref: {}\nSel: {}\n\n".format(ref, sel))
    os.makedirs(systeminfo.datapath+'/rdf', exist_ok=True)
    selectdict = {}
    if os.path.isfile("leaflet_assignment.dat"): 
        res_in_leaf = []
        with open("leaflet_assignment.dat", "r") as leas:
            leas.readline()
            for line in leas:
                cols = line.split()
                if int(cols[1]):
                    res_in_leaf.append(cols[0])
        select_one_leaflet = 'and resid '+' '.join(res_in_leaf)
    else:
        select_one_leaflet = 'and z<4'
    for selection in (('ref', ref), ('sel', sel)):
        selectstring = '{} {}'.format(selection[1], select_one_leaflet)
        select_fname = '{}/select{}_{}'.format(systeminfo.temppath, selection[0], selection[1])
        selectdict.update({selection[1]:select_fname})
        with open(select_fname, "w") as selfile:
            print("Selection for {} is:\n{}\n".format(selection[0], selectstring))
            print(selectstring, file=selfile)
    outputfile = '{}/rdf/rdf_{}-{}.xvg'.format(systeminfo.datapath, ref, sel).replace(" ", "")
    outputfile_cn = '{}/rdf/nr_{}-{}.xvg'.format(systeminfo.datapath, ref, sel).replace(" ", "")
    g_rdf_arglist = [gmx_exec, 'rdf', '-xy', '-xvg', 'none',
                     '-f', systeminfo.trjpath, '-s', systeminfo.tprpath,
                     '-o', outputfile, '-ref', '-sf', selectdict[ref],
                     '-sel', '-sf', selectdict[sel], '-cn', outputfile_cn,
                      ]
    out, err = exec_gromacs(g_rdf_arglist)
    rdf_log = 'rdf_{}-{}.log'.format(ref, sel).replace(" ", "")
    with open(rdf_log, "w") as logfile:
        print('STDERR:\n{}\n\nSTDOUT:\n{}'\
                .format(err.decode(), out.decode()), file=logfile)

def calculate_distance(self):
    print("\n____Calculating distances____\n")
    neiblist=self.get_neighbor_dict()
    os.makedirs(self.datapath+'/distcalc', exist_ok=True)
    for i in range(1,self.NUMBEROFPARTICLES+1):
        all_N_of_res=list(set([neibs for t in neiblist[i].keys() for neibs in neiblist[i][t]]))
        
        print("Working on lipid ",i,'...',end='\r')
        neighbors=[str(x) for x in all_N_of_res]
        host=str(i)
        selectstring='resid '+host+' and name P O3;\n'
        for x in neighbors:
            selectstring+=''.join(['resid ',x,' and name P O3;\n'])
        selectionfile=''.join([self.temppath,'/seldist'])
        outputfile=''.join([self.datapath,'/distcalc','/dist_to_hostresid',str(i),'.xvg'])
        with open(selectionfile,"w") as sf:
            sf.write(selectstring)    
        g_pairdist_arglist=[gmx_exec,'-nobackup','pairdist','-f',self.trjpath,\
                            '-s',self.tprpath,'-sf',selectionfile,'-o',outputfile,\
                            ]
        out, err = self.exec_gromacs(g_pairdist_arglist)
        with open("gmx_pairdist.log","w") as logfile:
            logfile.write(err.decode())
            logfile.write(150*'_')
            logfile.write(out.decode())
            logfile.write(150*'_')
    with open("all_distances.dat","w") as all_distances:
        # print("Time \t Host \t Neib \t distance/nm",file=all_distances)
        print("Time \t Host \t Neib \t distance/nm",file=all_distances)
        for resid in range(1,self.NUMBEROFPARTICLES+1):
            
            #print("Working on residue {}".format(resid), end="\r")
            distancefile = self.datapath+'/distcalc/dist_to_hostresid'+str(resid)+'.xvg'
            with open(distancefile,"r") as dfile:
                res_to_row = {}
                for distanceline in dfile: #folderlayout is: time distance_resHost_resNeib ...
                    distanceline_cols = distanceline.split()
                    if '@ s' in distanceline:                     #creating a dict to know which column(energies) belong to which residue
                        rownumber = int(distanceline_cols[1][1:])+1                 #time is at row 0 !
                        resnumber = distanceline_cols[4]
                        res_to_row.update({resnumber:rownumber})
                    if '@' not in distanceline and '#' not in distanceline:     #pick correct energies from energyfile and print
                        time_distanceline = float(distanceline_cols[0])
                        neighbors_are = neiblist[resid][float(time_distanceline)]
                        for neib in neighbors_are:
                            distance = distanceline_cols[res_to_row[str(neib)]]
                            print("{} \t {} \t {} \t {}".format(time_distanceline, resid, neib, distance), file=all_distances)

class Neighbors():
    ''' All calculations regarding lipid neighbors '''

    def __init__(self, systeminfo):
        self.mysystem = systeminfo
        self.cutoff = systeminfo.cutoff
        self.resnames = ' '.join(systeminfo.molecules)
        self.head_atomnames = ' '.join(lipidmolecules.central_atom_of.values())


    def create_selectionfile_neighborsearch(self, resid, refatoms='P'):
        filename="{}/neighbors_of_residue{}".format(self.mysystem.temppath, resid)
        with open(filename,"w") as selection:
            if refatoms == 'P':
                print(\
                    'host =  resid {0} and (name {2});\n'
                    'neibs = (resname {3} and name {2} and not host) and within {1} of host;\n'
                    'neibs;'\
                    .format(resid, self.cutoff, self.head_atomnames, self.resnames), file=selection)
                    #'host =  resid {0} and (name P O3);\n'
                    #'allOAtoms = resname CHL1 and name O3 and not host;\n'
                    #'allPAtoms = resname DPPC DUPC and name P and not host;\n'
                    #'neibOs = allOAtoms and within {1} of host;\n'
                    #'neibPs = allPAtoms and within {1} of host;\n'
                    #'neibs = neibPs or neibOs;\n'
                    #'neibs;'\
                    #.format(resid, self.cutoff,), file=selection)
            elif refatoms == 'bothtails':
                print(
                    'host = resid {0} and name C34 C24 O3;\n'
                    'allOAtoms = resname CHL1 and name O3 and not host;\n'
                    'allTail1Atoms = resname DPPC DUPC and name C34 and not host;\n'
                    'allTail2Atoms = resname DPPC DUPC and name C24 and not host;\n'
                    'neibOs = allOAtoms and within {1} of host;\n'
                    'neibTail1 = allTail1Atoms and within {1} of host;\n'
                    'neibTail2 = allTail2Atoms and within {1} of host;\n'
                    'neibs = neibOs or neibTail1 or neibTail2;\n'
                    'neibs;'\
                    .format(resid, self.cutoff), file=selection)
            elif refatoms == 'onetail':
                raise ValueError("It's not working like this, as the index is compared and so no lipid is selected")
                #print(
                #    'host = resid {0} and name C34 C24 O3;\n'
                #    'allOAtoms = resname CHL1 and name O3 and not host;\n'
                #    'allTail1Atoms = resname DPPC DUPC and name C34 and not host;\n'
                #    'allTail2Atoms = resname DPPC DUPC and name C24 and not host;\n'
                #    'neibOs = allOAtoms and within {1} of host;\n'
                #    'neibTail1 = allTail1Atoms and within {1} of host;\n'
                #    'neibTail2 = allTail2Atoms and within {1} of host;\n'
                #    'neibs = neibOs or (neibTail1 and neibTail2);\n'
                #    'neibs;'\
                #    .format(resid, self.cutoff), file=selection)
            else:
                raise ValueError("Wrong input for refatoms with: {}".format(refatoms))
        return filename

    def get_neighbor_dict(self, neighbor_filename='neighbor_info', verbose='off'):
        ''' Returns a list of all neighbors being in the
         cutoff distance at least once in the trajectory.
        Neighborfile is !required! and is output of 
        "determine_neighbors()" '''
        neighborfile = neighbor_filename
        neibdict = {}
        with open(neighborfile,"r") as neibmap:
            neibmap.readline() 
            for line in neibmap:
                cols = line.split()
                resid = int(cols[0])
                time = float(cols[1])
                if resid not in neibdict.keys():
                    neibdict.update({resid:{}})
                try:
                    neiblist = [int(x) for x in cols[3].split(',')]
                    neiblist = list(set(neiblist))
                    neibdict[resid].update({time:neiblist})
                except IndexError:
                    neibdict[resid].update({time:[]})
                    if verbose == 'on':
                        print("No neighbors of residue {} at time {}.".format(cols[0],cols[1])) 
        return neibdict

    def create_indexfile(self):
        ''' Creates an indexfile containing all indices of atom of each residue in system (resid_X) and all indices of all atoms in system.   '''
        print("\n_____Creating index file____\n")
        #OUTPUT IS:    resindex_all.ndx     | in the cwd!
        resindex_all = open("resindex_all.ndx","w")
        for mol in range(1, self.mysystem.NUMBEROFMOLECULES+1):
        #for mol in range(1,2):
            print("Working on residue {}".format(mol), end='\r')
            selectionfile = self.mysystem.temppath+'/tmp_selectionfile'
            with open(selectionfile,"w") as sf:
                lipidtype = self.mysystem.resid_to_lipid[mol]
                if lipidtype  not in lipidmolecules.sterols: 
                    tailhalf12_l = [] ### to get half the tails
                    tailhalf22_l = []
                    for molpart in lipidmolecules.tail_atoms_of[lipidtype]:
                        tailhalf12_l.extend(molpart[:len(molpart)//2])
                        tailhalf22_l.extend(molpart[len(molpart)//2:])
                    tailhalf12 = ' '.join(tailhalf12_l)
                    tailhalf22 = ' '.join(tailhalf22_l)
                    #lentailbyfour=len(self.tail_atoms_of[lipidtype])//4                 ### Problem for the half tail t12/t22:
                    #tailhalf12=' '.join(self.tail_atoms_of[lipidtype][:lentailbyfour]+self.tail_atoms_of[lipidtype][2*lentailbyfour:3*lentailbyfour])
                    #tailhalf22=' '.join(self.tail_atoms_of[lipidtype][lentailbyfour:2*lentailbyfour]+self.tail_atoms_of[lipidtype][3*lentailbyfour:])

                    tailatomlist = lipidmolecules.tail_atoms_of[lipidtype]
                    tailatoms = [x for index in tailatomlist for x in index] ##unpacking
                    headatoms = lipidmolecules.head_atoms_of[lipidtype]
                    
                    methylatomslists = []
                    for tailindex in range(len(lipidmolecules.tailcarbons_of[lipidtype])):
                        handlinghydr = [iter(lipidmolecules.tailhydr_of[lipidtype][tailindex])]*2
                        methylgrouplist = [i for i in zip(lipidmolecules.tailcarbons_of[lipidtype][tailindex], zip(*handlinghydr))]# Getting a list of tuples like [C1,(H1,H2),....]
                        methylgrouplist = [i for tup in methylgrouplist for i in tup]# Unpacking this list
                        methylgrouplist_unp = []
                        for particle in methylgrouplist:# get rid of hydr tuples
                            if isinstance(particle, tuple):
                                methylgrouplist_unp.extend(particle)
                            else:
                                methylgrouplist_unp.append(particle)
                        methylatomslists.append(methylgrouplist_unp)
                    methylgroups = [[methylatomslists[0][i:i+3]\
                                    +methylatomslists[0][i+3:i+6]\
                                    +methylatomslists[1][i:i+3]\
                                    +methylatomslists[1][i+3:i+6]]\
                               for i in range(0, len(methylatomslists[0])-1, 6)]
                    methylstrings = [' '.join(t) for i in methylgroups for t in i]
                    selprefixes = [("", r'".*"'),\
                                   ("h", ' '.join(headatoms)),\
                                   ("t", ' '.join(tailatoms)),\
                                   ("t12", tailhalf12),\
                                   ("t22", tailhalf22),\
                                   ("C0", methylstrings[0]),\
                                   ("C1", methylstrings[1]),\
                                   ("C2", methylstrings[2]),\
                                   ("C3", methylstrings[3]),\
                                   ("C4", methylstrings[4]),\
                                   ("C5", methylstrings[5]),\
                                   ("C6", methylstrings[6]),\
                                   ]
                    selectionlist = []
                    for item in selprefixes:
                        if item[0] == "":
                            selectionlist += ['resid_{0}=resid {0} and resname {1} and name {2};\n'\
                                             .format(str(mol), lipidtype, item[1])]
                        else:
                            selectionlist += ['resid_{0}_{1}=resid {1} and resname {2} and name {3};\n'\
                                             .format(item[0], str(mol), lipidtype, item[1])]
                    
                    lastlineitems = ['resid_{};\n'.format(str(mol))]
                    lastlineitems += ['resid_{}_{};\n'.format(item[0], str(mol)) for item in selprefixes if item[0] != ""]
                    selectionlist += ''.join(lastlineitems)
                    selectionstring=''.join(selectionlist)
                    #selectionlist += [''.join(["resid_h_",\
                    #               str(i),"=resid ",str(i)," and resname ",\
                    #               lipidtype," and name ",' '.join(headatoms),\
                    #               ";\n"]),\
                    #               ]
                    print(selectionstring)
                    sf.write(selectionstring)
                elif self.mysystem.resid_to_lipid[mol] in lipidmolecules.sterols:
                    selectionstring = 'resid_{0}=resid {0}  and resname {1};\nresid_{0};'.format(str(mol), self.mysystem.resid_to_lipid[mol])
                    sf.write(selectionstring)
            outputindex = self.mysystem.indexpath+"/resid_"+str(mol)+".ndx"
            gmx_select_arglist = [gmx_exec, 'select', '-s', self.mysystem.gropath, '-sf',\
                                  selectionfile, '-on', outputindex,\
                                  ]
            out, err = exec_gromacs(gmx_select_arglist)
            with open("gmx_select.log","w") as logfile: 
                logfile.write(err.decode())
                logfile.write(150*'_')
                logfile.write(out.decode())
                logfile.write(150*'_')
            ##append resid.ndx to resindex_all.ndx
            with open(outputindex, "r") as output_index:
                filecontent = output_index.readlines()
                resindex_all.write(''.join(filecontent)+'\n\n')
        ### To have whole system indices in one group    
        make_ndx_output = self.mysystem.temppath+'/make_ndx_system.ndx' 
        gmx_make_ndx_arglist = [gmx_exec, 'make_ndx', '-f', self.mysystem.gropath, '-o', make_ndx_output] 
        inp_str = b'keep 0\nq\n'  
        out, err = exec_gromacs(gmx_make_ndx_arglist, inp_str) 
        with open("gmx_make_ndx.log","w") as logfile, open(make_ndx_output, "r") as output:
            logfile.write(err.decode())
            logfile.write(150*'_')
            logfile.write(out.decode())
            logfile.write(150*'_')
            filecontent = output.readlines()
            resindex_all.write(''.join(filecontent)+"\n\n")
        resindex_all.close() 

    def determine_neighbors(self, refatoms='P', overwrite=True):
        ''' Creates "neighbor_info" containing all information on lipid arrangement '''
        print("\n____Determining neighbors____\n")
        os.makedirs(self.mysystem.datapath+'/neighborfiles', exist_ok=True)
        with open("neighbor_info", "w") as outfile:
            outfile.write('Resid \t Time \t Number_of_neighbors \t List_of_Neighbors \n')
            for residue in range(1, self.mysystem.NUMBEROFMOLECULES+1):
                print(". . . Working on residue: {} . . .".format(residue), end="\r")
                selectionfile = self.create_selectionfile_neighborsearch(residue, refatoms=refatoms)
                indexoutput = '{}/neighbors_of_residue{}.ndx'.format(self.mysystem.indexpath, residue)
                datafileoutput = '{}/neighborfiles/neighbors_of_residue{}.dat'.format(self.mysystem.datapath, residue)
                if os.path.isfile(datafileoutput) and overwrite == False:
                    print("Neighbor file of residue {} already exists. Skipping.".format(residue))
                else:
                    cmdlist=[\
                        gmx_exec,'select','-s', self.mysystem.tprpath, '-f', self.mysystem.trjpath,\
                        '-sf',selectionfile,'-on',indexoutput,'-oi',datafileoutput,\
                        '-b',str(self.mysystem.t_start), '-e', str(self.mysystem.t_end),\
                        '-dt', str(self.mysystem.dt),\
                        ]
                    out, err = exec_gromacs(cmdlist)
                    with open("gmx_select.log","w") as logfile:
                        logfile.write(err.decode())
                        logfile.write(150*'_')
                        logfile.write(out.decode())
                        logfile.write(150*'_')
                with open(datafileoutput,"r") as datfile:
                    for line in datfile:
                        cols = line.split()
                        time = cols.pop(0)
                        nneibs = cols.pop(0)
                        neibindeces = [int(x) for x in cols]
                        #print("index: ", neibindeces)
                        #print("dict: ", self.mysystem.index_to_resid.keys())
                        neibresid = [self.mysystem.index_to_resid[x] for x in neibindeces]
                        residlist = ','.join([str(x) for x in neibresid])
                        print('{} \t {} \t {} \t {}'.format(residue,time,nneibs,residlist), file=outfile)

class Energy():
    ''' All about energy calculation '''

    DENOMINATOR = 40

    def __init__(self, systeminfo, parts, neighborfile_name='neighbor_info', resindex_all='resindex_all', overwrite=True):
        knownparts = ['complete', 'head-tail', 'head-tailhalfs', 'carbons']
        if parts not in knownparts:
            raise ValueError("Part keyword specified is not known.")
        self.neiblist = Neighbors(systeminfo).get_neighbor_dict()
        self.resindex_all = resindex_all
        self.overwrite = overwrite
        self.groupblocks = ()
        self.mysystem = systeminfo
        if parts == 'complete':
            self.molparts = ["resid_"]
            self.parts = ''
            self.denominator = self.DENOMINATOR
            self.molparts_short = [""]
            if neighborfile_name != 'neighbor_info':
                self.all_energies = 'all_energies_{}.dat'.format(neighborfile_name)
            else:
                self.all_energies = 'all_energies.dat'
            #self.interactions = ['']
        elif parts == 'head-tail':
            self.molparts = ["resid_h_", "resid_t_"]
            self.parts = parts
            self.denominator = int(self.DENOMINATOR/2)
            self.molparts_short = ["h_", "t_"]
            #self.interactions = ['head-tail', 'head-head', 'tail-tail']
            self.all_energies='all_energies_headtail.dat'
        elif parts == 'head-tailhalfs':
            self.molparts = ["resid_h_", "resid_t12_", "resid_t22_"]
            self.parts = parts
            self.denominator = int(self.DENOMINATOR/4)
            self.molparts_short = ["h_","t12_","t22_"]
            #self.interactions = ['head-tail12', 'tail12-tail12', 'head-tail22', 'tail22-tail22']
            self.all_energies = 'all_energies_headtailhalfs.dat'
        elif parts == 'carbons':
            self.molparts = ['resid_C{}_'.format(i) for i in range(7)]
            self.parts = parts
            self.denominator = int(self.DENOMINATOR/10)
            self.molparts_short = ['C{}_'.format(i) for i in range(7)]
            #self.interactions = ['C{0}-C{0}'.format(i) for i in range(7)]\
            #                    +['C{0}-C{1}'.format(i, i+1) for i in range(6)]\
            #                    +['C{0}-C{1}'.format(i, i-1) for i in range(1,7)]
            self.all_energies = "all_energies_carbons.dat"
        print('\n Calculating for energygroups:', self.molparts)

    def run_calculation(self, startres=-1, endres=-1):
        ''' Runs an energy calculation with settings from Energy() instance.
            For each residue the energy to all neighbors seen during MD is calculated
            and written to .edr files.
            Procedure is as follows:
            1. The neighbors are divided into fragments ("groupfragments")
            2. For each fragment:
                an mdp file is created (create_MDP)
                a tpr file is generated (create_TPR)
            3. The actual mdrun -rerun is performed (do_Energyrun)
            4. .xvg tables are generate from .edr files
            
        '''
        print('''\n____Rerunning MD for energyfiles,
         creating xvgtables with relevant energies.____\n
         Caution mdp-file must not have energy_grps indicated!\n''')
        if startres == -1:
            startres = 1
        if endres == -1:
            endres = self.mysystem.NUMBEROFMOLECULES
        for res in range(startres,endres+1): 
            print('\n',strftime("%H:%M:%S :", localtime()),'Working on lipid '+str(res)+'...')
            all_neibs_of_res = list(set([neibs for t in self.neiblist[res].keys() for neibs in self.neiblist[res][t]]))
            N_neibs = len(all_neibs_of_res)
            # ''' Dividing calculations to separate neighbors considered in each energyrun
            # as max(energygroups) ≈ 64 | not sure about exact number
            if N_neibs % self.denominator == 0:
                number_of_groupfragments = (N_neibs//self.denominator)
            else:
                number_of_groupfragments = (N_neibs//self.denominator)+1
            print("Needing {} energy run(s)".format(number_of_groupfragments)) 
            for groupfragment in range(number_of_groupfragments):
                print("On fragment",groupfragment)
                
                g_energy_output = ''.join([\
                    self.mysystem.energypath, '/xvgtables/energies_residue',\
                    str(res), '_', str(groupfragment), self.parts, '.xvg',\
                    ])
                groupblockstart = groupfragment*self.denominator
                groupblockend = (groupfragment+1)*self.denominator
                self.groupblocks = (groupblockstart, groupblockend)
                # File in-/outputs
                groupfragment=str(groupfragment) 
                mdpout = ''.join([self.mysystem.energypath, '/mdpfiles/energy_mdp_recalc_resid', str(res), '_', groupfragment, self.parts, '.mdp'])
                tprout = ''.join([self.mysystem.energypath, '/tprfiles/mdrerun_resid', str(res), '_', groupfragment, self.parts, '.tpr'])
                energyf_output = ''.join([self.mysystem.energypath, '/edrfiles/energyfile_resid', str(res), '_'+groupfragment, self.parts, '.edr'])
                xvg_out = ''.join([self.mysystem.energypath, '/xvgtables/energies_residue', str(res), '_', groupfragment, self.parts, '.xvg'])
                energygroups = self.gather_energygroups(res, all_neibs_of_res)
                relev_energies = self.get_relev_energies(res, all_neibs_of_res)
                # Run functions
                self.create_MDP(mdpout, energygroups)
                self.create_TPR(mdpout, tprout)
                if os.path.isfile(energyf_output) and self.overwrite == False:
                    print("Edrfile for lipid {} part {} already exists. Will skip this calculation.".format(res, groupfragment))
                else:
                    self.do_Energyrun(res, groupfragment, tprout, energyf_output)
                if os.path.isfile(g_energy_output) and self.overwrite == False:
                    print("Xvgtable for lipid {} part {} already exists. Will skip this calculation.".format(res, groupfragment))
                else:
                    self.write_XVG(energyf_output, tprout, relev_energies, xvg_out)
        return 'Done'

    def selfinteractions_edr_to_xvg(self):
        ''' Extracts all self interaction energy values from .edr files using gmx energy '''
        for res in range(1,self.mysystem.NUMBEROFMOLECULES+1):
            relev_energies = self.get_relev_self_interaction(res)
            tprout = ''.join([self.mysystem.energypath, '/tprfiles/mdrerun_resid', str(res), '_', '0', self.parts, '.tpr'])
            energyf_output = ''.join([self.mysystem.energypath, '/edrfiles/energyfile_resid', str(res), '_'+'0', self.parts, '.edr'])
            xvg_out = ''.join([self.mysystem.energypath, '/xvgtables/energies_residue', str(res), '_selfinteraction', self.parts, '.xvg'])
            self.write_XVG(energyf_output, tprout, relev_energies, xvg_out)
    def selfinteractions_xvg_to_dat(self):
        ''' Extracts all self interaction energy entries from xvg files
            and writes them to "selfinteractions.dat"
        '''
        with open("selfinteractions.dat", "w") as energyoutput:
            print(\
                  '{: <10}{: <10}'
                  '{: <20}{: <20}{: <20}{: <20}{: <20}{: <20}{: <20}'\
                  .format("Time", "Lipid",
                          "Etot", "VdWSR", "CoulSR", "VdW14", "Coul14", "VdWtot", "Coultot", ),
                  file=energyoutput)
            for resid in range(1, self.mysystem.NUMBEROFMOLECULES+1):
                #print("Working on residue {}".format(resid), end="\r")
                #residtype = self.mysystem.resid_to_lipid[resid]
                xvg_out = ''.join([self.mysystem.energypath, '/xvgtables/energies_residue', str(resid), '_selfinteraction', self.parts, '.xvg'])
                with open(xvg_out,"r") as xvgfile:
                    res_to_row = {}
                    for energyline in xvgfile: #folderlayout is: time Coul_resHost_resNeib LJ_resHost_resNeib ...
                        energyline_cols = energyline.split()
                        if '@ s' in energyline:                     #creating a dict to know which column(energies) belong to which residue
                            print(energyline_cols)
                            row = int(energyline_cols[1][1:])+1                 #time is at row 0 !
                            host = energyline_cols[3].split("resid_")[1][:-1]
                            energytype = energyline_cols[3].split(":")[0][1:]
                            #print("Host: {} Type: {} Row: {}".format(host, energytype, row))
                            res_to_row.update({(energytype, host):row})
                        elif '@' not in energyline and '#' not in energyline:     #pick correct energies from energyfile and print
                            time = float(energyline_cols[0])
                            #print("TIME", time)
                            if time % self.mysystem.dt != 0:
                                continue
                            try:
                                vdw_sr = energyline_cols[res_to_row[('LJ-SR', str(host))]]
                                vdw_14 = energyline_cols[res_to_row[('LJ-14', str(host))]]
                                coul_sr = energyline_cols[res_to_row[('Coul-SR', str(host))]]
                                coul_14 = energyline_cols[res_to_row[('Coul-14', str(host))]]
                                vdw_tot = float(vdw_14) + float(vdw_sr)
                                coul_tot = float(coul_14) + float(coul_sr)
                            except KeyError:
                                continue
                            Etot = float(vdw_sr)+float(coul_sr)+float(vdw_14)+float(coul_14)
                            print(\
                                  '{: <10}{: <10}{: <20.5f}'
                                  '{: <20}{: <20}{: <20}{: <20}{: <20.5f}{: <20.5f}'
                                  .format(time, resid,  Etot,
                                          vdw_sr, coul_sr, vdw_14, coul_14, vdw_tot, coul_tot,),\
                                  file=energyoutput)

    def gather_energygroups(self, res, all_neibs_of_res):
        ''' Set which part of molecule should be considered '''
        energygroup_indeces=[res]+all_neibs_of_res[self.groupblocks[0]:self.groupblocks[1]]
        energygroup_list=[]
        for index in energygroup_indeces:
            if self.mysystem.resid_to_lipid[index] in lipidmolecules.sterols:
                energygroup_list.append(''.join(["resid_",str(index)]))
            else:
                for part in self.molparts:
                    energygroup_list.append(''.join([part,str(index)]))
        energygroup_string=' '.join(energygroup_list)
        return energygroup_string

    def get_relev_energies(self, res, all_neibs_of_res):
        ''' Returns string that describes all entries
            needed to be extracted from energy file using gmx energy 
            This version is for lipid-lipid interaction
            for self interaction search function "get_relev_self_interaction"
        '''
        Etypes=["Coul-SR:", "LJ-SR:"]
        energyselection=[]
        for interaction in Etypes:
            counterhost=0 #for the cholesterol as it has just 1 molpart
            for parthost in self.molparts:
                if self.mysystem.resid_to_lipid[res]=='CHL1' and counterhost==0:
                    parthost="resid_"
                    counterhost+=1
                elif self.mysystem.resid_to_lipid[res]=='CHL1' and counterhost!=0:
                    continue
                for neib in all_neibs_of_res[self.groupblocks[0]:self.groupblocks[1]]:
                    counterneib=0
                    for partneib in self.molparts:
                        if self.mysystem.resid_to_lipid[neib]=='CHL1' and counterneib==0:
                            partneib='resid_'
                            counterneib+=1
                        elif self.mysystem.resid_to_lipid[neib]=='CHL1' and counterneib!=0:
                            continue
                        energyselection.append(''.join([interaction,parthost,str(res),"-",partneib,str(neib)]))
        all_relev_energies='\n'.join(energyselection+['\n'])
        return all_relev_energies

    def get_relev_self_interaction(self, res):
        ''' Returns string that describes all entries
            needed to be extracted from energy file using gmx energy 
            This version is for lipid self interaction
        '''
        Etypes=["Coul-SR:", "LJ-SR:", "Coul-14:", "LJ-14:"]
        energyselection=[]
        for interaction in Etypes:
            for parthost in self.molparts:
                energyselection.append(''.join([interaction,parthost,str(res),"-",parthost,str(res)]))
        all_relev_energies='\n'.join(energyselection+['\n'])
        return all_relev_energies

    def create_MDP(self, mdpout: str, energygroups: str):
        ''' Create Mdpfile '''
        os.makedirs(self.mysystem.energypath+'/mdpfiles', exist_ok=True)
        with open(mdpout,"w") as mdpfile_rerun:#, open(mdp_raw,"r") as mdpfile_raw:
            raw_mdp =[x.strip() for x in '''
            integrator              = md               
            dt                      = 0.002            
            nsteps                  =                  
            nstlog                  = 100000           
            nstxout                 = 0                
            nstvout                 = 0                
            nstfout                 = 0                
            nstcalcenergy           = 1000             
            nstenergy               = 100              
            cutoff-scheme           = Verlet           
            nstlist                 = 20               
            rlist                   = 1.2              
            coulombtype             = pme              
            rcoulomb                = 1.2              
            vdwtype                 = Cut-off          
            vdw-modifier            = Force-switch     
            rvdw_switch             = 1.0              
            rvdw                    = 1.2              
            tcoupl                  = Nose-Hoover      
            tau_t                   = 1.0               
            tc-grps                 = System        
            pcoupl                  = Parrinello-Rahman
            pcoupltype              = semiisotropic    
            tau_p                   = 5.0
            compressibility         = 4.5e-5  4.5e-5   
            ref_p                   = 1.0     1.0      
            constraints             = h-bonds          
            constraint_algorithm    = LINCS            
            continuation            = yes               
            nstcomm                 = 100              
            comm_mode               = linear           
            refcoord_scaling        = com              
            '''.split('\n')]
            raw_mdp.append('ref_t = '+str(self.mysystem.temperature))
            #mdp_raw_content = mdpfile_raw.readlines()
            #if not mdp_raw_content:
            #    raise EOFError(".mdp-file is empty")
            energygrpline = ''.join(['energygrps\t\t\t=', energygroups, '\n'])
            raw_mdp.append(energygrpline)
            mdpfile_rerun.write('\n'.join(raw_mdp)+'\n')

    def create_TPR(self, mdpoutfile: str, tprout: str):
        ''' Create TPRFILE with GROMPP '''
        print(strftime("%H:%M:%S :", localtime()), '...Creating .tpr-file...')
        os.makedirs(self.mysystem.energypath+'/tprfiles', exist_ok=True)
        grompp_arglist=[gmx_exec, 'grompp', '-f', mdpoutfile, '-p',\
                        self.mysystem.toppath, '-c', self.mysystem.gropath, '-o', tprout,\
                        '-n', self.resindex_all, '-po', mdpoutfile\
                        ]
        out, err = exec_gromacs(grompp_arglist)
        with open("gmx_grompp.log","a") as logfile:
            logfile.write(err.decode())
            logfile.write(100*'_')
            logfile.write(out.decode())
            logfile.write(100*'_')

    def do_Energyrun(self, res, groupfragment, tprrerun_in, energyf_out):
        ''' Create .edr ENERGYFILE with mdrun -rerun '''
        print(strftime("%H:%M:%S :", localtime()), '...Rerunning trajectory for energy calculation...')
        os.makedirs(self.mysystem.energypath+'/edrfiles', exist_ok=True)
        os.makedirs(self.mysystem.energypath+'/logfiles', exist_ok=True)
        logoutput_file = self.mysystem.energypath+'/logfiles'+'/mdrerun_resid'+str(res)+self.parts+'frag'+groupfragment+'.log'
        trajout = 'EMPTY.trr' # As specified in mdpfile, !NO! .trr-file should be written
        mdrun_arglist = [gmx_exec, 'mdrun', '-s', tprrerun_in,'-rerun', self.mysystem.trjpath,\
                        '-e', energyf_out, '-o', trajout,'-g', logoutput_file,\
                        ]#,'-nt','8']
        out, err = exec_gromacs(mdrun_arglist)
        with open("gmx_mdrun.log","a") as logfile:
            logfile.write(err.decode())
            logfile.write(100*'_')
            logfile.write(out.decode())
            logfile.write(100*'_')

    def write_XVG(self, energyf_in, tprrerun_in, all_relev_energies, xvg_out):
        ''' Create XVG-TABLE with all relevant energies '''
        print(strftime("%H:%M:%S :", localtime()),'...Extracting all relevant energies from .edr file...')
        os.makedirs(self.mysystem.energypath+'/xvgtables', exist_ok=True)
        g_energy_arglist=[gmx_exec,'energy','-f',energyf_in,\
                          '-s', tprrerun_in,'-o', xvg_out,\
                          ]
        inp_str=all_relev_energies.encode()
        out, err = exec_gromacs(g_energy_arglist, inp_str)
        with open("gmx_energy.log","a") as logfile:
            logfile.write(err.decode())
            logfile.write(100*'_')
            logfile.write(out.decode())
            logfile.write(100*'_')

    def write_energyfile(self):
        ''' Creates files: "all_energies_<interaction>.dat '''
        print('____ Create energy file ____')
        with open(self.all_energies, "w") as energyoutput:
            print(\
                  '{: <10}{: <10}{: <10}{: <20}'
                  '{: <20}{: <20}{: <20}'\
                  .format("Time", "Host", "Neighbor", "Molparts",\
                                           "VdW", "Coul", "Etot"),\
                  file=energyoutput)
            for resid in range(1, self.mysystem.NUMBEROFMOLECULES+1):
                print("Working on residue {}".format(resid), end="\r")
                residtype = self.mysystem.resid_to_lipid[resid]
                all_neibs_of_res = list(set([neibs for t in self.neiblist[resid].keys() for neibs in self.neiblist[resid][t]]))
                n_neibs = len(all_neibs_of_res)
                if n_neibs % self.denominator == 0:
                    number_of_groupfragments = (n_neibs//self.denominator)
                else:
                    number_of_groupfragments = (n_neibs//self.denominator)+1
                #print(number_of_groupfragments)
                print("Nneibs: {}\tNfrags: {}\t".format(n_neibs, number_of_groupfragments))
                processed_neibs = []
                for part in range(number_of_groupfragments):
                    #print(part)
                    #groupblockstart = part*self.denominator
                    #groupblockend = (part+1)*self.denominator
                    #neighbors_part_are = all_neibs_of_res[groupblockstart:groupblockend]
                    #print("NEIBSPART", neighbors_part_are)
                    xvgfilename = self.mysystem.energypath+'/xvgtables/energies_residue'+str(resid)+'_'+str(part)+self.parts+'.xvg'
                    with open(xvgfilename,"r") as xvgfile:
                        #print("LOOKING IN FILE:", xvgfilename)
                        res_to_row = {}
                        for energyline in xvgfile: #folderlayout is: time Coul_resHost_resNeib LJ_resHost_resNeib ...
                            energyline_cols = energyline.split()
                            if '@ s' in energyline:                     #creating a dict to know which column(energies) belong to which residue
                                row = int(energyline_cols[1][1:])+1                 #time is at row 0 !
                                neib = energyline_cols[3].split("resid_")[2][:-1]
                                host = energyline_cols[3].split("resid_")[1][:-1]
                                energytype = energyline_cols[3].split("-")[0][1:]
                                if energytype == 'LJ':
                                    processed_neibs.append(int(neib))
                                res_to_row.update({(energytype, host, neib):row})
                            elif '@' not in energyline and '#' not in energyline:     #pick correct energies from energyfile and print
                                time = float(energyline_cols[0])
                                if time % self.mysystem.dt != 0:
                                    continue
                                #for neib in neighbors_part_are:
                                for neib in all_neibs_of_res:
                                    if self.mysystem.system == 'dppc_dupc_chol25' and ((int(host) == 372 and neib == 242) or (int(host) == 242 and neib == 372)):
                                        continue
                                    neibtype = self.mysystem.resid_to_lipid[neib]
                                    counterhost = 0
                                    for parthost in self.molparts:
                                        parthost = parthost[6:]
                                        if residtype == 'CHL1' and counterhost == 0:
                                            parthost = ''
                                            counterhost += 1
                                        elif residtype == 'CHL1' and counterhost != 0:
                                            continue
                                        counterneib = 0
                                        for partneib in self.molparts:
                                            partneib = partneib[6:]
                                            if neibtype == 'CHL1' and counterneib == 0:
                                                partneib = ''
                                                counterneib += 1
                                            elif neibtype == 'CHL1' and counterneib != 0:
                                                continue
                                            if parthost[:-1] == '':
                                                interhost = 'w'
                                            else:
                                                interhost = parthost[:-1]
                                            if partneib[:-1] == '':
                                                interneib = 'w'
                                            else:
                                                interneib = partneib[:-1]
                                            inter = ''.join([interhost, '_', interneib])
                                            #print(('LJ', parthost+str(resid), partneib+str(neib)))
                                            try:
                                                vdw = energyline_cols[res_to_row[('LJ', parthost+str(resid), partneib+str(neib))]]
                                                coul = energyline_cols[res_to_row[('Coul', parthost+str(resid), partneib+str(neib))]]
                                            except KeyError:
                                                continue
                                            Etot = float(vdw)+float(coul)
                                            print(\
                                                  '{: <10}{: <10}{: <10}{: <20}'
                                                  '{: <20}{: <20}{: <20.5f}'
                                                  .format(time, resid, neib, inter,\
                                                                            vdw, coul, Etot),\
                                                  file=energyoutput)
                #print("PROC NEIBS", processed_neibs)
                for pneib in processed_neibs:
                    #print("PNEIB, ALLNEIBS", pneib, all_neibs_of_res)
                    all_neibs_of_res.remove(pneib)
                if len(all_neibs_of_res):
                    print("Missing neighbour-ids", all_neibs_of_res)
                    raise ValueError('Not all neighbours found in xvgfile')

    def check_exist_xvgs(self):
        ''' Checks if all .xvg-files containing lipid interaction exist '''
        all_okay = True
        with open("missing_xvgfiles.info", "w") as inffile:
            print("#Files missing", file=inffile)
            for resid in range(1, self.mysystem.NUMBEROFMOLECULES+1):
                all_neibs_of_res = list(set([neibs for t in self.neiblist[resid].keys() for neibs in self.neiblist[resid][t]]))
                n_neibs = len(all_neibs_of_res)
                if n_neibs % self.denominator == 0:
                    number_of_groupfragments = (n_neibs//self.denominator)
                else:
                    number_of_groupfragments = (n_neibs//self.denominator)+1
                for part in range(number_of_groupfragments):
                    xvgfilename = self.mysystem.energypath+'/xvgtables/energies_residue'+str(resid)+'_'+str(part)+self.parts+'.xvg'
                    if not os.path.isfile(xvgfilename):
                        print(xvgfilename, file=inffile)
                        all_okay = False
        if not all_okay:
            print('THERE ARE FILES MISSING')
            #raise ValueError('THERE ARE FILES MISSING')
        return all_okay
