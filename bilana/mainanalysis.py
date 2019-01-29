''' Set of different tools for analysis
    Specify more here:
        -
 '''
import os
#import sys
import pprint
import re
import bisect
import logging

import numpy as np
import MDAnalysis as mda

from time import localtime, strftime
from bilana.common import AutoVivification, GRO_format
from bilana.energyfilecreator import read_scdinput
from bilana import lipidmolecules
from bilana import gromacstoolautomator as gmxauto
from bilana.systeminfo import SysInfo

pp = pprint.PrettyPrinter()


logger = logging.getLogger()
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
ch = logging.StreamHandler()
logger.setLevel(logging.INFO)
ch.setLevel(logging.INFO)
#logger.setLevel(logging.DEBUG)
#ch.setLevel(logging.DEBUG)
ch.setFormatter(formatter)
#logger.addHandler(ch)
logger.debug("Initialized logger with handers %s", logger.handlers)



def is_neighbor_in_leaflet(systeminfo_inst):
    ''' Searches for interleaflet neighborhood '''
    leaflet_assign = systeminfo_inst.res_to_leaflet
    #nlipids = systeminfo_inst.NUMBEROFMOLECULES
    neiblist = gmxauto.Neighbors(systeminfo_inst).get_neighbor_dict()
    host_has_interleafletneib = []
    for host in range(*systeminfo_inst.MOLRANGE):
        neibs_times = neiblist[host]
        host_leaflet = leaflet_assign[host]
        for t in neibs_times:
            for neib in neibs_times[t]:
                neib_leaflet = leaflet_assign[neib]
                if neib_leaflet != host_leaflet:
                    host_has_interleafletneib.append([host, neib])
    pp.pprint(host_has_interleafletneib)

class Scd():
    ''' All about calculating the lipid Scd order parameter '''

    def __init__(self, systeminfo):
        self.systeminfo = systeminfo
        self.atomlist = lipidmolecules.scd_tail_atoms_of
        self.neiblist = gmxauto.Neighbors(self.systeminfo).get_neighbor_dict()
        self.components = self.systeminfo.molecules
        #print(self.atomlist)
    def create_scdfile(self, grofilename=None, outputfile='scd_distribution.dat'):
        # Check if trajectory is already converted to .gro format
        if grofilename is None:
            grofilename = 'traj_complete.gro'
            grofilepath = ''.join([self.systeminfo.datapath, '/grofiles/'])
            gmxauto.produce_gro(self.systeminfo, grofilename=grofilename)
        else:
            grofilepath = ''
        # Start reading grofile and calculating scd
        with open(grofilepath+grofilename,"r") as grofile, open(outputfile, "w") as scdfile:
            components = self.systeminfo.molecules
            print("{: <12}{: <10}{: <7}{: <15}".format("Time", "Residue", "Type", "Scd") + (len(components)*'{: ^7}').format(*components), file=scdfile)
            resid_old = 0
            time = None
            lipidtype_old = ''
            coorddict = {}
            print(strftime("%H:%M:%S :", localtime()),"... read data from .gro-file ...")
            for line in grofile:
                # In gro each frame starts with time annotation t=
                if 't=' in line:
                    re_time = re.compile(r'.*t=\s+(\d+\.\d+).*')
                    time_tmp = float(re_time.match(line).group(1))
                    if time is None:
                        time = time_tmp
                    if float(self.systeminfo.t_end) < time:
                        break
                regmatch = GRO_format.regexp.match(line)

                # If matched: Line is normal atom indicator -> See common/GROM_format
                if float(self.systeminfo.t_start) <= time and regmatch is not None:
                    grps = regmatch.groups()
                    atom = grps[2].strip()
                    lipidtype = grps[1].strip()

                    # Skip molecules that are not in bilayer
                    if lipidtype[:-2] not in lipidmolecules.TAIL_ATOMS_OF.keys()\
                        and lipidtype not in lipidmolecules.STEROLS\
                        and lipidtype not in lipidmolecules.PROTEINS:
                        logger.debug("Skip: %s", lipidtype)
                        continue

                    # Initialize lipidtype_old (last lipidtype that was used for calculation)
                    if not lipidtype_old:
                        lipidtype_old = lipidtype

                    all_atmlst = [atm for atmlst in self.atomlist(lipidtype) for atm in atmlst]
                    if atom in all_atmlst:# and lipidtype in self.systeminfo.molecules:
                        resid = int(grps[0].strip())
                        logger.debug("Groups %s", grps)
                        if not resid_old:
                            logger.debug("Setting resid_old from %s to %s", resid_old, resid)
                            resid_old = resid
                        if resid != resid_old and coorddict:
                            # New entry for molecule started in gro, now calculating for resid_old
                            logger.debug("Resid / resid_old: %s / %s", resid, resid_old)
                            logger.debug("Atomlist %s", self.atomlist(lipidtype_old))
                            logger.debug("coordinates: %s", coorddict)
                            scd_value = self.scd_of_res(coorddict, self.atomlist(lipidtype_old))
                            coorddict = {}
                            # Counting the number of neighbors with type
                            neibs = self.neiblist[resid_old][float(time)]
                            neib_comp_list = []
                            for lip in self.components:
                                ncomp = [self.systeminfo.resid_to_lipid[N] for N in neibs].count(lip)
                                neib_comp_list.append(ncomp)

                            # Print results to output file
                            print("{: <12.2f}{: <10}{: <7}{: <15.8}".format(
                                time, resid_old, lipidtype_old, scd_value)\
                                + (len(neib_comp_list)*'{: ^7}').format(*neib_comp_list),
                                file=scdfile)

                            # Reset all indicators
                            resid_old = resid
                            time = time_tmp
                            lipidtype_old = lipidtype
                        # Read new coordinates of <resid> and <atom>
                        coordinates = [float(x) for x in grps[4:7]]
                        coorddict[atom] = coordinates
                    else: # Not matched
                        continue
            # To calculate the last value
            scd_value = self.scd_of_res(coorddict, self.atomlist(lipidtype_old))
            coorddict = {}
            neibs = self.neiblist[resid_old][float(time)]
            neib_comp_list = []
            for lip in self.components:
                ncomp = [self.systeminfo.resid_to_lipid[N] for N in neibs].count(lip)
                neib_comp_list.append(ncomp)
            print("{: <12.2f}{: <10}{: <7}{: <15.8}".format(
                time, resid_old, lipidtype_old, scd_value)\
                + (len(neib_comp_list)*'{: ^7}').format(*neib_comp_list),
                file=scdfile)
        print(strftime("%H:%M:%S :", localtime()),"Finished reading.")
        return

    def scd_of_res(self, coorddict, atomlist):#, neibstraightness=0,neiblist=None):
        #print(atomlist)
        scds_of_atoms = []
        scds_of_tails = []
        #scds_of_tails_corrected=[]
        for tail in atomlist:
            #print(tail)
            for atomindex in range(len(tail)-1): ### -1 because last res is not taken (taking index len() implies "range+1")
                atm1, atm2 = tail[atomindex], tail[atomindex+1]    ### Attention: In the tail list only Scd-specific (every 2nd) atom is included!! Thus: atomindex-atomindex+1
                #print("CALC FOR ", atm1, atm2)
                coords_atm1, coords_atm2 = np.array(coorddict[atm1]),\
                                            np.array(coorddict[atm2])
                diffvector = coords_atm1 - coords_atm2
                normdiffvector = np.linalg.norm(diffvector)
                cos = np.dot(diffvector,[0,0,1])/normdiffvector
                scds_of_atoms += [0.5 * (3 * cos**2 - 1)]
            scds_of_tails += [sum(scds_of_atoms)/len(scds_of_atoms)]
            #print("TAIL:", scds_of_tails)
        totalscd = sum(scds_of_tails)/len(scds_of_tails)
        #print(totalscd)
        return totalscd

    def create_scd_histogram(self, scdfile):
        #grofile_output=self.temppath+'/calc_scd_for'+str(self.lipidmolecule)+'.gro'
        time_resid_to_scd={}
        with open(scdfile,"r") as sfile:
            sfile.readline()
            for line in sfile:
                cols = line.split()
                time = float(cols[0])
                res = int(cols[1])
                scd = float(cols[3])
                keytup = (time, res)
                time_resid_to_scd.update({keytup:scd})
        for lipid in self.systeminfo.molecules:
            if scdfile[-4:]=='.dat':
                scdfile=scdfile[:-4]
            data_output='{}_{}.dat'.format(scdfile, lipid)
            with open(data_output,"w") as outfile:
                print("{: <10} {: <10} {: <20}".format("Time","Lipid","Lipid_Scd"),file=outfile)
                endtime=int(float(time))
                for res in range(*self.systeminfo.MOLRANGE):
                    restype = self.systeminfo.resid_to_lipid[res]
                    if restype != lipid:
                        continue
                    print("Working on residue {} ".format(res),end="\r")
                    for t in range(self.systeminfo.t_start, endtime+1, self.systeminfo.dt):
                        time = float(t)
                        scd_host=float(time_resid_to_scd[(time, res)])
                        print("{: <10}{: <10}{: <20.5f}".format(time, res, scd_host),file=outfile)

def get_neighbor_of(hostres, time):
    'returns list of resids of neighbors of host resid; Inputs hostres: Resid host, time: Time as integer'
    time=float(time)
    with open("neighbor_info","r") as ninfo:
        ninfo.readline()
        for line in ninfo:
            cols = line.split()
            res = str(cols[0])
            t = cols[1]
            if res == str(hostres) and t == ''.join([str(time),'00']):
                try:
                    return cols[3].split(',')
                except IndexError:
                    return []
    print(time, "I should never get here...")

def create_leaflet_assignment_file(sysinfo_obj):
    ''' Creates a file with that assigns all lipids to upper or lower leaflet
                        !Attention!
            !Flip flops of Cholesterol are not considered! Though should it?
    '''
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
                logger.debug("Atmn/last_tail %s/%s", atomname, last_tail_atm)
                if old_resid != resid:
                    logger.debug("Resid/resname/atmname %s/%s/%s", resid, resname, atomname)
                    logger.debug("Coords head/base: %s/%s", coord_head, coord_base)
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
                if atomname in lipidmolecules.CENTRAL_ATOM_OF.values():
                    coord_head = np.array(coords)
                    logger.debug("Added head %s", atomname)
                if atomname == last_tail_atm:
                    coord_base = np.array(coords)
                    logger.debug("Added tail %s", atomname)
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
        print("UP:", sum_upper, "LOW", sum_lower)

def calc_avg_tilt(systeminfo, outputfname="bilayer_tilt.dat"):
    ''' Calculate the tilt of blocks of molecules (host + its neighbors) '''
    neiblist = gmxauto.Neighbors(systeminfo).get_neighbor_dict()
    uni = mda.Universe(systeminfo.tprpath, systeminfo.trjpath)
    len_traj = len(uni.trajectory)
    with open(outputfname, "w") as outf:
        print("{: <15}{: <5}{: <15}{: <15}{: <15}{: <9}{: <15}{: <15}{: <15}"\
            .format('time', 'res', 'scd', 'tiltangle', 'cos_tiltangle', 'leaflet', 'x' , 'y', 'z'), file=outf)
        for res in range(*systeminfo.MOLRANGE):
            for t in range(len_traj):
                time = uni.trajectory[t].time
                resname = systeminfo.resid_to_lipid[res]
                xyz = uni.select_atoms("resid {} and name {}".format(res -1, lipidmolecules.central_atom_of(resname))).positions[0]
                leaflet = systeminfo.res_to_leaflet[res]
                try:
                    neibs = neiblist[res][time]
                except KeyError:
                    continue

                #for i in neibs:
                #    resname = systeminfo.resid_to_lipid[i]
                #    if resname not in systeminfo.molecules\
                #        or lipidmolecules.is_sterol(resname)\
                #        or lipidmolecules.is_protein(resname):
                #        neibs.remove()

                lipid_block = [res] + neibs

                logger.debug("Resid %s at time %s and neibs %s", res, time, neibs)

                #if not leaflet:
                #    turn = 180
                #else:
                #    turn = 0
                cos_tiltangle = tilt_of_molecules(uni, lipid_block, time)
                tiltangle = np.arccos(cos_tiltangle)*(180/np.pi) #- turn
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

        logger.debug("resid %s\natms_res %s\nresname %s", resid, atms_res, resname)
        #if lipidmolecules.is_sterol(resname):
        #    continue


        # Choose atoms along which tilt is to be calculated
        head_atmnames = lipidmolecules.head_atoms_of(resname)
        tail_atmnames = lipidmolecules.tailcarbons_of(resname)
        logger.debug("Atomnames of tail: %s", tail_atmnames)
        tail_atmnames = [i for lt in tail_atmnames for i in lt]
        logger.debug("Atomnames of head: %s", head_atmnames)
        tail1 = tail_atmnames[int(np.round(len(tail_atmnames)/2)):]
        tail2 = tail_atmnames[:int(np.round(len(tail_atmnames)/2))]
        logger.debug("Names of tail1: %s\nand tail2: %s: ", tail1, tail2)

        # Get coords from mda_universe
        logger.debug("Selection matrix %s", get_indices(head_atmnames, atms_res))
        head_com  = selection[get_indices(head_atmnames, atms_res)].center_of_mass(pbc=True)
        tail1_com = selection[get_indices(tail1, atms_res)].center_of_mass(pbc=True)
        tail2_com = selection[get_indices(tail2, atms_res)].center_of_mass(pbc=True)
        logger.debug("COM of:\nhead:%s\ntail1:%s\ntail2:%s", head_com, tail1_com, tail2_com)

        # Calculate average tilt between vectors head -> tail1 and head -> tail2
        zvec = np.array([0, 0, 1])
        diff_ht1 = head_com-tail1_com
        diff_ht2 = head_com-tail2_com
        logger.debug("Diff1: %s, diff2: %s", diff_ht1, diff_ht2)
        cos1 = np.dot(diff_ht1, zvec) / (np.linalg.norm(diff_ht1))
        cos2 = np.dot(diff_ht2, zvec) / (np.linalg.norm(diff_ht2))
        logger.debug("cos1: %s, cos2 %s", cos1, cos2)
        avg_cos = (cos1 + cos2) / 2.0
        tilts.append(avg_cos)
    tiltangle = np.mean(tilts)
    logger.debug("Average tilt of group %s", tiltangle)
    return tiltangle

def calc_neighbor_distribution(minscd, maxscd, write='on', binwidth=0.2,
                               count_type='Ntot',
                               outfilename="p-N.dat",
                               neibfilename='neighbor_info',
                               scdfilename='scd_distribution.dat',
                               pair_type='single'):
    '''
        Calculates p(N) for lipids of type, exhibiting order around scdval+-0.1
    '''
    systeminfo = SysInfo('inputfile')
    timetoscd, maxtime = read_scdinput(scdfilename)
    resid2lipid = systeminfo.resid_to_lipid
    if pair_type == 'pair':
        neiblist = gmxauto.Neighbors(systeminfo).get_neighbor_dict()
    scdvallist = np.arange(minscd, maxscd+binwidth, binwidth)
    pofn = dict()
    typeset = set()
    with open(neibfilename, "r") as nfile:
        nfile.readline()
        for line in nfile:
            cols = line.split()
            resid = int(cols[0])
            time = float(cols[1])
            if time > maxtime:
                continue
            scdval = timetoscd[(time, resid)]
            ltype = resid2lipid[resid]
            if pair_type == 'single':
                typeset |= {ltype}
                if count_type == 'Ntot':
                    ntot = int(cols[2])
                else:
                    try:
                        neibs = cols[3]
                        ntot = int([resid2lipid[int(i)] for i in cols[3].split(',')].count(count_type))
                    except IndexError:
                        ntot = 0
                scdbin = round(scdvallist[bisect.bisect_right(scdvallist[:-1], scdval)], 2) - round(binwidth/2, 2)
                #print("COMPARE:",scdval, scdvallist, bisect.bisect_right(scdvallist[:-1], scdval), scdbin)
                #print(scdbin, scdval)
                if scdbin not in pofn.keys():
                    pofn[scdbin] = {}
                if ltype not in list(pofn[scdbin].keys()):
                    pofn[scdbin][ltype] = {}
                    #print("HERE", pofn[scdbin], ltype, pofn[scdbin][ltype])
                if ntot not in pofn[scdbin][ltype].keys():
                    pofn[scdbin][ltype][ntot] = 1
                else:
                    pofn[scdbin][ltype][ntot] += 1
            elif pair_type == 'pair':
                neibs = cols[3].split(',')
                for neib in neibs:
                    if neib <= resid:
                        continue
                    neibscd = timetoscd[(time, neib)]
                    neibtype = resid2lipid[neib]
                    avgscd = (scdval+neibscd)/2
                    pair = '{}_{}'.format(ltype, neibtype)
                    typeset |= {pair}
                    pair_neighbors = list(set(neibs) + set(neiblist[neib]) - set(neib) - set(resid))
                    if count_type == 'Ntot':
                        pair_ntot = len(pair_neighbors)
                    else:
                        pair_ntot = [resid2lipid[int(i)] for i in pair_neighbors].count(count_type)
                    scdbin = round(scdvallist[bisect.bisect_right(scdvallist, avgscd)], 2) - (binwidth/2)
                    if scdbin not in pofn.keys():
                        pofn[scdbin] = {}
                    if ltype not in pofn[scdbin].keys():
                        pofn[scdbin][pair] = {}
                    if ntot not in pofn[scdbin][pair].keys():
                        pofn[scdbin][pair][pair_ntot] = 1
                    else:
                        pofn[scdbin][pair][pair_ntot] += 1
    if write == 'on':
        write_pofn_file(outfilename, pofn, typeset, count_type=count_type)
    return pofn, typeset

def write_pofn_file(outfilename, pofn, typeset, averageing='off', count_type='Ntot'):
    with open(outfilename, "w") as outf:
        print("{: <10}{: <5}{: <20}".format('Scd', count_type, ' '.join(sorted(typeset))), file=outf)
        for scdval in sorted(pofn.keys()):
            #ltypes = pofn[scdval].keys()
            keyset = set()
            for ltype in sorted(typeset):
                try:
                    keyset |= set(pofn[scdval][ltype].keys())
                except KeyError:
                    continue
            for N in sorted(keyset):
                values = []
                for lipidtype in sorted(typeset):
                    try:
                        if averageing == 'off':
                            values += [str(pofn[scdval][lipidtype][N])]
                        else:
                            values += [str(round(np.mean(pofn[scdval][lipidtype][N])))]
                    except KeyError:
                        values += ['0']
                print("{: <10}{: <5}{: <20}".format(round(scdval, 2), N, ' '.join(values)), file=outf)

def pofn_Tavg(Membsystem, Trange, outfilename, count_type='Ntot'):
    ''' '''
    Tmin, Tmax, Tstep = Trange
    temperatures = [T for T in range(Tmin, Tmax+1, Tstep)]
    new_pofn = AutoVivification()
    for T in temperatures:
        os.chdir(Membsystem+'_'+str(T))
        pofn_T, typeset = calc_neighbor_distribution(0.0, 1.0, count_type=count_type, write='off')
        for scdval in sorted(pofn_T.keys()):
            keyset = set()
            for ltype in sorted(typeset):
                try:
                    keyset |= set(pofn_T[scdval][ltype].keys())
                except KeyError:
                    continue
            for N in sorted(keyset):
                for lipidtype in sorted(typeset):
                    try:
                        freq = pofn_T[scdval][lipidtype][N]
                    except KeyError:
                        freq = 0
                    #print(freq)
                    try:
                        new_pofn[scdval][lipidtype][N].append(freq)
                        print(new_pofn[scdval][lipidtype][N])
                    except AttributeError:
                        new_pofn[scdval][lipidtype][N] = []
                        print(new_pofn[scdval][lipidtype][N])
                        new_pofn[scdval][lipidtype][N].append(freq)
        os.chdir("../")
    write_pofn_file(outfilename, new_pofn, typeset, averageing='on')


def calc_averagedistance(self,distancefile):
    distdata={}
    neiblist=self.find_all_neighbors()
    with open(distancefile,"r") as df:
        df.readline()
        for line in df:
            col=line.split()
            time=float(col[0])
            if float(time)<float(self.t_start):
                continue
            host=int(col[1])
            print(host,time,end='\r')
            neib=int(col[2])
            if self.resid_to_lipid[host]!='DPPC' and self.resid_to_lipid[neib]!='DPPC':
                continue
            distance=float(col[3])
            try:
                neighbors=neiblist[host][time]
                nchol=[self.resid_to_lipid[int(i)] for i in neighbors].count('CHL1')
            except KeyError:
                print('Attention: Neighbor file seems not to be complete time "{}" is missing'.format(time))
            try:
                distdata[nchol]+=[distance]
            except KeyError:
                distdata.update({nchol:[distance]})
    with open("average_distance_DPPC.dat",'w') as outputf:
        for key in distdata:
            avg=sum(distdata[key])/len(distdata[key])
            print("{0: <5} {1: <20}".format(key,avg),file=outputf)

def write_neighbortype_distr(systeminfo, fname="neighbortype_distribution.dat"):
    '''
        Creates datafile < fname > with columns:
        < time >  < residue > < resname > < N comp 1 > < N comp 2 > ... < N comp i >
    '''
    neiblist = gmxauto.Neighbors(systeminfo).get_neighbor_dict()
    components = systeminfo.molecules
    with open(fname, "w") as outf:
        print("{: <12}{: <10}{: <7}".format("time", "resid", "lipidtype")\
            + (len(components)*'{: ^7}').format(*components),
            file=outf)
        for resid in range(*systeminfo.MOLRANGE):
            lipidtype = systeminfo.resid_to_lipid[resid]
            for time in range(systeminfo.t_start, systeminfo.t_end, systeminfo.dt):
                neibs = neiblist[resid][float(time)]
                neib_comp_list = []
                for lip in components:
                    ncomp = [systeminfo.resid_to_lipid[N] for N in neibs].count(lip)
                    neib_comp_list.append(ncomp)
                print("{: <12}{: <10}{: <7}".format(time, resid, lipidtype)\
                    + (len(neib_comp_list)*'{: ^7}').format(*neib_comp_list),
                    file=outf)
