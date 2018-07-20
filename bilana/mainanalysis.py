''' Set of different tools for analysis
    Specify more here:
        -
 '''
import os
import sys 
import numpy as np
import pprint
import re
from time import localtime, strftime
from bilana.common import AutoVivification, GRO_format
from bilana.energyfilecreator import read_scdinput
from bilana import lipidmolecules
from bilana import gromacstoolautomator as gmxauto
from bilana.systeminfo import SysInfo
import bisect

pp = pprint.PrettyPrinter()

def is_neighbor_in_leaflet(systeminfo_inst):
    ''' Searches for interleaflet neighborhood '''
    leaflet_assign = systeminfo_inst.res_to_leaflet
    nlipids = systeminfo_inst.NUMBEROFMOLECULES
    neiblist = gmxauto.Neighbors(systeminfo_inst).get_neighbor_dict()
    host_has_interleafletneib = []
    for host in range(1, nlipids+1):
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
        #print(self.atomlist)
    def create_scdfile(self, grofilename=None, outputfile='scd_distribution.dat'):
        if grofilename is None:
            grofilename = 'traj_complete.gro'
            grofilepath = ''.join([self.systeminfo.datapath, '/grofiles/'])
            gmxauto.produce_gro(self.systeminfo, grofilename=grofilename)
        with open(grofilepath+grofilename,"r") as grofile, open(outputfile, "w") as scdfile:
            print("Time     Residue     Type      Scd", file=scdfile)
            resid_old = 1
            time = 0
            lipidtype_old = ''
            coorddict = {}
            print(strftime("%H:%M:%S :", localtime()),"... read data from .gro-file ...")
            #========================RESID=======RESNAME======ATOMNAME=======INDEX===============X============Y============Z=======
            #regexp = re.compile(r'^([\s,\d]{5})([\w,\s]{5})([\d,\w,\s]{5})([\s*\d+]{5})(\s*-?\d+\.\d+\s*-?\d+\.\d+\s*-?\d+\.\d+).*')
            for line in grofile:
                if 't=' in line:
                    time_tmp = float(line[line.index('t=')+2:].strip())   #to get rid of obsolete decimals
                    if not time:
                        time = time_tmp
                    #print("...at time {}".format(time), flush=True, end="\r")
                    if float(self.systeminfo.t_end) < time:
                        #print("breaking at", time)
                        break
                #print("Match", regexp.match(line), line[:15], end='\n')
                #print("Time match is", float(systeminfo.t_start)<=time, end='\n')
                regmatch = GRO_format.regexp.match(line)
                if float(self.systeminfo.t_start) <= time and regmatch is not None:
                    grps = regmatch.groups()
                    #print(grps)
                    #print("Reading data at", time, end='\n')
                    #atom = line[9:15].strip()
                    #lipidtype = line[5:9]
                    atom = grps[2].strip()
                    lipidtype = grps[1].strip()
                    if not lipidtype_old:
                        lipidtype_old = lipidtype
                    all_atmlst = [atm for atmlst in self.atomlist[lipidtype] for atm in atmlst]
                    if atom in all_atmlst and lipidtype in self.systeminfo.molecules:
                        #resid = int(line[:5].strip())
                        resid = int(grps[0].strip())
                        if resid != resid_old and coorddict:
                            if not resid_old:
                                resid_old = resid
                            scd_value = self.scd_of_res(coorddict, self.atomlist[lipidtype_old])
                            coorddict = {}
                            print("{: <10}{: <8}{: <6}{: <20}".format(time, resid_old, lipidtype_old, scd_value), file=scdfile)
                            resid_old = resid
                            time = time_tmp
                            lipidtype_old = lipidtype
                        #coordinates = [float(x) for x in line[20:44].split()]
                        coordinates = [float(x) for x in grps[4:7]]
                        coorddict[atom] = coordinates
                    else:
                        continue
            scd_value = self.scd_of_res(coorddict, self.atomlist[lipidtype_old])
            coorddict = {}
            print("{: <10}{: <8}{: <6}{: <20}".format(time, resid_old, lipidtype_old, scd_value), file=scdfile)
        print(strftime("%H:%M:%S :", localtime()),"Finished reading.")
        return 
    def scd_of_res(self, coorddict, atomlist):#, neibstraightness=0,neiblist=None):
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
                for res in range(1, self.systeminfo.NUMBEROFMOLECULES+1):
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
            cols=line.split()
            res=str(cols[0])
            t=cols[1]   
            if res==str(hostres) and t==''.join([str(time),'00']):
                try:
                    return cols[3].split(',')
                except IndexError:
                    return []
    print(time,"I should never get here...")

def create_leaflet_assignment_file(sysinfo_obj):
    ''' Creates a file with that assigns all lipids to upper or lower leaflet
                        !Attention! 
            !Flip flops of Cholesterol are not considered! Though should it?
    '''
    outputdict = {}
    outputfilename = 'leaflet_assignment.dat'
    grofile_path = sysinfo_obj.gropath
    with open(grofile_path, "r") as gfile:
                #========================RESID=======RESNAME======ATOMNAME=======INDEX===============X============Y============Z=======
        regexp = re.compile(r'^([\s,\d]{5})([\w,\s]{5})([\d,\w,\s]{5})([\s*\d+]{5})\s*(-?\d+\.\d+\s*-?\d+\.\d+\s*-?\d+\.\d+).*')
        old_resid = 1
        sum_upper = 0
        sum_lower = 0
        for line in gfile:
            match = regexp.match(line)
            if match:
                resid = int(match.group(1).split()[0])
                atomname = match.group(3).split()[0]
                coords = [float(i) for i in match.group(5).split()]
                #print(resid, atomname, coords)
                if atomname in lipidmolecules.central_atom_of.values():
                    coord_head = np.array(coords)
                if atomname in [i[-1] for it in lipidmolecules.scd_tail_atoms_of.values() for i in it]:
                    coord_base = np.array(coords)
                if old_resid != resid:
                    if coord_head is None or coord_base is None:
                        continue
                    new_coords = coord_head - coord_base
                    cos = np.dot(new_coords, np.array([0.0,0.0,1.0]))/np.linalg.norm(new_coords)
                    if cos <= 0:
                        sum_upper += 1
                        outputdict[old_resid] = 0
                    else:
                        sum_lower += 1
                        outputdict[old_resid] = 1
                    #vecdir = np.arccos(cos)
                    #print(vecdir, coord_head, coord_base)
                    old_resid = resid
                    coord_head = coord_base = None
        print("UP:", sum_upper, "LOW", sum_lower)
    with open(outputfilename, "w") as outf:
        print("{: <7} {: <5}".format('resid', 'leaflet'), file=outf)
        for res in range(1, sysinfo_obj.NUMBEROFMOLECULES+1):
            print("{: <7} {: <5}".format(res, outputdict[res]), file=outf)

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

