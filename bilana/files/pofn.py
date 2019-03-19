import numpy as np
import bisect
from ..common import exec_gromacs

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


def pofn(systeminfo, neiblist, fname="neighbortype_distribution.dat"):
    '''
        Creates datafile < fname > with columns:
        < time >  < residue > < resname > < N comp 1 > < N comp 2 > ... < N comp i >
    '''
    components = systeminfo.molecules
    with open(fname, "w") as outf:
        print("{: <12}{: <10}{: <7}".format("time", "resid", "lipidtype")\
            + (len(components)*'{: ^7}').format(*components),
            file=outf)
        for resid in systeminfo.MOLRANGE:
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

# BROKEN
#def calc_neighbor_distribution(minscd, maxscd,
#    systeminfo,
#    write='on',
#    binwidth=0.2,
#    count_type='Ntot',
#    outfilename="p-N.dat",
#    neibfilename='neighbor_info',
#    scdfilename='scd_distribution.dat',
#    pair_type='single'):
#    '''
#        Calculates p(N) for lipids of type, exhibiting order around scdval+-0.1
#    '''
#    timetoscd, maxtime = read_scdinput(scdfilename)
#    resid2lipid = systeminfo.resid_to_lipid
#    if pair_type == 'pair':
#        neiblist = gmxauto.Neighbors(systeminfo).get_neighbor_dict()
#    scdvallist = np.arange(minscd, maxscd+binwidth, binwidth)
#    pofn = dict()
#    typeset = set()
#    with open(neibfilename, "r") as nfile:
#        nfile.readline()
#        for line in nfile:
#            cols = line.split()
#            resid = int(cols[0])
#            time = float(cols[1])
#            if time > maxtime:
#                continue
#            scdval = timetoscd[(time, resid)]
#            ltype = resid2lipid[resid]
#            if pair_type == 'single':
#                typeset |= {ltype}
#                if count_type == 'Ntot':
#                    ntot = int(cols[2])
#                else:
#                    try:
#                        neibs = cols[3]
#                        ntot = int([resid2lipid[int(i)] for i in cols[3].split(',')].count(count_type))
#                    except IndexError:
#                        ntot = 0
#                scdbin = round(scdvallist[bisect.bisect_right(scdvallist[:-1], scdval)], 2) - round(binwidth/2, 2)
#                #print("COMPARE:",scdval, scdvallist, bisect.bisect_right(scdvallist[:-1], scdval), scdbin)
#                #print(scdbin, scdval)
#                if scdbin not in pofn.keys():
#                    pofn[scdbin] = {}
#                if ltype not in list(pofn[scdbin].keys()):
#                    pofn[scdbin][ltype] = {}
#                    #print("HERE", pofn[scdbin], ltype, pofn[scdbin][ltype])
#                if ntot not in pofn[scdbin][ltype].keys():
#                    pofn[scdbin][ltype][ntot] = 1
#                else:
#                    pofn[scdbin][ltype][ntot] += 1
#            elif pair_type == 'pair':
#                neibs = cols[3].split(',')
#                for neib in neibs:
#                    if neib <= resid:
#                        continue
#                    neibscd = timetoscd[(time, neib)]
#                    neibtype = resid2lipid[neib]
#                    avgscd = (scdval+neibscd)/2
#                    pair = '{}_{}'.format(ltype, neibtype)
#                    typeset |= {pair}
#                    pair_neighbors = list(set(neibs) + set(neiblist[neib]) - set(neib) - set(resid))
#                    if count_type == 'Ntot':
#                        pair_ntot = len(pair_neighbors)
#                    else:
#                        pair_ntot = [resid2lipid[int(i)] for i in pair_neighbors].count(count_type)
#                    scdbin = round(scdvallist[bisect.bisect_right(scdvallist, avgscd)], 2) - (binwidth/2)
#                    if scdbin not in pofn.keys():
#                        pofn[scdbin] = {}
#                    if ltype not in pofn[scdbin].keys():
#                        pofn[scdbin][pair] = {}
#                    if ntot not in pofn[scdbin][pair].keys():
#                        pofn[scdbin][pair][pair_ntot] = 1
#                    else:
#                        pofn[scdbin][pair][pair_ntot] += 1
#    if write == 'on':
#        write_pofn_file(outfilename, pofn, typeset, count_type=count_type)
#    return pofn, typeset
#