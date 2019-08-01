'''
This module automates certain tasks in gromacs:
    1. Class Energy():   Calculating the interaction energies between lipids
    2. Class Neighbor(): Determine neighbors of each lipid in bilayer using cutoff and
                         reference atoms
    3. Class RDF():      Calculate radial distribution functions of sel atoms around ref atoms
    4. Class MSD():      Calculate the MSD of sel atoms
'''
import os
import numpy as np
import MDAnalysis as mda
from .. import log
from ..common import exec_gromacs, loop_to_pool, GMXNAME
from ..systeminfo import SysInfo
from ..definitions import lipidmolecules

LOGGER = log.LOGGER
LOGGER = log.create_filehandler("bilana_wrapgmx.log", LOGGER)


class Neighbors(SysInfo):
    '''
        Stores instance of class SysInfo() to get paths of all necessary MD-files
        Main task is to create file called "neighbor_info" with columns
            [<time> <resid> <N lipid neighbors> <oneighbor resids>]
        1. Use determine_neighbors() in two ways:
            A.) Create neighbor_info file using gromacs:
                1. Create a selectionfile using gmx select
                2. Use selectionfiles to create neighbor_info file (controlled in determine_neighbors())
            B.) Using 2D scheme and MDAnalysis to create neighbor_info
        2 Create index file for energy runs in create_indexfile()
    '''
    LOGGER = LOGGER

    def determine_neighbors(self, option="2D", parallel=True, **kwargs):
        ''' Controller for different determine neighbor functions '''
        LOGGER.info("Determining neighbors...")
        if option == "gmx":
            if parallel:
                self._determine_neighbors_parallel(**kwargs)
            else:
                self._determine_neighbors_serial(**kwargs)
        elif option == "2D":
            self._determine_neighbors_2d(**kwargs)
        else:
            raise ValueError("Invalid value for options.")

    def _determine_neighbors_2d(self, refatoms="name P O3", mode="atom", outputfilename="neighbor_info", _2D=True, **kwargs):
        ''' Creates neighbor_info file based on 2d distances (xy plane) '''
        traj_len = len(self.universe.trajectory)
        with open(outputfilename, "w") as outf:
            print("{: <20}{: <20}{: <20}{: <20}".format("Resid", "Time", "Number_of_neighbors", "List_of_Neighbors"), file=outf)
            for t in range(traj_len):
                time = self.universe.trajectory[t].time
                if self.t_end < time or self.t_start > time:
                    continue
                LOGGER.info("Time %s", time)
                refatomgrp = self.universe.select_atoms(refatoms)
                leaflets = self.get_ref_postions(mode, refatomgrp)
                for leaf in leaflets:

                    if _2D:
                        for i in range(len(leaf)): # Make it "2D"
                            leaf[i][1][2] = 0
                    leafoutput = {}

                    for res, coord in leaf:
                        res += 1
                        neiblist  = []
                        if res in leafoutput.keys():
                            neiblist = leafoutput[res]

                        pos_array = mda.lib.distances.distance_array(coord, np.array([j for i, j in leaf]), box=self.universe.dimensions)
                        neiblist_tmp = [leaf[i][0]+1 for i,j in enumerate(pos_array[0]) if j <= self.cutoff*10]
                        neiblist += [str(n) for n in neiblist_tmp if n != res]
                        neiblist = list(set(neiblist))
                        leafoutput[res] = neiblist

                    for res, neiblist in leafoutput.items():
                        n_neibs = len(neiblist)
                        print("{: <20}{: <20}{: <20}{: <20}".format(res, time, n_neibs, ','.join(neiblist)), file=outf)

    def get_ref_postions(self, mode, refatomgrp):
        ''' Possible modes atom, center, tails'''
        modes = ["atom", "center", "tails"]
        if mode == "atom":
            if len(refatomgrp.resids) != len(set(refatomgrp.resids)):
                raise ValueError("Refatoms string leads to more than one entry per molecule")
            center = refatomgrp.center_of_mass()
            leaf1 = [(resid, pos) for resid, pos  in enumerate(refatomgrp.positions) if pos[2] >= center[2]]
            leaf2 = [(resid, pos) for resid, pos  in enumerate(refatomgrp.positions) if pos[2] <  center[2]]
        elif mode == "center":
            center = refatomgrp.center_of_mass()
            positions = []
            for res in refatomgrp.residues:
                pos = res.atoms.center_of_mass()
                positions.append(pos)
            leaf1 = [(resid, pos) for resid, pos  in enumerate(positions) if pos[2] >= center[2]]
            leaf2 = [(resid, pos) for resid, pos  in enumerate(positions) if pos[2] <  center[2]]
        elif mode == "tails":
            center = refatomgrp.center_of_mass()
            positions = []
            reslist = []
            for atm in refatomgrp.atoms:
                resid = atm.resid - 1
                positions.append((resid, atm.position))
                reslist.append(resid)
            leaf1 = [(resid, pos) for resid, pos  in positions if pos[2] >= center[2]]
            leaf2 = [(resid, pos) for resid, pos  in positions if pos[2] <  center[2]]
        else:
            raise ValueError("Invalid mode, choose one of {}".format(modes))
        return (leaf1, leaf2)


    def _determine_neighbors_parallel(self, refatoms='P', overwrite=True, outputfilename="neighbor_info"):
        ''' Creates "neighbor_info" containing all information on lipid arrangement '''
        LOGGER.info("____Determining neighbors____\n")
        os.makedirs(self.datapath+'/neighborfiles', exist_ok=True)
        cmdlists = []
        datafileoutputs = []
        # Get command lists
        LOGGER.info("preparing files ...")
        for residue in  self.MOLRANGE:
            if residue not in self.resid_to_lipid.keys():
                continue
            selectionfile = self.create_selectionfile_neighborsearch(residue, refatoms=refatoms)
            if selectionfile is None:
                continue
            indexoutput = '{}/neighbors_of_residue{}.ndx'.format(self.indexpath, residue)
            datafileoutput = '{}/neighborfiles/neighbors_of_residue{}.dat'.format(self.datapath, residue)
            datafileoutputs.append(datafileoutput)
            if os.path.isfile(datafileoutput) and not overwrite:
                print("Neighbor file of residue {} already exists. Skipping.".format(residue))
            else:
                cmdlist=[
                    GMXNAME, 'select', '-s', self.tprpath, '-f', self.trjpath,
                    '-sf', selectionfile,'-on', indexoutput,'-oi', datafileoutput,
                    '-b', str(self.t_start), '-e', str(self.t_end),
                    '-dt', str(self.dt),
                    ]
                cmdlists.append(cmdlist)
        LOGGER.info("Running gromacs processes ...")
        if cmdlists:
            outputlogs = loop_to_pool(exec_gromacs, cmdlists)
            with open("gmx_select_determineneighbors.log","w") as logfile:
                for out, err in outputlogs:
                    logfile.write(out)
                    logfile.write(err)
        LOGGER.info("Read output")
        neibfiledata = loop_to_pool(_process_neighborfileoutput, datafileoutputs)
        neibfiledata = [i for tup in neibfiledata for i in tup]
        neibfiledata.sort()
        LOGGER.info("Print output to file %s", outputfilename)
        with open(outputfilename, "w") as outfile:
            print('{: <10}{: <25}{: <10}{: >25}'.format('Resid', 'Time', 'Number_of_neighbors', 'List_of_Neighbors'), file=outfile)
            for datatup in neibfiledata:#datatup is (residue, time, nneibs, neibindeces)
                LOGGER.debug("Writing %s", datatup)
                residue  = datatup[0]
                time     = datatup[1]
                nneibs   = datatup[2]
                neibindeces = [int(x) for x in datatup[3]]
                neibresids  = [self.index_to_resid[x] for x in neibindeces]
                neibstring  = ','.join([str(x) for x in neibresids])
                print('{: <10}{: <25}{: <10}{: >25}'.format(residue, time, nneibs, neibstring), file=outfile)


    def _determine_neighbors_serial(self, refatoms='P', overwrite=True, outputfilename="neighbor_info"):
        ''' Creates "neighbor_info" containing all information on lipid arrangement '''
        LOGGER.info("\n____Determining neighbors____\n")
        os.makedirs(self.datapath+'/neighborfiles', exist_ok=True)
        with open(outputfilename, "w") as outfile:
            print('{: <10}{: <20}{: <10}{: >25}'.format('Resid', 'Time', 'Number_of_neighbors', 'List_of_Neighbors'), file=outfile)
            for residue in  self.MOLRANGE:
                if residue not in self.resid_to_lipid.keys():
                    continue
                selectionfile = self.create_selectionfile_neighborsearch(residue, refatoms=refatoms)
                if selectionfile is None:
                    continue
                LOGGER.info(". . . Working on residue: %s . . .", residue)
                indexoutput = '{}/neighbors_of_residue{}.ndx'.format(self.indexpath, residue)
                datafileoutput = '{}/neighborfiles/neighbors_of_residue{}.dat'.format(self.datapath, residue)
                if os.path.isfile(datafileoutput) and not overwrite:
                    print("Neighbor file of residue {} already exists. Skipping.".format(residue))
                else:
                    cmdlist=[
                        GMXNAME, 'select', '-s', self.tprpath, '-f', self.trjpath,
                        '-sf', selectionfile,'-on', indexoutput,'-oi', datafileoutput,
                        '-b', str(self.t_start), '-e', str(self.t_end),
                        '-dt', str(self.dt),
                        ]
                    out, err = exec_gromacs(cmdlist)
                    with open("gmx_select.log","w") as logfile:
                        logfile.write(err)
                        logfile.write(out)
                with open(datafileoutput,"r") as datfile:
                    for line in datfile:
                        cols = line.split()
                        time = cols.pop(0)
                        nneibs = cols.pop(0)
                        neibindeces = [int(x) for x in cols]
                        neibresid = [self.index_to_resid[x] for x in neibindeces]
                        residlist = ','.join([str(x) for x in neibresid])
                        print('{: <10}{: <20}{: <10}{: >25}'.format(residue, time, nneibs, residlist), file=outfile)

    def create_selectionfile_neighborsearch(self, resid, refatoms='P'):
        ''' Create a selectionfile to get an input for gmx select
            This function is used by determine_neighbors()
            Selection depends on choice of
                refatom
                cutoff
            returns name of created selectionfile
        '''
        filename = "{}/selection_resid{}".format(self.temppath, resid)
        hosttype = self.resid_to_lipid[resid]
        if hosttype not in self.molecules:
            return
        ref_atm = lipidmolecules.central_atom_of
        # hoststring creates selection for all atoms of lipid with resid
        hoststring = "resid {0} and (name {1})"\
                     .format(resid, ref_atm(hosttype))
        # same as hoststring but for neighbors of host
        neibstring_parts = ' or '.join(
            ["(resname {} and name {})"\
             .format(i, ref_atm(i))\
             for i in self.RESNAMES],
            )
        neibstring = "(({0}) and not host) and within {1} of host"\
                     .format(neibstring_parts, self.cutoff)
        with open(filename,"w") as selection:
            if refatoms == 'P':
                print(\
                    'host =  {};\n'
                    'neibs = {};\n'
                    'neibs;'\
                    .format(hoststring, neibstring), file=selection)
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
            elif refatoms == 'glycerol':
                print(
                    'host = resid {0} and name C31 C21 O3;\n'
                    'allOAtoms = resname CHL1 and name O3 and not host;\n'
                    'allTail1Atoms = resname DPPC DUPC and name C31 and not host;\n'
                    'allTail2Atoms = resname DPPC DUPC and name C21 and not host;\n'
                    'neibOs = allOAtoms and within {1} of host;\n'
                    'neibTail1 = allTail1Atoms and within {1} of host;\n'
                    'neibTail2 = allTail2Atoms and within {1} of host;\n'
                    'neibs = neibOs or neibTail1 or neibTail2;\n'
                    'neibs;'\
                    .format(resid, self.cutoff), file=selection)
            elif refatoms == 'tails_com':
                tail_atm = lipidmolecules.tailcarbons_of
                tailstr_l = []
                for resn in self.RESNAMES:
                    for tailn in [0, 1]:
                        tailstr = "tail{}=(resname {} and name {});\n"\
                                  .format(tailn, resn, ' '.join(tail_atm(resn)[tailn]))
                        tailstr_l.append(tailstr)
                tailstr = ''.join(tailstr_l)
                print(
                    '{2}'
                    'host1 = (resid {0} and (tail0 or name O3));\n'
                    'host2 = (resid {0} and (tail1 or name O3));\n'
                    'allOAtoms = resname CHL1 and name O3 and not (host1 or host2);\n'
                    'neibOs = allOAtoms and (within {1} of com of host1 or within {1} of com of host2);\n'
                    'neibTail = (tail0 or tail1) and (within {1} of com of host1 or within {1} of com of host2);\n'
                    'neibs = neibOs or neibTail;\n'
                    'neibs;'\
                    .format(resid, self.cutoff, tailstr), file=selection)
            else:
                raise ValueError("Wrong input for refatoms with: {}"\
                                 .format(refatoms))
        return filename

    def create_indexfile(self):
        '''
            Creates an indexfile containing all indices
            of atom of each residue in system (resid_X) and all indices
            of all atoms in system.
            This index file is used for the energy calculations and stores definitions of the energygroups
        '''
        LOGGER.info("\n_____Creating index file____\n")
        resindex_all = open("resindex_all.ndx", "w")
        for mol in self.MOLRANGE:
            LOGGER.info("Working on residue %s", mol)
            selectionfilename = self.temppath+'/tmp_selectionfile'
            self.create_selectionfile_indexcreation(selectionfilename, mol)
            outputindex = self.indexpath+"/resid_"+str(mol)+".ndx"
            gmx_select_arglist = [GMXNAME, 'select', '-s',
                                  self.tprpath, '-sf',
                                  selectionfilename, '-on', outputindex,
                                  ]
            out, err = exec_gromacs(gmx_select_arglist)
            with open("gmx_select.log", "w") as logfile:
                logfile.write(err)
                logfile.write(out)
            with open(outputindex, "r") as output_index:
                filecontent = output_index.readlines()
                resindex_all.write(''.join(filecontent)+'\n\n')
        ### To have whole system indices in one group
        make_ndx_output = self.temppath+'/make_ndx_system.ndx'
        gmx_make_ndx_arglist = [GMXNAME, 'make_ndx',
                                '-f', self.gropath,
                                '-o', make_ndx_output]
        inp_str = b'keep 0\nq\n'
        out, err = exec_gromacs(gmx_make_ndx_arglist, inp_str)
        with open("gmx_make_ndx.log","w") as logfile:
            logfile.write(err)
            logfile.write(out)
        with open(make_ndx_output, "r") as output:
            filecontent = output.readlines()
            resindex_all.write(''.join(filecontent)+"\n\n")
        resindex_all.close()

    def create_selectionfile_indexcreation(self, selectionfilename, mol):
        '''
            Creates a selectionfile that is used by create_indexfile()
                Selection choose different parts of lipid chains via the atomnames
                Selections are written to file <selectionfilename>
            Example selectionfile:
                resid_204=resid 204 and resname DPPC and name ".*";         # From variable resid_defs
                resid_h_204=resid 204 and resname DPPC and name C11 C12 C13 C14 C15 N P O11 O12 O13 O14 C1 C2 O21 C21 O22 C3 O31 O32 C31 H13A H13B H13C H14A H14B H14C H15A H15B H15C H12A H12B H11A H11B HA HB HS HX HY;
                resid_204; # From lastlineitems
                resid_h_204;
        '''
        with open(selectionfilename,"w") as sf:
            lipidtype = self.resid_to_lipid[mol]
            if lipidtype not in self.molecules:
                return
            if lipidtype  not in lipidmolecules.STEROLS+lipidmolecules.PROTEINS:
                tailatoms, headatoms = self.get_atoms_tail(lipidtype)
                tailhalf12, tailhalf22 = self.get_atoms_tailhalfs(lipidtype)
                methylatomstrings = self.get_atoms_tailcarbons(lipidtype)

                selprefixes = [("", r'".*"'),\
                                ("h", ' '.join(headatoms)),\
                                ("t", ' '.join(tailatoms)),\
                                ("t12", tailhalf12),\
                                ("t22", tailhalf22),\
                                ("C0", methylatomstrings[0]),\
                                ("C1", methylatomstrings[1]),\
                                ("C2", methylatomstrings[2]),\
                                ("C3", methylatomstrings[3]),\
                                ("C4", methylatomstrings[4]),\
                                ("C5", methylatomstrings[5]),\
                                ("C6", methylatomstrings[6]),\
                                ]
                resid_defs = self.define_selection_for_resids(mol, lipidtype, selprefixes)

                lastlineitems = ['resid_{};\n'.format(str(mol))]
                lastlineitems += ['resid_{}_{};\n'.format(prefname, str(mol)) for prefname, atomnames in selprefixes if prefname != ""]
                lastlineitems = ''.join(lastlineitems)

                sf.write(resid_defs + lastlineitems)

            elif self.resid_to_lipid[mol] in lipidmolecules.STEROLS+lipidmolecules.PROTEINS:
                selectionstring = 'resid_{0}=resid {0} and resname {1};\n'\
                                    'resid_{0};'\
                                    .format(
                                        str(mol),
                                        self.resid_to_lipid[mol],
                                        )
                sf.write(selectionstring)

    def define_selection_for_resids(self, mol, lipidtype, descriptor):
        ''' This creates the lines where selections of each resid are defined '''
        selectionlist = []
        for prefname, atomnames in descriptor:
            if prefname == "":
                selectionlist += ['resid_{0}=resid {0} and resname {1} and name {2};\n'\
                                 .format(str(mol), lipidtype, atomnames)]
            else:
                selectionlist += ['resid_{0}_{1}=resid {1} and resname {2} and name {3};\n'\
                                 .format(prefname, str(mol), lipidtype, atomnames)]
        return ''.join(selectionlist)
    def get_atoms_tail(self, lipidtype):
        ''' Returns list of names of atoms in tail and head '''
        tailatomlist = lipidmolecules.tail_atoms_of(lipidtype)
        tailatoms = [x for chain in tailatomlist for x in chain] ##unpacking
        headatoms = lipidmolecules.head_atoms_of(lipidtype)
        return tailatoms, headatoms
    def get_atoms_tailhalfs(self, lipidtype):
        ''' Returns a list of atoms for the halfs of each lipid chain '''
        tailhalf12_l = [] ### to get half the tails
        tailhalf22_l = []
        for molpart in lipidmolecules.tail_atoms_of(lipidtype):
            tailhalf12_l.extend(molpart[:len(molpart)//2])
            tailhalf22_l.extend(molpart[len(molpart)//2:])
        tailhalf12 = ' '.join(tailhalf12_l)
        tailhalf22 = ' '.join(tailhalf22_l)
        return tailhalf12, tailhalf22
    def get_atoms_tailcarbons(self, lipidtype):
        ''' Returns a list of carbon atom names '''
        methylatomslists = []
        for tailindex in range(len(lipidmolecules.tailcarbons_of(lipidtype))):
            handlinghydr = [iter(lipidmolecules.tailhydr_of(lipidtype)[tailindex])]*2
            methylgrouplist = [i for i in zip(lipidmolecules.tailcarbons_of(lipidtype)[tailindex], zip(*handlinghydr))]# Getting a list of tuples like [C1,(H1,H2),....]
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
        methylatomstrings = [' '.join(t) for i in methylgroups for t in i]
        return methylatomstrings


def get_neighbor_dict(neighbor_filename='neighbor_info', verbose='off'):
    ''' Returns a list of all neighbors being in the
        cutoff distance at least once in the trajectory.
        Neighborfile is required and is output of determine_neighbors()

        Dict layout is:
        neibdict[resid][time] --> [neibs]
    '''
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
                    print("No neighbors of residue {} at time {}."\
                            .format(cols[0], cols[1]))
    return neibdict

def _process_neighborfileoutput(datafileoutput):
    ''' Definition on module level to be picklable for parallelization '''
    residue = int(os.path.basename(datafileoutput.replace("neighbors_of_residue", "")).replace(".dat", ""))
    data = []
    with open(datafileoutput,"r") as datfile:
        for line in datfile:
            cols = line.split()
            time = cols.pop(0)
            nneibs = cols.pop(0)
            neibindeces = cols[:]
            #neibindeces = [int(x) for x in cols]
            #neibresid = [in2res[x] for x in neibindeces]
            #neibindexlist = ','.join([str(x) for x in neibresid])
            data.append((residue, time, nneibs, neibindeces))
    return data
