'''
This module automates certain tasks in gromacs:
    1. Class Energy():   Calculating the interaction energies between lipids
    2. Class Neighbor(): Determine neighbors of each lipid in bilayer using cutoff and
                         reference atoms
    3. Class RDF():      Calculate radial distribution functions of sel atoms around ref atoms
    4. Class MSD():      Calculate the MSD of sel atoms
'''
import os, sys
import numpy as np
import pandas as pd
import MDAnalysis as mda
from .. import log
from . import protein
from ..common import exec_gromacs, loop_to_pool, GMXNAME
from ..systeminfo import SysInfo
from ..definitions import lipidmolecules
from .leaflets import molecule_leaflet_orientation

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
            self._determine_neighbors_2d(parallel=parallel, **kwargs)
        else:
            raise ValueError("Invalid value for options.")

    def _determine_neighbors_2d(self, mode="atom", outputfilename="neighbor_info", _2D=True, parallel=True, **kwargs):
        ''' Creates neighbor_info file based on 2d distances (xy plane) '''
        DEBUG = False
        if DEBUG:
            LOGGER.setLevel("DEBUG")

        refatoms = self.reference_atom_selection

        traj_len = len(self.universe.trajectory)
        inpargs = []
        LOGGER.info("Collect inpargs")
        outp = []
        leaflets = [[[],], [[],]]
        for t in range(traj_len):
            time = self.universe.trajectory[t].time
            if self.t_end < time or self.t_start > time:
                continue
            if time % self.dt != 0:
                continue

            refatomgrp = self.universe.select_atoms(refatoms)
            LOGGER.debug("Found %s atoms: %s ", len(refatomgrp.atoms), refatomgrp.atoms)
            leaflets_tmp = leaflets
            leaflets = self.get_ref_positions(mode, refatomgrp) # leaflets=[(resid1, pos1), ...]
            LOGGER.debug("Dimension of leaflets %s", np.array(leaflets).shape)
            if len(leaflets[0]) != len(leaflets[1]) and np.array(leaflets_tmp).shape[-1] != 0:
                #raise ValueError("Number of found ref positions differs: {} vs {}".format(len(leaflets[0]), len(leaflets[1])))
                LOGGER.warning("frame %s: Number of found ref positions differs: %s vs %s", t, len(leaflets[0]), len(leaflets[1]))
                LOGGER.warning("Difference1: %s", set([i[0] for i in leaflets[0]]).symmetric_difference(set([i[0] for i in leaflets_tmp[0]])))

            # universe.dimensions has to be copied!! otherwise reference will be changed and always last boxdimensions are used
            inp = (time, leaflets, self.universe.dimensions.copy(), self.cutoff, _2D)
            if parallel:
                inpargs.append(inp)
            else:
                outp.append(self._get_outp_line(*inp))

        if parallel:
            LOGGER.info("Sending jobs to pool")
            outp = loop_to_pool(self._get_outp_line, inpargs, maxtasknum=8)
            outp = [i for l in outp for i in l]
        else:
            outp = [i for l in outp for i in l]

        outp.sort()
        LOGGER.info("Write output file")
        with open(outputfilename, "w") as outf:
            print("{: <20}{: <20}{: <20}{: <20}".format("Resid", "Time", "Number_of_neighbors", "List_of_Neighbors"), file=outf)
            for line in outp:
                print("{: <20}{: <20}{: <20}{: <20}".format(*line), file=outf)

    @staticmethod
    def _get_outp_line(time, leaflets, boxdim, cutoff, twodimensional):
        ''' Looks for all neighbors of hostid with host_pos within cutoff '''
        outp = []
        for all_coords_per_leaflet in leaflets:

            if twodimensional:
                for i in range(len(all_coords_per_leaflet)): # Make it "2D"
                    all_coords_per_leaflet[i][1][2] = 0      # By setting z values to 0

            for hostid, host_pos in all_coords_per_leaflet:
                position_array =  np.array([pos for resid, pos in all_coords_per_leaflet]) # Get all positions in leaflet selection
                dist_array = mda.lib.distances.distance_array(host_pos, position_array, box=boxdim)[0] # output is [ [[dist1], [dist2], ...] ]

                neiblist = []
                for resndx, distance in enumerate(dist_array):


                    if distance <= cutoff*10.0:
                        LOGGER.debug("host %s, resndx %s, neib %s", hostid, resndx, all_coords_per_leaflet[resndx][0])
                        LOGGER.debug("hostpos vs neibpos: %s vs %s", host_pos, all_coords_per_leaflet[resndx][1] )
                        LOGGER.debug("box %s", boxdim)
                        LOGGER.debug("dist %s and dist_self %s", distance, dist_helper(host_pos,  all_coords_per_leaflet[resndx][1], boxdim))
                        neiblist.append(all_coords_per_leaflet[resndx][0])

                    #else:
                        #LOGGER.debug("REJECTED: host%s  neib%s", hostid, all_coords_per_leaflet[resndx][0])

                neiblist = list(set(neiblist)) # delete duplicates
                neiblist.sort()
                neiblist = [ str(n) for n in neiblist if n != hostid ] # delete host entry
                n_neibs = len(neiblist)
                outp.append( (hostid, time, n_neibs, ','.join(str(i) for i in neiblist)) )
        LOGGER.info("finished at %s", time)
        return outp

    def get_ref_positions(self, mode, refatomgrp):
        ''' Possible modes atom, center, tails'''
        modes = ["atom", "center", "tails"]
        if mode == "atom":
            if len(refatomgrp.resids) != len(set(refatomgrp.resids)):
                raise ValueError("Refatoms string leads to more than one entry per molecule: {}-{}".format(refatomgrp, refatomgrp.resids))
            orientations = {}
            for residue in refatomgrp.residues:
                headsel = 'resname {} and name {}'.format(residue.resname, ' '.join( lipidmolecules.head_atoms_of(residue.resname) ) )
                tailsel = 'resname {} and name {}'.format(residue.resname, ' '.join( [taillist[-1] for taillist in lipidmolecules.tailcarbons_of(residue.resname) ] ) )
                head_pos = residue.atoms.select_atoms( headsel ).center_of_geometry()
                tail_pos = residue.atoms.select_atoms( tailsel ).center_of_geometry()
                orientations[residue.resid] =  molecule_leaflet_orientation( head_pos, tail_pos )
            #center = refatomgrp.center_of_mass()
            #leaf1 = [(atm.resid, atm.position) for atm  in refatomgrp.atoms if atm.position[2] >= center[2]]
            #leaf2 = [(atm.resid, atm.position) for atm  in refatomgrp.atoms if atm.position[2] <  center[2]]
            leaf1 = np.array( [(atm.resid, atm.position) for atm  in refatomgrp.atoms if orientations[ atm.resid ] ] )
            leaf2 = np.array( [(atm.resid, atm.position) for atm  in refatomgrp.atoms if not orientations[ atm.resid ] ] )
            LOGGER.debug("dim leaf1 %s and leaf2 %s ", leaf1.shape, leaf2.shape)

        elif mode == "tails":
            cnt = 0
            orientations = {}
            for residue in refatomgrp.residues:
                headsel = 'resname {} and name {}'.format(residue.resname, ' '.join( lipidmolecules.head_atoms_of(residue.resname) ) )
                tailsel = 'resname {} and name {}'.format(residue.resname, ' '.join( np.array( lipidmolecules.tailcarbons_of(residue.resname) )[:,-1] ) )
                LOGGER.debug("Head selection: %s", headsel)
                LOGGER.debug("Tail selection: %s", tailsel)
                head_pos = residue.atoms.select_atoms( headsel ).center_of_mass()
                tail_pos = residue.atoms.select_atoms( tailsel ).center_of_mass()
                orientations[residue.resid] =  molecule_leaflet_orientation( head_pos, tail_pos )
            leaf1 = np.array( [(atm.resid, atm.position) for atm  in refatomgrp.atoms if orientations[ atm.resid ] ] )
            leaf2 = np.array( [(atm.resid, atm.position) for atm  in refatomgrp.atoms if not orientations[ atm.resid ] ] )
            LOGGER.debug("Leaf1 is\n%s", leaf1)
            for resid  in refatomgrp.resids:
                resid += cnt
                mask1 = np.where( leaf1[:,0] == resid )[0] # Index of duplicate per resid
                mask2 = np.where( leaf2[:,0] == resid )[0] # are stored here
                if mask1.shape[0] == 2:
                    LOGGER.debug("Index of second tail: %s", np.where( leaf1[:,0] == resid ))
                    cnt += 1 # count number of additional resids
                    leaf1[ leaf1[:,0] > resid ] += 1 # increase resid number above processed resid
                    leaf2[ leaf2[:,0] > resid ] += 1 # for each duplicate
                    LOGGER.debug("leaf1 now %s", leaf1)
                    leaf1[ mask1[1] ][0] += 1  # Increase resid of second item in resid duplicate
                    LOGGER.debug("leaf1 and then %s", leaf1)
                elif mask2.shape[0] == 2: # Same here as for mask1
                    LOGGER.debug("Index of second tail: %s", np.where( leaf2[:,0] == resid ))
                    cnt += 1
                    leaf1[ leaf1[:,0] > resid ] += 1
                    leaf2[ leaf2[:,0] > resid ] += 1
                    LOGGER.debug("leaf2 now %s", leaf2)
                    leaf2[ mask2[1] ][0] += 1
                    LOGGER.debug("leaf2 and then %s", leaf2)
            leaf1, leaf2 = list(leaf1), list(leaf2)


        else:
            raise ValueError("Invalid mode, choose one of {}".format(modes))
        return (leaf1, leaf2)

    def determine_neighbors_protein(self, mode="atom", prot_cutoff=1.2, outputfilename="neighbor_info_protein", overwrite=False, **kwargs):
        '''
            Determine all neighbors of a protein, based on cutoff distance
            neighbor_info_prot then looks like:
            < time > < resid >
        '''
        DEBUG = False
        if DEBUG:
            LOGGER.setLevel("DEBUG")

        refatoms = self.reference_atom_selection

        traj_len = len(self.universe.trajectory)
        inpargs = []
        LOGGER.info("Collect inpargs")
        outp = []
        leaflets = [[[],], [[],]]
        if not overwrite and os.path.isfile(outputfilename):
            LOGGER.info("File {} exists, will not overwrite".format(outputfilename))
        with open(outputfilename, "w") as outf:
            print("{: <20}{: <20}{: <20}{: <20}".format("Resid", "Time", "Number_of_neighbors", "List_of_Neighbors"), file=outf)
            for t in range(traj_len):
                time = self.universe.trajectory[t].time
                if self.t_end < time or self.t_start > time:
                    continue

                refatomgrp = self.universe.select_atoms(refatoms)
                LOGGER.debug("Found %s atoms: %s ", len(refatomgrp.atoms), refatomgrp.atoms)
                leaflets_tmp = leaflets
                leaflets = self.get_ref_positions(mode, refatomgrp) # leaflets=( [(resid1, pos1), ...], [(residN, posN), ...] )

                # Get correct protein position
                ref_res_prot = protein.get_reference_resids() # [pos_leaf1, pos,leaf2]
                sel1 = "resid {}".format(' '.join([str(i) for i in ref_res_prot[0]]))
                sel2 = "resid {}".format(' '.join([str(i) for i in ref_res_prot[1]]))
                prot_pos_leaf1 = self.universe.select_atoms( sel1 ).center_of_mass()
                prot_pos_leaf2 = self.universe.select_atoms( sel2 ).center_of_mass()
                leaflets = [ [(0, prot_pos_leaf1), *leaflets[0]], [(0, prot_pos_leaf2), *leaflets[1]] ]

                if len(leaflets[0]) != len(leaflets[1]):
                    #raise ValueError("Number of found ref positions differs: {} vs {}".format(len(leaflets[0]), len(leaflets[1])))
                    LOGGER.warning("frame %s: Number of found ref positions differs: %s vs %s", t, len(leaflets[0]), len(leaflets[1]))
                    LOGGER.warning("Difference1: %s", set([i[0] for i in leaflets[0]]).symmetric_difference(set([i[0] for i in leaflets_tmp[0]])))

                boxdim = self.universe.dimensions.copy()
                cutoff = self.cutoff

                for all_coords_per_leaflet in leaflets:

                    for i in range(len(all_coords_per_leaflet)): # Make it "2D"
                        all_coords_per_leaflet[i][1][2] = 0      # By setting z values to 0

                    prot_res, prot_pos = all_coords_per_leaflet[0]
                    position_array =  np.array([pos for resid, pos in all_coords_per_leaflet]) # Get all positions in leaflet selection
                    dist_array = mda.lib.distances.distance_array(prot_pos, position_array, box=boxdim)[0] # output is [ [[dist1], [dist2], ...] ]

                    neiblist = []
                    for resndx, distance in enumerate(dist_array):
                        if distance <= prot_cutoff*10.0:
                            LOGGER.debug("prot %s, resndx %s, neib %s", prot_res, resndx, all_coords_per_leaflet[resndx][0])
                            LOGGER.debug("protpos vs neibpos: %s vs %s", prot_pos, all_coords_per_leaflet[resndx][1] )
                            LOGGER.debug("box %s", boxdim)
                            LOGGER.debug("dist %s and dist_self %s", distance, dist_helper(prot_pos,  all_coords_per_leaflet[resndx][1], boxdim))
                            neiblist.append(all_coords_per_leaflet[resndx][0])

                    neiblist = list(set(neiblist)) # delete duplicates
                    neiblist.sort()
                    neiblist = [ str(n) for n in neiblist if n != prot_res] # delete prot entry
                    n_neibs = len(neiblist)
                    line =  (prot_res, time, n_neibs, ','.join(str(i) for i in neiblist))

                    print("{: <20}{: <20}{: <20}{: <20}".format(*line), file=outf)


    #def _determine_neighbors_parallel(self, refatoms='P', overwrite=True, outputfilename="neighbor_info"):
    #    ''' Creates "neighbor_info" containing all information on lipid arrangement '''
    #    LOGGER.info("____Determining neighbors____\n")
    #    os.makedirs(self.datapath+'/neighborfiles', exist_ok=True)
    #    cmdlists = []
    #    datafileoutputs = []
    #    # Get command lists
    #    LOGGER.info("preparing files ...")
    #    LOGGER.debug("for mols: %s", self.MOLRANGE)
    #    for residue in  self.MOLRANGE:
    #        if residue not in self.resid_to_lipid.keys():
    #            continue
    #        selectionfile = self.create_selectionfile_neighborsearch(residue, refatoms=refatoms)
    #        if selectionfile is None:
    #            continue
    #        indexoutput = '{}/neighbors_of_residue{}.ndx'.format(self.indexpath, residue)
    #        datafileoutput = '{}/neighborfiles/neighbors_of_residue{}.dat'.format(self.datapath, residue)
    #        datafileoutputs.append(datafileoutput)
    #        if os.path.isfile(datafileoutput) and not overwrite:
    #            print("Neighbor file of residue {} already exists. Skipping.".format(residue))
    #            LOGGER.debug("Skipping %s", residue)
    #        else:
    #            cmdlist=[
    #                GMXNAME, 'select', '-s', self.tprpath, '-f', self.trjpath,
    #                '-sf', selectionfile,'-on', indexoutput,'-oi', datafileoutput,
    #                '-b', str(self.t_start), '-e', str(self.t_end),
    #                '-dt', str(self.dt),
    #                ]
    #            cmdlists.append(cmdlist)
    #    LOGGER.info("Running gromacs processes ...")
    #    LOGGER.debug("CMDS %s", cmdlists)
    #    if cmdlists:
    #        outputlogs = loop_to_pool(exec_gromacs, cmdlists)
    #        with open("gmx_select_determineneighbors.log","w") as logfile:
    #            for out, err in outputlogs:
    #                logfile.write(out)
    #                logfile.write(err)
    #    LOGGER.info("Read output")
    #    neibfiledata = loop_to_pool(_process_neighborfileoutput, datafileoutputs)
    #    neibfiledata = [i for tup in neibfiledata for i in tup]
    #    neibfiledata.sort()
    #    LOGGER.info("Print output to file %s", outputfilename)
    #    with open(outputfilename, "w") as outfile:
    #        print('{: <10}{: <25}{: <10}{: >25}'.format('Resid', 'Time', 'Number_of_neighbors', 'List_of_Neighbors'), file=outfile)
    #        for datatup in neibfiledata:#datatup is (residue, time, nneibs, neibindeces)
    #            LOGGER.debug("Writing %s", datatup)
    #            residue  = datatup[0]
    #            time     = datatup[1]
    #            nneibs   = datatup[2]
    #            neibindeces = [int(x) for x in datatup[3]]
    #            neibresids  = [self.index_to_resid[x] for x in neibindeces]
    #            neibstring  = ','.join([str(x) for x in neibresids])
    #            print('{: <10}{: <25}{: <10}{: >25}'.format(residue, time, nneibs, neibstring), file=outfile)


    #def _determine_neighbors_serial(self, refatoms='P', overwrite=True, outputfilename="neighbor_info"):
    #    ''' Creates "neighbor_info" containing all information on lipid arrangement '''
    #    LOGGER.info("\n____Determining neighbors____\n")
    #    os.makedirs(self.datapath+'/neighborfiles', exist_ok=True)
    #    with open(outputfilename, "w") as outfile:
    #        print('{: <10}{: <20}{: <10}{: >25}'.format('Resid', 'Time', 'Number_of_neighbors', 'List_of_Neighbors'), file=outfile)
    #        for residue in  self.MOLRANGE:
    #            if residue not in self.resid_to_lipid.keys():
    #                continue
    #            selectionfile = self.create_selectionfile_neighborsearch(residue, refatoms=refatoms)
    #            if selectionfile is None:
    #                continue
    #            LOGGER.info(". . . Working on residue: %s . . .", residue)
    #            indexoutput = '{}/neighbors_of_residue{}.ndx'.format(self.indexpath, residue)
    #            datafileoutput = '{}/neighborfiles/neighbors_of_residue{}.dat'.format(self.datapath, residue)
    #            if os.path.isfile(datafileoutput) and not overwrite:
    #                print("Neighbor file of residue {} already exists. Skipping.".format(residue))
    #            else:
    #                cmdlist=[
    #                    GMXNAME, 'select', '-s', self.tprpath, '-f', self.trjpath,
    #                    '-sf', selectionfile,'-on', indexoutput,'-oi', datafileoutput,
    #                    '-b', str(self.t_start), '-e', str(self.t_end),
    #                    '-dt', str(self.dt),
    #                    ]
    #                out, err = exec_gromacs(cmdlist)
    #                with open("gmx_select.log","w") as logfile:
    #                    logfile.write(err)
    #                    logfile.write(out)
    #            with open(datafileoutput,"r") as datfile:
    #                for line in datfile:
    #                    cols = line.split()
    #                    time = cols.pop(0)
    #                    nneibs = cols.pop(0)
    #                    neibindeces = [int(x) for x in cols]
    #                    neibresid = [self.index_to_resid[x] for x in neibindeces]
    #                    residlist = ','.join([str(x) for x in neibresid])
    #                    print('{: <10}{: <20}{: <10}{: >25}'.format(residue, time, nneibs, residlist), file=outfile)

    #def create_selectionfile_neighborsearch(self, resid, refatoms='P'):
    #    ''' Create a selectionfile to get an input for gmx select
    #        This function is used by determine_neighbors()
    #        Selection depends on choice of
    #            refatom
    #            cutoff
    #        returns name of created selectionfile
    #    '''
    #    filename = "{}/selection_resid{}".format(self.temppath, resid)
    #    hosttype = self.resid_to_lipid[resid]
    #    if hosttype not in self.molecules:
    #        return
    #    ref_atm = lipidmolecules.central_atom_of
    #    # hoststring creates selection for all atoms of lipid with resid
    #    hoststring = "resid {0} and (name {1})"\
    #                 .format(resid, ref_atm(hosttype))
    #    # same as hoststring but for neighbors of host
    #    neibstring_parts = ' or '.join(
    #        ["(resname {} and name {})"\
    #         .format(i, ref_atm(i))\
    #         for i in self.RESNAMES],
    #        )
    #    neibstring = "(({0}) and not host) and within {1} of host"\
    #                 .format(neibstring_parts, self.cutoff)
    #    with open(filename,"w") as selection:
    #        if refatoms == 'P':
    #            print(\
    #                'host =  {};\n'
    #                'neibs = {};\n'
    #                'neibs;'\
    #                .format(hoststring, neibstring), file=selection)
    #        elif refatoms == 'bothtails':
    #            print(
    #                'host = resid {0} and name C34 C24 O3;\n'
    #                'allOAtoms = resname CHL1 and name O3 and not host;\n'
    #                'allTail1Atoms = resname DPPC DUPC and name C34 and not host;\n'
    #                'allTail2Atoms = resname DPPC DUPC and name C24 and not host;\n'
    #                'neibOs = allOAtoms and within {1} of host;\n'
    #                'neibTail1 = allTail1Atoms and within {1} of host;\n'
    #                'neibTail2 = allTail2Atoms and within {1} of host;\n'
    #                'neibs = neibOs or neibTail1 or neibTail2;\n'
    #                'neibs;'\
    #                .format(resid, self.cutoff), file=selection)
    #        elif refatoms == 'glycerol':
    #            print(
    #                'host = resid {0} and name C31 C21 O3;\n'
    #                'allOAtoms = resname CHL1 and name O3 and not host;\n'
    #                'allTail1Atoms = resname DPPC DUPC and name C31 and not host;\n'
    #                'allTail2Atoms = resname DPPC DUPC and name C21 and not host;\n'
    #                'neibOs = allOAtoms and within {1} of host;\n'
    #                'neibTail1 = allTail1Atoms and within {1} of host;\n'
    #                'neibTail2 = allTail2Atoms and within {1} of host;\n'
    #                'neibs = neibOs or neibTail1 or neibTail2;\n'
    #                'neibs;'\
    #                .format(resid, self.cutoff), file=selection)
    #        elif refatoms == 'tails_com':
    #            tail_atm = lipidmolecules.tailcarbons_of
    #            tailstr_l = []
    #            for resn in self.RESNAMES:
    #                for tailn in [0, 1]:
    #                    tailstr = "tail{}=(resname {} and name {});\n"\
    #                              .format(tailn, resn, ' '.join(tail_atm(resn)[tailn]))
    #                    tailstr_l.append(tailstr)
    #            tailstr = ''.join(tailstr_l)
    #            print(
    #                '{2}'
    #                'host1 = (resid {0} and (tail0 or name O3));\n'
    #                'host2 = (resid {0} and (tail1 or name O3));\n'
    #                'allOAtoms = resname CHL1 and name O3 and not (host1 or host2);\n'
    #                'neibOs = allOAtoms and (within {1} of com of host1 or within {1} of com of host2);\n'
    #                'neibTail = (tail0 or tail1) and (within {1} of com of host1 or within {1} of com of host2);\n'
    #                'neibs = neibOs or neibTail;\n'
    #                'neibs;'\
    #                .format(resid, self.cutoff, tailstr), file=selection)
    #        else:
    #            raise ValueError("Wrong input for refatoms with: {}"\
    #                             .format(refatoms))
    #    return filename

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
        self.add_water_groups_to_index(self.SOLVENT)
        #self.add_leaflet_groups_to_index()

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
            #handlinghydr = [iter(lipidmolecules.tailhydr_of(lipidtype)[tailindex])]*2
            #methylgrouplist = [i for i in zip(lipidmolecules.tailcarbons_of(lipidtype)[tailindex], zip(*handlinghydr))]# Getting a list of tuples like [C1,(H1,H2),....]
            methylgrouplist = [i for i in zip(lipidmolecules.tailcarbons_of(lipidtype)[tailindex])]# Getting a list of tuples like [C1,(H1,H2),....]
            print(methylgrouplist)
            methylgrouplist = [i for tup in methylgrouplist for i in tup]# Unpacking this list
            print(methylgrouplist)
            methylgrouplist_unp = []
            for particle in methylgrouplist:# get rid of hydr tuples
                if isinstance(particle, tuple):
                    methylgrouplist_unp.extend(particle)
                else:
                    methylgrouplist_unp.append(particle)
            methylatomslists.append(methylgrouplist_unp)
            print(methylatomslists)
        methylgroups = [[methylatomslists[0][i:i+3]\
                        +methylatomslists[0][i+3:i+6]\
                        +methylatomslists[1][i:i+3]\
                        +methylatomslists[1][i+3:i+6]]\
                   for i in range(0, len(methylatomslists[0])-1, 6)]
        methylatomstrings = [' '.join(t) for i in methylgroups for t in i]
        print(methylatomstrings)
        return methylatomstrings

    def add_water_groups_to_index(self, watername, add_grp_to="resindex_all.ndx"):
        ''' Make index file using gmx select and append index group to <add_grp_to> '''
        selectionstr = 'solv=resname {}; solv;'.format(watername)
        outputsel = self.temppath + "/tmp_leaflet.ndx"
        cmd = [
            "gmx", "select", "-f", self.gropath,
            "-s", self.tprpath, "-select", selectionstr,
            "-on", outputsel
        ]
        out, err = exec_gromacs(cmd)
        print(out, err)
        with open(outputsel, "r") as selectionf, open(add_grp_to, "a") as ndxf:
            for line in selectionf:
                ndxf.write(line)

    def add_leaflet_groups_to_index(self, leaflet_assignment_fname, add_grp_to="resindex_all.ndx"):
        ''' Make index file using gmx select and append index group to <add_grp_to> '''
        leafdat = pd.read_table(leaflet_assignment_fname, delim_whitespace=True)
        resid_list = [[], []]
        outputsel = self.temppath + "/tmp_leaflet.ndx"
        for leafndx, fr in leafdat.groupby("leaflet"):
            resid_list[leafndx] += list(fr.resid)
        selectionstr = 'leaflet0=resname {0} and resid {1}; leaflet1=resname {0} and resid {2}; leaflet0; leaflet1;'.format(
                ' '.join(self.molecules),
                ' '.join( [ str(i) for i in resid_list[0] ] ),
                ' '.join( [ str(i) for i in resid_list[1] ] ),
            )
        print(selectionstr)
        cmd = [
            "gmx", "select", "-f", self.gropath,
            "-s", self.tprpath, "-select", selectionstr,
            "-on", outputsel
        ]
        out, err = exec_gromacs(cmd)
        print(out, err)
        with open(outputsel, "r") as selectionf, open(add_grp_to, "a") as ndxf:
            for line in selectionf:
                ndxf.write(line)


def get_neighbor_dict(neighborfilename='neighbor_info',):
    ''' Returns a list of all neighbors being in the
        cutoff distance at least once in the trajectory.
        Neighborfile is required and is output of determine_neighbors()

        Dict layout is:
        neibdict[resid][time] --> [neibs]
    '''
    neibdict = {}
    reswithnoneib = []
    with open(neighborfilename,"r") as neibmap:
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
                reswithnoneib.append((resid, time))
                neibdict[resid].update({time:[]})
    LOGGER.info("No neighbors for following residues at time %s", reswithnoneib)
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
            data.append((residue, time, nneibs, neibindeces))
    return data


def dist_helper(v1, v2, box):
    img_vectors = [
        (0,  1),
        (1,  0),
        (1,  1),
        (0, -1),
        (-1, 0),
        (-1, -1),
        (1, -1),
        (-1, 1),
        (0, 0),
        ]
    ns = []
    for img in img_vectors:
        img = (*img, 0)
        newv = v2 + (box[:3]*img)
        ns.append(np.linalg.norm(v1-newv))
    return np.array(ns).min()
