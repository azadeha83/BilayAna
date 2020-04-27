'''
    This module stores all functions that are needed to calculate the interaction energy of lipids or its parts
'''
import os
import subprocess
from . import neighbors
from .. import log
from ..common import exec_gromacs, GMXNAME, write_submitfile
from ..systeminfo import SysInfo
from ..definitions import lipidmolecules
from ..command_line import submit_missing_energycalculation

LOGGER = log.LOGGER
LOGGER = log.create_filehandler("bilana_energy.log", LOGGER)


class Energy(SysInfo):
    '''
        This class stores all information that is needed to automate the calculate the lipid interaction energy
        The main function is
            run_calculation
        For more information see its docstring
    '''

    DENOMINATOR = 40
    LOGGER = LOGGER

    def __init__(self,
        part,
        inputfilename="inputfile",
        neighborfilename='neighbor_info',
        resindex_all='resindex_all',
        overwrite=True,
        verbosity="INFO",
        ):
        super().__init__(inputfilename)
        log.set_verbosity(verbosity)
        knownparts = ['complete', 'head-tail', 'head-tailhalfs', 'carbons']
        if part not in knownparts:
            raise ValueError("Part keyword specified is not known.")
        self.neiblist = neighbors.get_neighbor_dict()
        self.resindex_all = resindex_all
        self.overwrite = overwrite
        self.groupblocks = ()
        self.part = part
        if part == 'complete':
            self.molparts = ["resid_"]
            self.part = ''
            self.denominator = self.DENOMINATOR
            self.molparts_short = [""]
            if neighborfilename != 'neighbor_info':
                self.all_energies = 'all_energies_{}.dat'.format(neighborfilename)
            else:
                self.all_energies = 'all_energies.dat'
        elif part == 'head-tail':
            self.molparts = ["resid_h_", "resid_t_"]
            self.denominator = int(self.DENOMINATOR/2)
            self.molparts_short = ["h_", "t_"]
            self.all_energies='all_energies_headtail.dat'
        elif part == 'head-tailhalfs':
            self.molparts = ["resid_h_", "resid_t12_", "resid_t22_"]
            self.denominator = int(self.DENOMINATOR/4)
            self.molparts_short = ["h_","t12_","t22_"]
            self.all_energies = 'all_energies_headtailhalfs.dat'
        elif part == 'carbons':
            self.molparts = ['resid_C{}_'.format(i) for i in range(7)]
            self.denominator = int(self.DENOMINATOR/10)
            self.molparts_short = ['C{}_'.format(i) for i in range(7)]
            self.all_energies = "all_energies_carbons.dat"
        print('\n Calculating for energygroups:', self.molparts)

    def run_calculation(self, resids):
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
        LOGGER.info('Rerunning MD for energyfiles')
        for res in resids:
            LOGGER.info('Working on lipid %s ...', res)
            all_neibs_of_res = list(set([neibs for t in self.neiblist[res].keys() for neibs in self.neiblist[res][t]]))
            nneibs = len(all_neibs_of_res)
            if nneibs % self.denominator == 0:
                number_of_groupfragments = (nneibs//self.denominator)
            else:
                number_of_groupfragments = (nneibs//self.denominator)+1
            LOGGER.info("Needing %s energy run(s)", number_of_groupfragments)
            for groupfragment in range(number_of_groupfragments):
                LOGGER.info("... on fragment %s ...", groupfragment)

                g_energy_output = ''.join([\
                    self.energypath, '/xvgtables/energies_residue',\
                    str(res), '_', str(groupfragment), self.part, '.xvg',\
                    ])
                groupblockstart = groupfragment*self.denominator
                groupblockend = (groupfragment+1)*self.denominator
                self.groupblocks = (groupblockstart, groupblockend)

                # File in-/outputs
                groupfragment=str(groupfragment)
                mdpout = ''.join([self.energypath, 'mdpfiles/energy_mdp_recalc_resid', str(res), '_', groupfragment, self.part, '.mdp'])
                tprout = ''.join([self.energypath, 'tprfiles/mdrerun_resid', str(res), '_', groupfragment, self.part, '.tpr'])
                energyf_output = ''.join([self.energypath, 'edrfiles/energyfile_resid', str(res), '_'+groupfragment, self.part, '.edr'])
                xvg_out = ''.join([self.energypath, 'xvgtables/energies_residue', str(res), '_', groupfragment, self.part, '.xvg'])
                energygroups = self.gather_energygroups(res, all_neibs_of_res)
                relev_energies = self.get_relev_energies(res, all_neibs_of_res)

                # Run functions
                self.create_MDP(mdpout, energygroups)
                self.create_TPR(mdpout, tprout)
                if os.path.isfile(energyf_output) and not self.overwrite:
                    LOGGER.info("Edrfile for lipid %s part %s already exists. Will skip this calculation.", res, groupfragment)
                else:
                    self.do_Energyrun(res, groupfragment, tprout, energyf_output)
                if os.path.isfile(g_energy_output) and not self.overwrite:
                    LOGGER.info("Xvgtable for lipid %s part %s already exists. Will skip this calculation.", res, groupfragment)
                else:
                    self.write_XVG(energyf_output, tprout, relev_energies, xvg_out)
        return 1

    def run_lip_leaflet_interaction(self, resids):
        ''' '''
        LOGGER.info('Rerunning MD for energyfiles')
        self.groupblocks = [None, None]
        for res in resids:
            LOGGER.info('Working on lipid %s ...', res)
            leaflet = self.res_to_leaflet[res]

            res_other_leaflet = []
            for nres in self.MOLRANGE:
                leaf = self.res_to_leaflet[nres]
                if leaf != leaflet:
                    res_other_leaflet.append(nres)

            outputndx = self.temppath + "/energy_leaflet_res{}.ndx".format(res)
            self.resindex_all = outputndx
            selectionstr = 'System=all; leaflet=resname {0} and resid {1}; resid_{2}= resid {2}; resid_{2}; leaflet; System;'.format(
                ' '.join(self.molecules),
                ' '.join([str(i) for i in res_other_leaflet]),
                res,
                )

            print(selectionstr)
            cmd = [
                "gmx", "select", "-f", self.gropath,
                "-s", self.tprpath, "-select", selectionstr,
                "-on", outputndx,
            ]
            out, err = exec_gromacs(cmd)
            print(out, err)

            g_energy_output = ''.join([\
                self.energypath, '/xvgtables/energies_residue',\
                str(res), '_leaflet', '.xvg',\
                ])
            mdpout = ''.join([self.energypath, 'mdpfiles/energy_mdp_recalc_resid', str(res), '_leaflet', '.mdp'])
            tprout = ''.join([self.energypath, 'tprfiles/mdrerun_resid', str(res), "_leaflet", '.tpr'])
            energyf_output = ''.join([self.energypath, 'edrfiles/energyfile_resid', str(res), '_leaflet', '.edr'])
            xvg_out = ''.join([self.energypath, 'xvgtables/energies_residue', str(res), '_leaflet', '.xvg'])
            energygroups = "resid_{} leaflet".format(res)
            relev_energies = '\n'.join( ["Coul-SR:resid_{}-leaflet".format(res), "LJ-SR:resid_{}-leaflet".format(res), '\n'] )

            # Run functions
            self.create_MDP(mdpout, energygroups)
            self.create_TPR(mdpout, tprout)
            if os.path.isfile(energyf_output) and not self.overwrite:
                LOGGER.info("Edrfile for lipid %s part %s already exists. Will skip this calculation.", res)
            else:
                self.do_Energyrun(res, 0, tprout, energyf_output)
            if os.path.isfile(g_energy_output) and not self.overwrite:
                LOGGER.info("Xvgtable for lipid %s part %s already exists. Will skip this calculation.", res)
            else:
                self.write_XVG(energyf_output, tprout, relev_energies, xvg_out)
        return 1

    def create_lipid_water_interaction_file(self):
        ''' Create a file with entries of
            interaction of resid at time to solvent
            <Time> <resid> <resname> <Etot> <Evdw> <Ecoul>
        '''

    def create_lipid_leaflet_interaction_file(self):
        ''' Create a file with entries of
            interaction of resid at time to leaflet0 and leaflet1
            <Time> <resid> <resname> <host_leaflet> <leaflet_index> <Etot> <Evdw> <Ecoul>
        '''

    def gather_selfinteractions(self):
        '''
            Create selfinteraction.dat file from existing *.edr calculation
            !!! Attention !!! Can only be done after actual energy run.
        '''
        self.selfinteractions_edr_to_xvg()
        self.selfinteractions_xvg_to_dat()

    def selfinteractions_edr_to_xvg(self):
        ''' Extracts all self interaction energy values from .edr files using gmx energy '''
        missing_edr = []
        for res in self.MOLRANGE:
            relev_energies = self.get_relev_self_interaction(res)
            tprout = ''.join([self.energypath, 'tprfiles/mdrerun_resid', str(res), '_', '0', self.part, '.tpr'])
            energyf_output = ''.join([self.energypath, 'edrfiles/energyfile_resid', str(res), '_'+'0', self.part, '.edr'])
            if not os.path.isfile(energyf_output):
                missing_edr.append(energyf_output)
                continue
            xvg_out = ''.join([self.energypath, 'xvgtables/energies_residue', str(res), '_selfinteraction', self.part, '.xvg'])
            self.write_XVG(energyf_output, tprout, relev_energies, xvg_out)
        if missing_edr:
            raise FileNotFoundError("Following files are missing: {}".format(missing_edr))

    def selfinteractions_xvg_to_dat(self):
        ''' Extracts all self interaction energy entries from xvg files
            and writes them to "selfinteractions.dat"
        '''
        with open("selfinteractions.dat", "w") as energyoutput:
            print(\
                  '{: <10}{: <10}{: <10}'
                  '{: <20}{: <20}{: <20}{: <20}{: <20}{: <20}{: <20}'\
                  .format("Time", "resid", "resname",
                          "Etot", "VdWSR", "CoulSR", "VdW14", "Coul14", "VdWtot", "Coultot", ),
                  file=energyoutput)

            for resid in self.MOLRANGE:
                xvg_out = ''.join([self.energypath, 'xvgtables/energies_residue', str(resid), '_selfinteraction', self.part, '.xvg'])
                restype = self.resid_to_lipid[resid]
                with open(xvg_out,"r") as xvgfile:
                    res_to_row = {}

                    for energyline in xvgfile: #folderlayout is: time Coul_resHost_resNeib LJ_resHost_resNeib ...
                        energyline_cols = energyline.split()

                        if '@ s' in energyline:                     #creating a dict to know which column(energies) belong to which residue
                            row = int(energyline_cols[1][1:])+1                 #time is at row 0 !
                            host = energyline_cols[3].split("resid_")[1][:-1]
                            energytype = energyline_cols[3].split(":")[0][1:]
                            #print("Host: {} Type: {} Row: {}".format(host, energytype, row))
                            res_to_row.update({(energytype, host):row})

                        elif '@' not in energyline and '#' not in energyline:     #pick correct energies from energyfile and print
                            time = float(energyline_cols[0])
                            if time % self.dt != 0:
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
                            print(
                                  '{: <10}{: <10}{: <10}{: <20.5f}'
                                  '{: <20}{: <20}{: <20}{: <20}{: <20.5f}{: <20.5f}'
                                  .format(time, resid, restype, Etot,
                                          vdw_sr, coul_sr, vdw_14, coul_14, vdw_tot, coul_tot,),
                                  file=energyoutput)

    def gather_energygroups(self, res, all_neibs_of_res):
        ''' Set which part of molecule should be considered '''
        energygroup_resids = [res] + all_neibs_of_res[ self.groupblocks[0]:self.groupblocks[1] ]
        energygroup_list = []
        for resid in energygroup_resids:
            if self.resid_to_lipid[resid] in lipidmolecules.STEROLS+lipidmolecules.PROTEINS:
                energygroup_list.append(''.join(["resid_",str(resid)]))
            else:
                for part in self.molparts:
                    energygroup_list.append(''.join([part,str(resid)]))
            # Add leaflet and water groups
        #leaflet = self.res_to_leaflet[res]
        #if leaflet:
        #    energygroup_list += ["leaflet0", "solv"]
        #else:
        #    energygroup_list += ["leaflet1", "solv"]
        energygroup_list += ["solv"]

        energygroup_string = ' '.join(energygroup_list)
        return energygroup_string

    def get_relev_energies(self, res, all_neibs_of_res):
        '''
            Returns string that describes all entries
            needed to be extracted from energy file using gmx energy
            This version is for lipid-lipid interaction
            for self interaction search function "get_relev_self_interaction"
        '''
        Etypes=["Coul-SR:", "LJ-SR:"]
        energyselection=[]
        for interaction in Etypes:
            counterhost = 0 #for cholesterol as it has just 1 molpart

            for parthost in self.molparts:

                # If lipid is a sterol take just one part, --> counterhost=0 -> write | counter>0 -> skip
                if self.resid_to_lipid[res] in (lipidmolecules.STEROLS + lipidmolecules.PROTEINS) and counterhost == 0:
                    parthost="resid_"
                    counterhost += 1
                elif self.resid_to_lipid[res] in (lipidmolecules.STEROLS + lipidmolecules.PROTEINS) and counterhost != 0:
                    continue

                for neib in all_neibs_of_res[self.groupblocks[0]:self.groupblocks[1]]:
                    counterneib = 0

                    for partneib in self.molparts:
                        if self.resid_to_lipid[neib] in (lipidmolecules.STEROLS + lipidmolecules.PROTEINS) and counterneib == 0:
                            partneib='resid_'
                            counterneib+=1
                        elif self.resid_to_lipid[neib] in (lipidmolecules.STEROLS + lipidmolecules.PROTEINS) and counterneib != 0: ## Because self.part is changed to empty, for complete case 0:
                            continue

                        energyselection.append(''.join([interaction, parthost, str(res), "-", partneib,str(neib)]))

        all_relev_energies = '\n'.join(energyselection+['\n'])
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
        ''' Create mdpfile '''
        os.makedirs(self.energypath+'/mdpfiles', exist_ok=True)
        with open(mdpout,"w") as mdpfile_rerun:
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
            raw_mdp.append('ref_t = '+str(self.temperature))
            energygrpline = ''.join(['energygrps\t\t\t=', energygroups, '\n'])
            raw_mdp.append(energygrpline)
            mdpfile_rerun.write('\n'.join(raw_mdp)+'\n')

    def create_TPR(self, mdpoutfile: str, tprout: str):
        ''' Create TPRFILE with GROMPP '''
        os.makedirs(self.energypath+'tprfiles', exist_ok=True)
        grompp_arglist = [GMXNAME, 'grompp', '-f', mdpoutfile, '-p',\
                        self.toppath, '-c', self.gropath, '-o', tprout,\
                        '-n', self.resindex_all, '-po', mdpoutfile\
                        ]
        out, err = exec_gromacs(grompp_arglist)
        with open("gmx_grompp.log","a") as logfile:
            logfile.write(err)
            logfile.write(out)

    def do_Energyrun(self, res, groupfragment, tprrerun_in, energyf_out):
        ''' Create .edr ENERGYFILE with mdrun -rerun '''
        LOGGER.info('...Rerunning trajectory for energy calculation...')
        os.makedirs(self.energypath+'edrfiles', exist_ok=True)
        os.makedirs(self.energypath+'logfiles', exist_ok=True)
        logoutput_file = self.energypath+'logfiles/'+'mdrerun_resid'+str(res)+self.part+'frag'+str(groupfragment)+'.log'
        trajout = 'EMPTY.trr' # As specified in mdpfile, !NO! .trr-file should be written
        
        mdrun_arglist = [GMXNAME, 'mdrun', '-s', tprrerun_in, '-rerun', self.trjpath,
                    '-e', energyf_out, '-o', trajout,'-g', logoutput_file,
                    ]    
        # if self.ff == 'all_atom':
        #     mdrun_arglist = [GMXNAME, 'mdrun', '-s', tprrerun_in, '-rerun', self.trjpath,
        #                 '-e', energyf_out, '-o', trajout,'-g', logoutput_file,
        #                 ]    
        # else:
        #     print('coarse')
        #     mdrun_arglist = [GMXNAME, 'mdrun', '-s', tprrerun_in, '-rerun', self.trjpath_energy,
        #                     '-e', energyf_out, '-o', trajout,'-g', logoutput_file,
        #                     ]
        out, err = exec_gromacs(mdrun_arglist)
        with open("gmx_mdrun.log","a") as logfile:
            logfile.write(err)
            logfile.write(out)

    def write_XVG(self, energyf_in, tprrerun_in, all_relev_energies, xvg_out):
        ''' Create XVG-TABLE with all relevant energies '''
        os.makedirs(self.energypath+'xvgtables', exist_ok=True)
        g_energy_arglist=[GMXNAME, 'energy', '-f', energyf_in,
                          '-s', tprrerun_in,'-o', xvg_out,
                          ]
        inp_str=all_relev_energies.encode()
        out, err = exec_gromacs(g_energy_arglist, inp_str)
        with open("gmx_energy.log","a") as logfile:
            logfile.write(err)
            logfile.write(out)

    def write_energyfile(self, submit_missing_data=True):
        ''' Creates files: "all_energies_<interaction>.dat
            NOTE: This function is too long. It should be separated into smaller parts.
        '''
        missing_energydata = []
        LOGGER.info('Create energy file')
        with open(self.all_energies, "w") as energyoutput:
            print(
                  '{: <10}{: <10}{: <10}{: <20}'
                  '{: <20}{: <20}{: <20}'\
                  .format("Time", "Host", "Neighbor", "Molparts",\
                                           "VdW", "Coul", "Etot"),\
                  file=energyoutput)

            for resid in self.MOLRANGE:
                LOGGER.info("Working on residue %s ...", resid)

                # Get neighborhood of resid
                residtype        = self.resid_to_lipid[resid]
                all_neibs_of_res = list(set([neibs for t in self.neiblist[resid].keys() for neibs in self.neiblist[resid][t]]))
                n_neibs          = len(all_neibs_of_res)
                LOGGER.debug("All neibs of res %s are %s", resid, all_neibs_of_res)

                # Get number of fragments (how many runs per residue)
                if n_neibs % self.denominator == 0:
                    number_of_groupfragments = (n_neibs//self.denominator)
                else:
                    number_of_groupfragments = (n_neibs//self.denominator)+1
                LOGGER.debug("Nneibs: %s Nfrags: %s", n_neibs, number_of_groupfragments)

                all_neibs_of_res = [ all_neibs_of_res[ ( i*self.denominator ):( (i+1)*self.denominator ) ] for i in range(number_of_groupfragments) ]

                for part in range(number_of_groupfragments):
                    LOGGER.debug("At part %s", part)
                    xvgfilename = self.energypath+'xvgtables/energies_residue'+str(resid)+'_'+str(part)+self.part+'.xvg'

                    with open(xvgfilename,"r") as xvgfile:
                        res_to_rowindex = {}

                        for energyline in xvgfile: #folderlayout is: <time> <Coul_resHost_resNeib> <LJ_resHost_resNeib> ...
                            energyline_cols = energyline.split()

                            if '@ s' in energyline:                     #creating a dict to know which column(energies) belong to which residue
                                rowindex  = int(energyline_cols[1][1:])+1 # time is at row 0 !
                                neib = energyline_cols[3].split("resid_")[2][:-1]
                                host = energyline_cols[3].split("resid_")[1][:-1]
                                energytype = energyline_cols[3].split("-")[0][1:]
                                LOGGER.debug("Hostid: %s, Neibid: %s", host, neib)

                                res_to_rowindex[(energytype, host, neib)] = rowindex
                                LOGGER.debug("Adding to dict: Etype %s, host %s, neib %s", energytype, host, neib)

                            elif '@' not in energyline and '#' not in energyline: #pick correct energies from energyfile and print
                                time = float(energyline_cols[0])

                                for neib in all_neibs_of_res[part]:

                                    # This if clause is due to a broken simulation... In future it should be removed
                                    if self.system == 'dppc_dupc_chol25' and ((int(host) == 372 and neib == 242) or (int(host) == 242 and neib == 372)):
                                        continue

                                    neibtype = self.resid_to_lipid[neib]

                                    host_sterol_processed = False # We need that as sterols are not separated into different parts
                                    # so we just process sterols only once
                                    for parthost in self.molparts:
                                        parthost = parthost[6:] # remove "resid" from parthost string

                                        if residtype in lipidmolecules.STEROLS and not host_sterol_processed:
                                            parthost = ''
                                            host_sterol_processed = True
                                        elif residtype in lipidmolecules.STEROLS and host_sterol_processed:
                                            continue

                                        neib_sterol_processed = False
                                        for partneib in self.molparts:
                                            partneib = partneib[6:] # remove "resid" from parthost string

                                            if neibtype in lipidmolecules.STEROLS and not neib_sterol_processed:
                                                partneib = ''
                                                neib_sterol_processed = True
                                            elif neibtype in lipidmolecules.STEROLS and neib_sterol_processed:
                                                continue

                                            # if partstring is empty take whole lipid
                                            if parthost.replace("_", "") == '':
                                                interhost = 'w'
                                            else:
                                                interhost = parthost.replace("_", "")

                                            if partneib.replace("_", "")== '':
                                                interneib = 'w'
                                            else:
                                                interneib = partneib.replace("_", "")
                                            inter = ''.join([interhost, '_', interneib])

                                            # Get energies from dict energyline_cols
                                            try:
                                                vdw  = energyline_cols[ res_to_rowindex[ ('LJ',   parthost+str(resid), partneib+str(neib)) ] ]
                                                coul = energyline_cols[ res_to_rowindex[ ('Coul', parthost+str(resid), partneib+str(neib)) ] ]
                                            except KeyError:
                                                LOGGER.warning("Data not found for %s - %s", parthost+str(resid), partneib+str(neib) )
                                                if resid not in missing_energydata:
                                                    missing_energydata.append(resid)

                                            Etot = float(vdw) + float(coul)
                                            print(\
                                                  '{: <12.1f}{: <10}{: <10}{: <20}'
                                                  '{: <20}{: <20}{: <20.5f}'
                                                  .format(time, resid, neib, inter,
                                                                            vdw, coul, Etot),
                                                  file=energyoutput)

        if missing_energydata:
            LOGGER.warning("Missing energydata: %s", missing_energydata)
            if not self.part:
                self.part = "complete"
            if submit_missing_data:
                for missing_resid in missing_energydata:
                    submit_missing_energycalculation(missing_resid, self.part, self.system, self.temperature)
            if os.path.isfile(self.all_energies):
                os.remove(self.all_energies)
            raise RuntimeError("There were inconsistencies in the data. See log files for further information.")

        LOGGER.info("File %s written successfully", self.all_energies)

    def check_exist_xvgs(self, check_len=False):
        ''' Checks if all .xvg-files containing lipid interaction exist

            check_len can be set to simulation length and all .xvg files that differ from that length are included to missing_xvgfile.info
        '''
        def read_lastline_only(fname):
            ''' Reads last line of file without loop. Attention: Crashes if file is empty'''
            MAXCOUNT = 100000
            cnt = 0
            with open(fname, "rb") as f:
                f.seek(-2, os.SEEK_END)     # Jump to the second last byte.
                while f.read(1) != b"\n" and f.tell() != 1:   # Until EOL is found... if f.tell() == 0 means cursor is at the beginning of file
                    f.seek(-2, os.SEEK_CUR) # ...jump back the read byte plus one more.
                    cnt += 1
                    if cnt > MAXCOUNT:
                        raise RuntimeError("Reach max byte count in file {}. Is it corrupted?".format(fname))

                if f.tell() == 1:
                    raise RuntimeError("File {} has no EOL character. Is it corrupted?".format(fname))

                last = f.readline()         # Read last line.
            return last

        missing_files = []

        time_ok = True
        missing_res = False
        for resid in self.MOLRANGE:
            all_neibs_of_res = list(set([neibs for t in self.neiblist[resid].keys() for neibs in self.neiblist[resid][t]]))
            n_neibs = len(all_neibs_of_res)

            if n_neibs % self.denominator == 0:
                number_of_groupfragments = (n_neibs//self.denominator)
            else:
                number_of_groupfragments = (n_neibs//self.denominator)+1

            for part in range(number_of_groupfragments):
                xvgfilename = self.energypath+'/xvgtables/energies_residue'+str(resid)+'_'+str(part)+self.part+'.xvg'

                if not os.path.isfile(xvgfilename):
                    # Check wether file exists
                    missing_res = True
                    missing_files.append(xvgfilename)

                elif os.stat(xvgfilename).st_size == 0:
                    # Check wether file is empty
                    missing_files.append(xvgfilename)

                elif check_len:
                    # Check wether whole trajectory was used
                    line = read_lastline_only(xvgfilename)
                    time = float(line.decode().split()[0])
                    if time != check_len:
                        time_ok = False
                        missing_files.append(xvgfilename)

        if missing_files:
            resubmitted_residues = [] # Need this list, because can only resubmit all files(/fragment) for residue
            for fname in missing_files:
                res = int(fname.split("/")[-1].replace("energies_residue", "").split("_")[0])
                if res not in resubmitted_residues:

                    if not self.part: ## Because self.part is changed to empty, for complete case
                        part = 'complete'
                    else:
                        part = self.part

                    submit_missing_energycalculation(res, part, self.system, self.temperature)
                    LOGGER.debug("Would submit %s", res)

                    resubmitted_residues.append(res)

                else:
                    continue
            LOGGER.warning("Submitted following energy calculations for residues: %s", resubmitted_residues)
            LOGGER.warning("Issues found:\nThere were missing residues: %s\nLength of trajectory not the same as indicated in inputfile: %s", missing_res, not time_ok)

            return False
        return True
