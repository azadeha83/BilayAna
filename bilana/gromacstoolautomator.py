''' Module for automating analysis using Gromacs '''

import subprocess
import re
import os
import sys 
from time import localtime, strftime

from bilana import lipidmolecules

from bilana.systeminfo import mysystem
global mysystem

#inputfile = sys.argv[1]
#mysystem = systeminfo.SysInfo(inputfile)


gmx_exec = 'gmx' #'gmx_5.0.4_mpi'
os.environ["GMX_MAXBACKUP"] = "-1"

def exec_gromacs(cmd,inp_str=None): 
    '''arglist (cmd) is list of arguments like "['gmx cmd','-f','tprfile','-e','en.edr']"
        inp_str must be encoded!
    '''
    if inp_str is None:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = proc.communicate()
    else:
        proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = proc.communicate(inp_str)
    proc.wait()
    proc.stdout.close()
    proc.stderr.close()
    return out, err

def trajectory_to_gro(overwrite='off', atomlist=None, lipids='all'):
    os.makedirs(mysystem.datapath+'/grofiles')
    getcoords = {}
    print('Converting trajectory-file to structure-file...\n')
    if lipids != 'all':
        molecules = [lipids]
    else:
        if "CHOL" in mysystem.molecules:
            molecules_new = mysystem.molecules.copy()
            molecules_new.remove("CHOL")
            #print(molecules_new)
        else:
            molecules_new = mysystem.molecules
        molecules = molecules_new
    for lipidmolecule in molecules:
        grofile_output = ''.join([mysystem.datapath, '/grofiles', '/calc_scd_for', str(lipidmolecule), '.gro'])
        sys.stdout.flush()
        if atomlist == None:
            atomlist = lipidmolecules.tail_atoms_of[lipidmolecule]+['P', 'O3']
        print(strftime("%H:%M:%S :", localtime()),"Processing {} ...".format(lipidmolecule))
        inp_str = str(lipidmolecule).encode()
        if os.path.isfile(grofile_output) and overwrite=='off':
            pass
        else:
            gmx_traj_arglist = [\
                gmx_exec,'trjconv','-s',mysystem.tprpath, '-f',mysystem.trjpath,\
                '-o',grofile_output,\
                '-b', str(mysystem.t_start), '-e', str(mysystem.t_end),\
                '-dt', str(mysystem.dt),\
                '-pbc', 'whole',\
                ]
            out,err = exec_gromacs(gmx_traj_arglist,inp_str)      #Creates a gro file containing {timeframe:resid:atoms:xyz}
            with open("gmx_traj.log","w") as logfile:
                logfile.write(err.decode())
                logfile.write(150*'_')
                logfile.write(out.decode())
                logfile.write(150*'_')
        with open(grofile_output,"r") as grofile:
            print(strftime("%H:%M:%S :", localtime()),"...read data output...")
            regexp = re.compile(r'[\s]*\d+'+lipidmolecule)
            for line in grofile:
                if 't=' in line:
                    time = float(line[line.index('t=')+2:].strip())   #to get rid of obsolete decimals
                    print("...at time {}".format(time), end="\r")
                    if float(mysystem.t_end) < time:
                        print("breaking at", time)
                        break
                #print("Match",regexp.match(line),line[:15],end='\r')
                #print("Time match is",float(self.t_start)<=time,end='\r')
                if float(mysystem.t_start) <= time and regexp.match(line) != None:
                    #print("Reading data at",time,end='\r')
                    atom = line[9:15].strip()
                    lipidtype = line[5:9]
                    if atom not in atomlist and lipidtype not in molecules: continue
                    resid = line[:5].strip()
                    coordinates = [float(x) for x in line[20:44].split()]
                    keytuple = (str(lipidtype),time,int(resid),str(atom))
                    getcoords.update({keytuple:coordinates})
    return getcoords, time-1000.0

def radialdistribution(self):
    print("\n_____Calculating radial distribution function ____\n")
    os.makedirs(self.datapath+'/rdf', exist_ok=True)
    selectstringP = 'name P;\nname P;\nname O3;'
    selectstringO3 = 'name O3;\nname P;\nname O3;'
    selectionfileP = self.temppath+'/selrdfP'
    selectionfileO3 = self.temppath+'/selrdfO3'
    outputfileP = self.datapath+'/rdf'+'/rdfhostP.xvg'
    outputfileO3 = self.datapath+'/rdf'+'/rdfhostO3.xvg'
    print("...preparing selectionfiles...")
    with open(selectionfileP,"w") as sfP,open(selectionfileO3,"w") as sfO3:
        sfP.write(selectstringP)
        sfO3.write(selectstringO3)
    g_rdf_arglistP = [gmx_exec,'-nobackup','rdf','-f',self.trjpath,'-s',\
                      self.tprpath,'-sf',selectionfileP,'-o',\
                      outputfileP,'-xy',\
                      ]
    print("...first selection...")
    out, err = self.exec_gromacs(g_rdf_arglistP)
    with open("gmx_rdfP.log","w") as logfile:
        logfile.write(err.decode())
        logfile.write(150*'_')
        logfile.write(out.decode())
        logfile.write(150*'_')
    g_rdf_arglistO3 = [gmx_exec,'-nobackup','rdf','-f',self.trjpath,\
                       '-s',self.tprpath,'-sf',selectionfileO3,'-o',\
                       outputfileO3,'-xy',\
                       ]
    print("...second selection...")
    out, err = self.exec_gromacs(g_rdf_arglistO3)
    with open("gmx_rdfO3.log","w") as logfile:
        logfile.write(err.decode())
        logfile.write(150*'_')
        logfile.write(out.decode())
        logfile.write(150*'_')

def calculate_distance(self):
    print("\n____Calculating distances____\n")
    neiblist=self.get_neighbor_dict()
    os.makedirs(self.datapath+'/distcalc', exist_ok=True)
    for i in range(1,self.NUMBEROFPARTICLES+1):
        all_N_of_res=list(set([neibs for t in neiblist[i].keys() for neibs in neiblist[i][t]]))
        sys.stdout.flush()
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
            sys.stdout.flush()
            print("Working on residue {}".format(resid), end="\r")
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

    #global mysystem

    def create_selectionfile_neighborsearch(self,resid):
        filename="{}/neighbors_of_residue{}".format(mysystem.temppath, resid)
        with open(filename,"w") as selection:
            print(\
                'host =  resid {} and (name P O3);\n'\
                'allOAtoms = resname CHL1 and name O3 and not host;\n'\
                'allPAtoms = resname DPPC DUPC and name P and not host;\n'\
                'neibOs = allOAtoms and within 1.0 of host;\n'\
                'neibPs = allPAtoms and within 1.0 of host;\n'\
                'neibs = neibPs or neibOs;\n'\
                'neibs;'\
                .format(resid), file=selection)
        return filename

    def get_neighbor_dict(self,verbose='off'):
        ''' Returns a list of all neighbors being in the
         cutoff distance at least once in the trajectory.
        Neighborfile is !required! and is output of 
        "determine_neighbors()" '''
        neighborfile = "neighbor_info"
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
        for i in range(1, mysystem.NUMBEROFMOLECULES+1):
        #for i in range(1,2):
            print("Working on residue {}".format(i), end='\r')
            selectionfile = mysystem.temppath+'/tmp_selectionfile'
            with open(selectionfile,"w") as sf:
                lipidtype = mysystem.resid_to_lipid[i]
                if lipidtype  != 'CHL1': 
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
                    carbons = ['C{} C{} C{} C{}'.format(tailatomlist[0][i],\
                                                        tailatomlist[0][i+1],\
                                                        tailatomlist[1][i],\
                                                        tailatomlist[1][i+1])\
                               for i in range(0, len(tailatomlist), 2)]
                    selprefixes = [("", r'".*"'),\
                                   ("h", headatoms),\
                                   ("t", tailatoms),\
                                   ("t12", tailhalf12),\
                                   ("t22", tailhalf22),\
                                   ("C1", carbons[0]),\
                                   ("C2", carbons[1]),\
                                   ("C3", carbons[2]),\
                                   ("C4", carbons[3]),\
                                   ("C5", carbons[4]),\
                                   ("C6", carbons[5]),\
                                   ("C7", carbons[6]),\
                                   ("C8", carbons[7]),\
                                   ]
                    selectionlist = []
                    for item in selprefixes:
                        selectionlist += 'resid_{0}{1}=resid {1} and resname {2} and name {3};\n'\
                                        .format(item[0], str(i), lipidtype, item[1])
                    
                    lastlineitems = ['resid_{}{};\n'.format(item[0], str(i)) for item in selprefixes]
                    selectionlist += ''.join(lastlineitems)
                    selectionstring=''.join(selectionlist)
                    #selectionlist += [''.join(["resid_h_",\
                    #               str(i),"=resid ",str(i)," and resname ",\
                    #               lipidtype," and name ",' '.join(headatoms),\
                    #               ";\n"]),\
                    #               ]
                    sf.write(selectionstring)
                elif mysystem.resid_to_lipid[i] == 'CHL1':
                    selectionstring = 'resid_{0}=resid {0}  and resname CHL1;\nresid_{0};'.format(str(i))
                    sf.write(selectionstring)
            outputindex = mysystem.indexpath+"/resid_"+str(i)+".ndx"
            gmx_select_arglist = [gmx_exec, 'select', '-s', mysystem.gropath, '-sf',\
                                  selectionfile, '-on', outputindex,\
                                  ]
            out, err = exec_gromacs(gmx_select_arglist)
            with open("gmx_select.log","w") as logfile: 
                logfile.write(err.decode())
                logfile.write(150*'_')
                logfile.write(out.decode())
                logfile.write(150*'_')
            ##append resid.ndx to resindex_all.ndx
            with open(outputindex,"r") as output_index:
                filecontent = output_index.readlines()
                resindex_all.write(''.join(filecontent)+'\n\n')
        ### To have whole system indices in one group    
        make_ndx_output = mysystem.temppath+'/make_ndx_system.ndx' 
        gmx_make_ndx_arglist = [gmx_exec, 'make_ndx', '-f', mysystem.gropath, '-o', make_ndx_output] 
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

    def determine_neighbors(self, overwrite=True):
        print("\n____Determining neighbors____\n")
        os.makedirs(mysystem.datapath+'/neighborfiles', exist_ok=True)
        with open("neighbor_info","w") as outfile:
            outfile.write('Resid \t Time \t Number_of_neighbors \t List_of_Neighbors \n')
            for residue in range(1,mysystem.NUMBEROFMOLECULES+1):
                print(". . . Working on residue: {} . . .".format(residue), end="\r")
                sys.stdout.flush()
                selectionfile = self.create_selectionfile_neighborsearch(residue)
                indexoutput = '{}/neighbors_of_residue{}.ndx'.format(mysystem.indexpath, residue)
                datafileoutput = '{}/neighborfiles/neighbors_of_residue{}.dat'.format(mysystem.datapath, residue)
                if os.path.isfile(datafileoutput) and overwrite == False:
                    print("Neighbor file of residue {} already exists. Skipping.".format(residue))
                else:
                    cmdlist=[\
                        gmx_exec,'select','-s',mysystem.tprpath,'-f',mysystem.trjpath,\
                        '-sf',selectionfile,'-on',indexoutput,'-oi',datafileoutput,\
                        '-b',str(mysystem.t_start), '-e', str(mysystem.t_end),\
                        '-dt', str(mysystem.dt),\
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
                        neibresid = [mysystem.index_to_resid[x] for x in neibindeces]
                        residlist = ','.join([str(x) for x in neibresid])
                        print('{} \t {} \t {} \t {}'.format(residue,time,nneibs,residlist), file=outfile)

class Energy():
    ''' All about energy calculation '''

    DENOMINATOR = 40

    def __init__(self, resindex_all, mdp_raw, overwrite=True, startres=1, endres=-1, parts='complete'):
        #self.mysystem = mysystem#
        knownparts = ['complete', 'head-tail', 'head-tailhalfs', 'carbons']
        if parts not in knownparts:
            raise ValueError("Part keyword specified is not known.")
        self.neiblist = Neighbors().get_neighbor_dict()
        self.resindex_all = resindex_all
        self.mdp_raw = mdp_raw
        self.overwrite = overwrite
        self.groupblocks = ()
        self.startres = startres
        if endres == -1:
            self.endres = mysystem.NUMBEROFMOLECULES
        if parts == 'complete':
            self.molparts = ["resid_"]
            self.parts = ''
            self.denominator = self.DENOMINATOR
            self.molparts_short = [""]
            self.all_energies = 'all_energies'
            self.interactions = ['']
        elif parts == 'head-tail':
            self.molparts = ["resid_h_","resid_t_"]
            self.parts = parts
            self.denominator = int(self.DENOMINATOR/2)
            self.molparts_short = ["h_", "t_"]
            self.interactions = ['head-tail', 'head-head', 'tail-tail']
            self.all_energies='all_energies_headtail.dat'
        elif parts == 'head-tailhalfs':
            self.molparts = ["resid_h_", "resid_t12_", "resid_t22_"]
            self.parts = parts
            self.denominator = int(self.DENOMINATOR/4)
            self.molparts_short = ["h_","t12_","t22_"]
            self.interactions = ['head-tail12', 'tail12-tail12', 'head-tail22', 'tail22-tail22']
            self.all_energies = 'all_energies_headtailhalfs.dat'
        elif parts == 'carbons':
            self.molparts = ['resid_C{}_'.format(i) for i in range(lipidmolecules.shortestchain/2)]
            self.parts = parts
            self.denominator = int(self.DENOMINATOR/10)
            self.molparts_short = ['C{}_'.format(i) for i in range(lipidmolecules.shortestchain/2)]
            self.interactions = ['C{0}-C{0}'.format(i) for i in range(lipidmolecules.shortestchain/2)]
            self.all_energies = ["all_energies_carbons.dat"]
        print('\n Calculating for energygroups:', self.molparts)

    def run_calculation(self):
        ''' Runs a complete energy calculation with settings from Energy() instance '''
        print('''\n____Rerunning MD for energyfiles,
         yielding xvgtables with relevant energies.____\n
         Caution mdp-file must not have energy_grps indicated!\n''')
        for res in range(self.startres,self.endres+1): 
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
                sys.stdout.flush()
                g_energy_output = ''.join([\
                    mysystem.energypath, '/xvgtables/energies_residue',\
                    str(res), '_', str(groupfragment), self.parts, '.xvg',\
                    ])
                if os.path.isfile(g_energy_output) and self.overwrite == False:
                    print("Xvgtable for lipid {} part {} already exists. Will skip this calculation.".format(res, groupfragment))
                    continue
                groupblockstart = groupfragment*self.denominator
                groupblockend = (groupfragment+1)*self.denominator
                self.groupblocks = (groupblockstart, groupblockend)
                # File in-/outputs
                groupfragment=str(groupfragment) 
                mdpout = ''.join([mysystem.energypath, '/mdpfiles/energy_mdp_recalc_resid', str(res), '_', groupfragment, self.parts, '.mdp'])
                tprout = ''.join([mysystem.energypath, '/tprfiles/mdrerun_resid', str(res), '_', groupfragment, self.parts, '.tpr'])
                energyf_output = ''.join([mysystem.energypath, '/edrfiles/energyfile_resid', str(res), '_'+groupfragment, self.parts, '.edr'])
                xvg_out = ''.join([mysystem.energypath, '/xvgtables/energies_residue', str(res), '_', groupfragment, self.parts, '.xvg'])
                energygroups = self.gather_energygroups(res, all_neibs_of_res)
                relev_energies = self.get_relev_energies(res, all_neibs_of_res)
                # Run functions
                self.create_MDP(self.mdp_raw, mdpout, energygroups)
                self.create_TPR(mdpout, tprout)
                self.do_Energyrun(res, groupfragment, tprout, energyf_output)
                self.write_XVG(energyf_output, tprout, relev_energies, xvg_out)
        return 'Done'

    def gather_energygroups(self, res, all_neibs_of_res):
        ''' Set which part of molecule should be considered '''
        energygroup_indeces=[res]+all_neibs_of_res[self.groupblocks[0]:self.groupblocks[1]]
        energygroup_list=[]
        for index in energygroup_indeces:
            if mysystem.resid_to_lipid[index]=='CHL1':
                energygroup_list.append(''.join(["resid_",str(index)]))
            else:
                for part in self.molparts:
                    energygroup_list.append(''.join([part,str(index)]))
        energygroup_string=' '.join(energygroup_list)
        return energygroup_string

    def get_relev_energies(self, res, all_neibs_of_res):
        Etypes=["Coul-SR:","LJ-SR:"]
        energyselection=[]
        for interaction in Etypes:
            counterhost=0 #for the cholesterol as it has just 1 molpart
            for parthost in self.molparts:
                if mysystem.resid_to_lipid[res]=='CHL1' and counterhost==0:
                    parthost="resid_"
                    counterhost+=1
                elif mysystem.resid_to_lipid[res]=='CHL1' and counterhost!=0:
                    continue
                for neib in all_neibs_of_res[self.groupblocks[0]:self.groupblocks[1]]:
                    counterneib=0
                    for partneib in self.molparts:
                        if mysystem.resid_to_lipid[neib]=='CHL1' and counterneib==0:
                            partneib='resid_'
                            counterneib+=1
                        elif mysystem.resid_to_lipid[neib]=='CHL1' and counterneib!=0:
                            continue
                        energyselection.append(''.join([interaction,parthost,str(res),"-",partneib,str(neib)]))
                #select_energies_coulomb='\n'.join(["Coul-SR:resid_"+str(res)+"-resid_"+str(x) for x in all_N_of_res[groupblockstart:groupblockend]])
                #select_energies_LJ='\n'.join(["LJ-SR:resid_"+str(res)+"-resid_"+str(x) for x in all_N_of_res[groupblockstart:groupblockend]])
                #select_all_relevant_energies=select_energies_coulomb+"\n"+select_energies_LJ
        all_relev_energies='\n'.join(energyselection+['\n'])
        return all_relev_energies

    def create_MDP(self, mdp_raw: str, mdpout: str, energygroups: str):
        ''' Create Mdpfile '''
        os.makedirs(mysystem.energypath+'/mdpfiles', exist_ok=True)
        with open(mdpout,"w") as mdpfile_rerun, open(mdp_raw,"r") as mdpfile_raw:
            mdp_raw_content = mdpfile_raw.readlines()
            energygrpline = ''.join(['energygrps\t\t\t=', energygroups, '\n'])
            mdp_raw_content.append(energygrpline)
            mdpfile_rerun.write('\n'.join(mdp_raw_content)+'\n')

    def create_TPR(self, mdpoutfile: str, tprout: str):
        ''' Create TPRFILE with GROMPP '''
        print(strftime("%H:%M:%S :", localtime()),'...Creating .tpr-file...')
        os.makedirs(mysystem.energypath+'/tprfiles', exist_ok=True)
        grompp_arglist=[gmx_exec, 'grompp', '-f', mdpoutfile, '-p',\
                        mysystem.toppath, '-c', mysystem.gropath, '-o', tprout,\
                        '-n', self.resindex_all, '-po', mdpoutfile\
                        ]
        out, err = exec_gromacs(grompp_arglist)
        with open("gmx_grompp.log","a") as logfile:
            logfile.write(err.decode())
            logfile.write(100*'_')
            logfile.write(out.decode())
            logfile.write(100*'_')

    def do_Energyrun(self, res, groupfragment, tprrerun_in, energyf_out):
        ''' Create ENERGYFILE with mdrun -rerun '''
        print(strftime("%H:%M:%S :", localtime()),'...Rerunning trajectory for energy calculation...')
        os.makedirs(mysystem.energypath+'/edrfiles', exist_ok=True)
        os.makedirs(mysystem.energypath+'/logfiles', exist_ok=True)
        logoutput_file = mysystem.energypath+'/logfiles'+'/mdrerun_resid'+str(res)+self.parts+'frag'+groupfragment+'.log'
        trajout = 'EMPTY.trr' # As specified in mdpfile, NO .trr-file should be written
        mdrun_arglist = [gmx_exec, 'mdrun', '-s', tprrerun_in,'-rerun', mysystem.trjpath,\
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
        os.makedirs(mysystem.energypath+'/xvgtables', exist_ok=True)
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
        with open(self.all_energies,"w") as energyoutput:
            print(\
                  '{: <10}{: <10}{: <10}{: <20}'
                  '{: <20}{: <20}{: <20}'\
                  .format("Time", "Host", "Neighbor", "Molparts",\
                                           "VdW", "Coul", "Etot"),\
                  file=energyoutput)
            for resid in range(1, mysystem.NUMBEROFMOLECULES+1):
                print("Working on residue {}".format(resid), end="\r")
                residtype = mysystem.resid_to_lipid[resid]
                all_neibs_of_res = list(set([neibs for t in self.neiblist[resid].keys() for neibs in self.neiblist[resid][t]]))
                n_neibs = len(all_neibs_of_res)
                if n_neibs % self.denominator == 0:
                    number_of_groupfragments = (n_neibs//self.denominator)
                else:
                    number_of_groupfragments = (n_neibs//self.denominator)+1
                #print(number_of_groupfragments)
                for part in range(number_of_groupfragments):
                    groupblockstart = part*self.denominator
                    groupblockend = (part+1)*self.denominator
                    neighbors_part_are = all_neibs_of_res[groupblockstart:groupblockend]
                    xvgfilename = mysystem.energypath+'/xvgtables/energies_residue'+str(resid)+'_'+str(part)+self.parts+'.xvg'
                    with open(xvgfilename,"r") as xvgfile:
                        res_to_row = {}
                        for energyline in xvgfile: #folderlayout is: time Coul_resHost_resNeib LJ_resHost_resNeib ...
                            energyline_cols = energyline.split()
                            if '@ s' in energyline:                     #creating a dict to know which column(energies) belong to which residue
                                row = int(energyline_cols[1][1:])+1                 #time is at row 0 !
                                neib = energyline_cols[3].split("resid_")[2][:-1]
                                host = energyline_cols[3].split("resid_")[1][:-1]
                                energytype = energyline_cols[3].split("-")[0][1:]
                                res_to_row.update({(energytype, host, neib):row})
                            elif '@' not in energyline and '#' not in energyline:     #pick correct energies from energyfile and print
                                time = float(energyline_cols[0])
                                for neib in neighbors_part_are:
                                    neibtype = mysystem.resid_to_lipid[neib]
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
                                            vdw = energyline_cols[res_to_row[('LJ', parthost+str(resid), partneib+str(neib))]]
                                            coul = energyline_cols[res_to_row[('Coul', parthost+str(resid), partneib+str(neib))]]
                                            Etot = float(vdw)+float(coul)
                                            print(\
                                                  '{: <10},{: <10},{: <10},{: <20},'
                                                  '{: <20},{: <20},{: <20.5f}'
                                                  .format(time, resid, neib, inter,\
                                                                            vdw, coul, Etot),\
                                                  file=energyoutput)

    def check_exist_xvgs(self):
        all_okay = True
        for resid in range(1, mysystem.NUMBEROFMOLECULES+1):
            all_neibs_of_res = list(set([neibs for t in self.neiblist[resid].keys() for neibs in self.neiblist[resid][t]]))
            n_neibs = len(all_neibs_of_res)
            if n_neibs % self.denominator == 0:
                number_of_groupfragments = (n_neibs//self.denominator)
            else:
                number_of_groupfragments = (n_neibs//self.denominator)+1
            for part in range(number_of_groupfragments):
                xvgfilename = mysystem.energypath+'/xvgtables/energies_residue'+str(resid)+'_'+str(part)+self.parts+'.xvg'
                if not os.path.isfile(xvgfilename):
                    print('File is missing:', xvgfilename)
                    all_okay = False
        print('All okay?', all_okay)
