#!/bin/env python3

import os, sys
import subprocess
from bilana.systeminfo import SysInfo

divisor = 10
startdir = os.getcwd()
systemname = sys.argv[1]
tempstart = sys.argv[2]
tempend = sys.argv[3]
lipidpart = sys.argv[4]
jobname = sys.argv[5]
Temperatures = [T for T in range(tempstart, tempend+1, 10)]
systems_to_calculate_for = ['./{}_{}'.format(systemname,T) for T in Temperatures]

for systemdir in systems_to_calculate_for:
    os.chdir(systemdir)
    mysystem = SysInfo('inputfile')
    size = mysystem.determine_systemsize_and_number_of_lipids()
    lipids_per_partfloor=size//divisor
    for jobpart in range(divisor):
        start_res = str(jobpart*lipids_per_partfloor)
        end_res = str((jobpart+1)*lipids_per_partfloor)
        jobscriptname = 'exec_energycalc'+str(jobpart)+'.py'
        with open(jobscriptname, "w") as jobf:
            print('import os, sys\n'
                    'from bilana import gromacstoolautomator as gmx\n'
                    'myenergy = gmx.Energy("resindex_all", "energy_recalculation.mdp", '
                        'startres='+str(start_res)+', endres='+str(end_res)+', parts="'+str(lipidpart)+'")\n'
                    'myenergy.run_calculation()\n'
                    'os.remove(sys.argv[0])\n',\
                file=jobf)
        cmd=['sbatch','-J',systemdir+'_'+jobpart+jobname,'submit.sh','python',jobscriptname]
        proc=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        out, err = proc.communicate()
        print(out.decode(),err.decode())
    lipids_per_partmod=size % divisor
    if lipids_per_partmod != 0:
        with open(jobscriptname, "w") as jobf:
            print('import os, sys\n'
                    'from bilana import gromacstoolautomator as gmx\n'
                    'myenergy = gmx.Energy("resindex_all", "energy_recalculation.mdp", '
                        'startres='+str(end_res+1)+', endres='+str(end_res+1+lipids_per_partmod)+', '
                        'parts="'+str(lipidpart)+'")\n'
                    'myenergy.run_calculation()\n'
                    'os.remove(sys.argv[0]\n',\
                file=jobf)
    cmd=['sbatch','-J',systemdir+'_-1'+jobname,'submit.sh','python',jobscriptname]
    proc=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    out, err = proc.communicate()
    print(out.decode(),err.decode())
    os.chdir(startdir)
    
    
    
