import os, sys
import subprocess
from bilana.systeminfo import SysInfo
from bilana import gromacstoolautomator as gmx
from bilana import mainanalysis

def write_submitfile(submitout, jobname, ncores=2, mem='2G'):
    with open(submitout, "w") as sfile:
        print('#!/bin/bash'
              '\n#SBATCH -A q0heuer'
              '\n#SBATCH -p short'
              '\n#SBATCH --output='+jobname+'.out'
              '\n#SBATCH --error='+jobname+'.err'
              '\n#SBATCH --mail-type=fail'
              '\n#SBATCH --mail-user=f_kell07@wwu.de'
              '\n#SBATCH --time=2-00:00:00'
              '\n#SBATCH --ntasks=1'
              '\n#SBATCH --nodes=1'
              '\n#SBATCH --cpus-per-task='+str(ncores)+
              '\n#SBATCH --mem='+str(mem)+
              '\nexport OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK'\
              '\nsrun $@', file=sfile)

def submit_energycalcs():
    ''' Divides energy run in smaller parts for faster computation '''
    divisor = 40
    startdir = os.getcwd()
    if len(sys.argv) != 6:
        print('Invalid number of input arguments. Specify:\n'
              '<systemname> <lowest T> <highest T> <part of lipid> <jobname>')
        sys.exit()
    systemname = sys.argv[1]
    tempstart = int(sys.argv[2])
    tempend = int(sys.argv[3])
    lipidpart = sys.argv[4]
    jobname = lipidpart+'_'+sys.argv[5]
    Temperatures = [T for T in range(tempstart, tempend+1, 10)]
    systems_to_calculate_for = ['./{}_{}'.format(systemname, T) for T in Temperatures]
    for systemdir in systems_to_calculate_for:
        os.chdir(systemdir)
        mysystem = SysInfo('inputfile')
        size = mysystem.number_of_lipids
        lipids_per_partfloor = size//divisor
        for jobpart in range(divisor):
            start_res = str((jobpart*lipids_per_partfloor)+1)
            end_res = str(((jobpart+1)*lipids_per_partfloor))
            job_filenames = str(start_res)+'-'+str(end_res)+'_'+jobname
            jobscriptname = 'exec_energycalc'+str(job_filenames)+'.py'
            with open(jobscriptname, "w") as jobf:
                print('import os, sys\n'
                        'from bilana import gromacstoolautomator as gmx\n'
                        'from bilana.systeminfo import SysInfo \n'
                        'mysystem = SysInfo("inputfile")\n'
                        'myenergy = gmx.Energy("resindex_all", "energy_recalculation.mdp", mysystem, parts="'+str(lipidpart)+'")\n'
                        'myenergy.run_calculation(startres='+str(start_res)+', endres='+str(end_res)+')\n'
                        'os.remove(sys.argv[0])\n',\
                    file=jobf)
            write_submitfile('submit.sh', job_filenames)
            cmd=['sbatch', '-J', systemdir[2:]+'_'+str(start_res)+'-'+str(end_res)+'_'+jobname, 'submit.sh','python', jobscriptname]
            proc=subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = proc.communicate()
            print(out.decode(), err.decode())
        lipids_per_partmod = size % divisor
        if lipids_per_partmod != 0:
            start_res = str(int(end_res)+1)
            end_res = str(mysystem.number_of_lipids)
            job_filenames = str(start_res)+'-'+str(end_res)+'_'+jobname
            jobscriptname = 'exec_energycalc'+str(job_filenames)+'.py'
            with open(jobscriptname, "w") as jobf:
                print('import os, sys\n'
                        'from bilana import gromacstoolautomator as gmx\n'
                        'from bilana.systeminfo import SysInfo \n'
                        'mysystem = SysInfo("inputfile")\n'
                        'myenergy = gmx.Energy("resindex_all", "energy_recalculation.mdp", mysystem, parts="'+str(lipidpart)+'")\n'
                        'myenergy.run_calculation(startres='+str(start_res)+', endres='+str(end_res)+')\n'
                        'os.remove(sys.argv[0])\n',\
                    file=jobf)
            write_submitfile('submit.sh', job_filenames)
            cmd = ['sbatch', '-J', systemdir[2:]+'_'+str(start_res)+'-'+str(end_res)+'_'+jobname, 'submit.sh', 'python', jobscriptname]
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = proc.communicate()
            print(out.decode(), err.decode())
        os.chdir(startdir)
        
def initialize_system():
    ''' Creates all core files like neighbor_info, resindex_all, scd_distribution '''
    mysystem = SysInfo('inputfile')
    gmx.Neighbors(mysystem).determine_neighbors()
    gmx.Neighbors(mysystem).create_indexfile()
    mainanalysis.Scd(mysystem).create_scdfile()
    #gmx.trajectory_to_gro(mysystem) Already is included in create_scdfile
    
    