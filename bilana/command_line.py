import os, sys
import subprocess
import re
from bilana.systeminfo import SysInfo
#from src import gromacstoolautomator as gmx
#from src import mainanalysis

def write_submitfile(submitout, jobname, ncores=2, mem='4G', prio=False):
    username = os.environ['LOGNAME']
    if not prio:
        queue = 'short'
    else:
        queue = 'prio'
    with open(submitout, "w") as sfile:
        print('#!/bin/bash'
              '\n#SBATCH -A q0heuer'
              '\n#SBATCH -p {queue}'
              '\n#SBATCH --output={jobname}.out'
              '\n#SBATCH --mail-type=fail'
              '\n#SBATCH --mail-user={username}@wwu.de'
              '\n#SBATCH --time=2-00:00:00'
              '\n#SBATCH --ntasks=1'
              '\n#SBATCH --nodes=1'
              '\n#SBATCH --cpus-per-task={ncores}'
              '\n#SBATCH --mem={mem}'
              '\nexport OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK'\
              '\nsrun $@'.format(**locals()), file=sfile)

def write_jobscript(scriptname, lipidpart, start_res, end_res, overwrite=True):
    with open(scriptname, "w") as jobf:
        print('import os, sys\n'
                'from bilana import gromacstoolautomator as gmx\n'
                'from bilana.systeminfo import SysInfo \n'
                'mysystem = SysInfo("inputfile")\n'
                'myenergy = gmx.Energy(mysystem, "'+str(lipidpart)+'",overwrite='+str(overwrite)+')\n'
                'myenergy.run_calculation(startres='+str(start_res)+', endres='+str(end_res)+')\n'
                'os.remove(sys.argv[0])\n',\
            file=jobf)


def get_minmaxdiv(startdiv, numerator, direction=-1):
    print(startdiv)
    startdiv = int(startdiv)
    if startdiv < 2 and direction == 1:
        startdiv += 1
    #elif startdiv > numerator//2 and direction == 1:
    #    raise ValueError("Divisor became too high.")
    if direction == 1:
        rangelist = [startdiv, 2*startdiv, 1]
    elif direction == -1:
        rangelist = [startdiv, startdiv//2, -1]
    else:
        raise ValueError("Wrong input for direction. Choose either -1 or 1.")
    if numerator % startdiv != 0:
        for new_div in range(*rangelist):
            if numerator % new_div == 0:
                print("found it", new_div)
                return new_div
        new_startdiv = int(startdiv*(2*direction))
        print('No value found. Restarting with {} as start value'.format(new_startdiv))
        return get_minmaxdiv(new_startdiv, numerator, direction)
    else:
        return startdiv

def submit_energycalcs():
    ''' Divide energyruns into smaller parts for faster computation '''
    startdir = os.getcwd()
    if len(sys.argv) < 7:
        print('Invalid number of input arguments. Specify:\n'
              '<systemname> <lowest T> <highest T> <part of lipid> <jobname> <max divisor> <overwrite: on/off>(opt.)')
        sys.exit()
    systemname = sys.argv[1]
    tempstart = int(sys.argv[2])
    tempend = int(sys.argv[3])
    lipidpart = sys.argv[4]
    jobname = lipidpart+'_'+sys.argv[5]
    startdivisor = sys.argv[6]
    overwrite = True
    try:
        if sys.argv[7] == 'off':
            overwrite = False
        elif sys.argv[7] == 'on':
            overwrite = True
        else:
            raise ValueError("Value for overwrite key invalid. Choose either 'off' or 'on' ")
    except IndexError:
        pass
    Temperatures = [T for T in range(tempstart, tempend+1, 10)]
    systems_to_calculate_for = ['./{}_{}'.format(systemname, T) for T in Temperatures]
    for systemdir in systems_to_calculate_for:
        os.chdir(systemdir)
        mysystem = SysInfo('inputfile')
        systemsize = mysystem.number_of_lipids
        divisor = get_minmaxdiv(startdivisor, systemsize)
        lipids_per_part = systemsize//divisor
        for jobpart in range(divisor):
            start_res = str((jobpart*lipids_per_part)+1)
            end_res = str(((jobpart+1)*lipids_per_part))
            jobfile_name = str(start_res)+'-'+str(end_res)+'_'+jobname
            jobscript_name = 'exec_energycalc'+str(jobfile_name)+'.py'
            write_jobscript(jobscript_name, lipidpart, start_res, end_res, overwrite=overwrite)
            write_submitfile('submit.sh', jobfile_name)
            cmd = ['sbatch', '-J', systemdir[2:]+'_'+str(start_res)+'-'+str(end_res)+'_'+jobname, 'submit.sh','python3', jobscript_name]
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = proc.communicate()
            print(out.decode(), err.decode())
        os.chdir(startdir)

def initialize_system():
    ''' Creates all core files like neighbor_info, resindex_all, scd_distribution '''
    startdir = os.getcwd()
    if len(sys.argv) < 4 or len(sys.argv) > 6:
        print('Invalid number of input arguments. Specify:\n'
              '<systemname> <lowest T> <highest T> <jobname> <refatoms=P>')
        sys.exit()
    systemname = sys.argv[1]
    tempstart = int(sys.argv[2])
    tempend = int(sys.argv[3])
    jobname = sys.argv[4]
    try:
        refatoms = sys.argv[5]
    except IndexError:
        refatoms = 'P'
    temperatures = [T for T in range(tempstart, tempend+1, 10)]
    systems_to_calculate_for = ['./{}_{}'.format(systemname, T) for T in temperatures]
    for systemdir in systems_to_calculate_for:
        os.chdir(systemdir)
        scriptfilename = 'exec'+systemdir[2:]+jobname+'.py'
        jobfilename = systemdir[2:]+jobname
        with open(scriptfilename, 'w') as scriptf:
            print(\
                'import os, sys'
                '\nimport subprocess'
                '\nfrom bilana import gromacstoolautomator as gmx'
                '\nfrom bilana import mainanalysis'
                '\nfrom bilana.systeminfo import SysInfo'
                '\nmysystem = SysInfo("inputfile")'
                '\ngmx.Neighbors(mysystem).determine_neighbors(refatoms="{0}")'
                '\ngmx.Neighbors(mysystem).create_indexfile()'
                '\ngmx.produce_gro(mysystem)'
                '\nmainanalysis.Scd(mysystem).create_scdfile()'
                '\nmainanalysis.create_leaflet_assignment_file(mysystem)'
                '\nif os.path.isfile("initialize.sh"):'
                '\n    subprocess.call("./initialize.sh")'
                '\nos.remove(sys.argv[0])'.format(refatoms),\
                file=scriptf)
            write_submitfile('submit.sh', jobfilename, mem='16G', prio=False)
            cmd = ['sbatch', '-J', jobfilename, 'submit.sh','python3', scriptfilename]
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = proc.communicate()
            print(out.decode(), err.decode())
        os.chdir(startdir)


def check_and_write():
    ''' Check if all energy files exist and write table with all energies '''
    startdir = os.getcwd()
    if len(sys.argv) < 6:
        print('Invalid number of input arguments. Specify:\n'
              '<systemname> <lowest T> <highest T> <part of lipid> <jobname> <neibfile>')
        sys.exit()
    systemname = sys.argv[1]
    tempstart = int(sys.argv[2])
    tempend = int(sys.argv[3])
    lipidpart = sys.argv[4]
    jobname = lipidpart+'_'+sys.argv[5]
    try:
        neibfilename = sys.argv[6]
    except IndexError:
        neibfilename = 'neighbor_info'
    temperatures = [T for T in range(tempstart, tempend+1, 10)]
    systems_to_calculate_for = ['./{}_{}'.format(systemname, T) for T in temperatures]
    for systemdir in systems_to_calculate_for:
        os.chdir(systemdir)
        temperature = systemdir[-3:]
        scriptfilename = 'exec'+systemdir[2:]+jobname+'.py'
        jobfilename = systemdir[2:]+jobname
        with open(scriptfilename, 'w') as scriptf:
            print(\
                'import os, sys'
                '\nimport subprocess'
                '\nfrom bilana import gromacstoolautomator as gmx'
                '\nfrom bilana.systeminfo import SysInfo'
                '\nfrom bilana import energyfilecreator as ef'
                '\nmysystem = SysInfo("inputfile")'
                #'\nmyenergystate = gmx.Energy(mysystem,"'+lipidpart+',neighborfile_name="neigbor_info"')'
                '\nmyenergystate = gmx.Energy(mysystem,"'+lipidpart+'",neighborfile_name="'+neibfilename+'")'
                '\nif myenergystate.check_exist_xvgs():'
                '\n    myenergystate.write_energyfile()'
                '\n    efstate = ef.EofScd(mysystem,"'+lipidpart+'", myenergystate.all_energies, "scd_distribution.dat")'
                '\n    efstate.create_eofscdfile()'
                '\nelse:'
#                '\n    raise ValueError("There are .edr files missing.")'
                '\n    print("Submitting missing energy calculations")'
                '\n    os.chdir("'+startdir+'")'
                '\n    cmd = ["submit_energycalcs", "'+systemdir[2:-4]+'", "'+temperature+'",'
                                '"'+temperature+'", "'+lipidpart+'", "en", "40", "off"]'
                '\n    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)'
                '\n    out, err = proc.communicate()'
                '\n    proc.wait()'
                '\n    proc.stdout.close()'
                '\n    proc.stderr.close()'
                '\n    print(out.decode(), err.decode())'
                '\n    os.chdir("'+systemdir+'")'
                '\nos.remove(sys.argv[0])',\
                file=scriptf)
            write_submitfile('submit.sh', jobfilename, mem='16G')
            cmd = ['sbatch', '-J', jobfilename, 'submit.sh','python3', scriptfilename]
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = proc.communicate()
            print(out.decode(), err.decode())
        os.chdir(startdir)

def write_eofscd():
    ''' Check if all energy files exist and write table with all energies '''
    startdir = os.getcwd()
    if len(sys.argv) != 6:
        print('Invalid number of input arguments. Specify:\n'
              '<systemname> <lowest T> <highest T> <part of lipid> <jobname> <neibfilename>')
        sys.exit()
    systemname = sys.argv[1]
    tempstart = int(sys.argv[2])
    tempend = int(sys.argv[3])
    lipidpart = sys.argv[4]
    jobname = lipidpart+'_'+sys.argv[5]
    try:
        neibfilename = sys.argv[6]
    except IndexError:
        neibfilename = 'neighbor_info'
    Temperatures = [T for T in range(tempstart, tempend+1, 10)]
    systems_to_calculate_for = ['./{}_{}'.format(systemname, T) for T in Temperatures]
    for systemdir in systems_to_calculate_for:
        os.chdir(systemdir)
        temperature = systemdir[-3:]
        scriptfilename = 'exec'+systemdir[2:]+jobname+'.py'
        jobfilename = systemdir[2:]+jobname
        with open(scriptfilename, 'w') as scriptf:
            print(\
                'import os, sys'
                '\nimport subprocess'
                '\nfrom bilana import gromacstoolautomator as gmx'
                '\nfrom bilana.systeminfo import SysInfo'
                '\nfrom bilana import energyfilecreator as ef'
                '\nmysystem = SysInfo("inputfile")'
                '\nmyenergystate = gmx.Energy(mysystem,"'+lipidpart+'",neighborfile_name="'+neibfilename+'")'
                '\nif myenergystate.check_exist_xvgs():'
#                '\n    myenergystate.write_energyfile()'
                '\n    efstate = ef.EofScd(mysystem,"'+lipidpart+'", myenergystate.all_energies, "scd_distribution.dat")'
                '\n    efstate.create_eofscdfile()'
                '\nelse:'
#                '\n    raise ValueError("There are .edr files missing.")'
                '\n    print("Submitting missing energy calculations")'
                '\n    os.chdir("'+startdir+'")'
                '\n    cmd = ["submit_energycalcs", "'+systemdir[2:-4]+'", "'+temperature+'",'
                                '"'+temperature+'", "'+lipidpart+'", "en", "40", "off"]'
                '\n    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)'
                '\n    out, err = proc.communicate()'
                '\n    proc.wait()'
                '\n    proc.stdout.close()'
                '\n    proc.stderr.close()'
                '\n    print(out.decode(), err.decode())'
                '\n    os.chdir("'+systemdir+'")'
                '\nos.remove(sys.argv[0])',\
                file=scriptf)
            write_submitfile('submit.sh', jobfilename, mem='16G', prio=True)
            cmd = ['sbatch', '-J', jobfilename, 'submit.sh','python3', scriptfilename]
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = proc.communicate()
            print(out.decode(), err.decode())
        os.chdir(startdir)

def write_nofscd():
    ''' Check if all energy files exist and write table with all energies '''
    startdir = os.getcwd()
    if len(sys.argv) < 5:
        print('Invalid number of input arguments. Specify:\n'
              '<systemname> <lowest T> <highest T> <jobname> <neighborfile>')
        sys.exit()
    systemname = sys.argv[1]
    tempstart = int(sys.argv[2])
    tempend = int(sys.argv[3])
    jobname = sys.argv[4]
    try:
        neibfile_name = sys.argv[5]
    except IndexError:
        neibfile_name = 'neighbor_info'
    temperatures = [T for T in range(tempstart, tempend+1, 10)]
    systems_to_calculate_for = ['./{}_{}'.format(systemname, T) for T in temperatures]
    for systemdir in systems_to_calculate_for:
        os.chdir(systemdir)
        scriptfilename = 'exec'+systemdir[2:]+jobname+'.py'
        jobfilename = systemdir[2:]+jobname
        with open(scriptfilename, 'w') as scriptf:
            print(\
                '\nimport os, sys'
                '\nfrom bilana.systeminfo import SysInfo'
                '\nfrom bilana import energyfilecreator as ef'
                '\nmysystem = SysInfo("inputfile")'
                '\nefstate = ef.NofScd(mysystem)'
                '\nefstate.create_NofScd_input("scd_distribution.dat", "'+neibfile_name+'")'
                '\nos.remove(sys.argv[0])',\
                file=scriptf)
            write_submitfile('submit.sh', jobfilename, mem='8G', prio=True)
            cmd = ['sbatch', '-J', jobfilename, 'submit.sh','python3', scriptfilename]
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = proc.communicate()
            print(out.decode(), err.decode())
        os.chdir(startdir)
