import os, sys
import subprocess
import re
from bilana.systeminfo import SysInfo
from bilana import gromacstoolautomator as gmx
from bilana import mainanalysis

def write_submitfile(submitout, jobname, ncores=2, mem='4G'):
    with open(submitout, "w") as sfile:
        print('#!/bin/bash'
              '\n#SBATCH -A q0heuer'
              '\n#SBATCH -p short'
              '\n#SBATCH --output='+jobname+'.out'
#              '\n#SBATCH --error='+jobname+'.err'
              '\n#SBATCH --mail-type=fail'
              '\n#SBATCH --mail-user=f_kell07@wwu.de'
              '\n#SBATCH --time=2-00:00:00'
              '\n#SBATCH --ntasks=1'
              '\n#SBATCH --nodes=1'
              '\n#SBATCH --cpus-per-task='+str(ncores)+
              '\n#SBATCH --mem='+str(mem)+
              '\nexport OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK'\
              '\nsrun $@', file=sfile)

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
    mysystem = SysInfo('inputfile')
    gmx.Neighbors(mysystem).determine_neighbors()
    gmx.Neighbors(mysystem).create_indexfile()
    mainanalysis.Scd(mysystem).create_scdfile()
    #gmx.trajectory_to_gro(mysystem) Already is included in create_scdfile

def mend_energyruns():
    ''' Run missing parts of energyruns '''
    def find_exceeded_jobs(dir_to_search):
        exceeded_jobs = []
        regex_filenames = re.compile(r'^(\d+)-(\d+)_(\w+)_en\.err$')
        for file in os.listdir(dir_to_search):
            regmatch = regex_filenames.match(file)
            if regmatch:
                startlipid = regmatch.group(1)
                endlipid = regmatch.group(2)
                lipidpart = regmatch.group(3)
                with open(file, "r") as errfile:
                    for line in errfile:
                        if 'DUE TO TIME LIMIT' in line:
                            exceeded_jobs.append((file[:-3]+'out', lipidpart, startlipid, endlipid))
                            break
        return exceeded_jobs
    def find_last_lipid_calculated(filename):
        calculated_lipids = []
        regex_lipidline = re.compile(r'Working on lipid (\d+)')
        with open(filename, "r") as jobf:
            for line in jobf:
                regexmatch = regex_lipidline.search(line)
                if regexmatch:
                    calculated_lipids.append(regexmatch.group(1))
        return calculated_lipids[-1]
    startdir = os.getcwd()
    if len(sys.argv) != 5:
        print('Invalid number of input arguments. Specify:\n'
              '<systemname> <lowest T> <highest T> <jobname>')
        sys.exit()
    systemname = sys.argv[1]
    tempstart = int(sys.argv[2])
    tempend = int(sys.argv[3])
    Temperatures = [T for T in range(tempstart, tempend+1, 10)]
    systems_to_calculate_for = ['./{}_{}'.format(systemname, T) for T in Temperatures]
    for systemdir in systems_to_calculate_for:
        print(systemdir)
        os.chdir(systemdir)
        exceeded_jobs = find_exceeded_jobs(os.getcwd())
        for job in exceeded_jobs:
            lipidpart = job[1]
            jobname = lipidpart+'_'+'rest'+sys.argv[4]
            jobfile = job[0]
            correct_last_lipid = int(job[3])
            actual_last_lipid = int(find_last_lipid_calculated(jobfile))
            Nmissing_lipids = correct_last_lipid-actual_last_lipid+1
            Nlipids_calculated = actual_last_lipid-int(job[2]) 
            #print("Last lipid was {}. Need to calculate energy for {} lipids and calculated {}. {}".format(actual_last_lipid, Nmissing_lipids, Nlipids_calculated, jobs))
            if Nmissing_lipids > Nlipids_calculated:
                Nruns = get_minmaxdiv(Nmissing_lipids//Nlipids_calculated, Nmissing_lipids, 1)
            else:
                Nruns = 1
            #endres = actual_last_lipid
            lipids_per_run = Nmissing_lipids//Nruns
            for start_res in range(actual_last_lipid, correct_last_lipid+1, lipids_per_run):
                #start_res = endres+jobpart*lipids_per_run
                end_res = start_res+lipids_per_run-1
                print(start_res, end_res)
                jobfile_name = str(start_res)+'-'+str(end_res)+'_'+jobname
                jobscript_name = 'exec_energycalc'+str(jobfile_name)+'.py'
                write_jobscript(jobscript_name, lipidpart, start_res, end_res)
                write_submitfile('submit.sh', jobfile_name)
                cmd = ['sbatch', '-J', systemdir[2:]+'_'+str(start_res)+'-'+str(end_res)+'_'+jobname, 'submit.sh','python3', jobscript_name]
                proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                out, err = proc.communicate()
                print(out.decode(), err.decode())
        os.chdir(startdir)

#===============================================================================
# def temporal_rerun():
#     def find_fkn_resten(dir_to_search):
#         restens = []
#         regex_filenames = re.compile(r'^exec_energycalc(\d+)-(\d+)_(\w+)_resten\.py$')
#         for file in os.listdir(dir_to_search):
#             regmatch = regex_filenames.match(file)
#             if regmatch:
#                 startlipid = regmatch.group(1)
#                 endlipid = regmatch.group(2)
#                 lipidpart = regmatch.group(3)
#                 restens.append((file[:-3]+'out', lipidpart, startlipid, endlipid))
#         return restens
#     startdir = os.getcwd()
#     if len(sys.argv) != 5:
#         print('Invalid number of input arguments. Specify:\n'
#               '<systemname> <lowest T> <highest T> <jobname>')
#         sys.exit()
#     systemname = sys.argv[1]
#     tempstart = int(sys.argv[2])
#     tempend = int(sys.argv[3])
#     Temperatures = [T for T in range(tempstart, tempend+1, 10)]
#     systems_to_calculate_for = ['./{}_{}'.format(systemname, T) for T in Temperatures]
#     for systemdir in systems_to_calculate_for:
#         print(systemdir)
#         os.chdir(systemdir)
#         exceeded_jobs = find_fkn_resten(os.getcwd())
#         for job in exceeded_jobs:
#             lipidpart = job[1]
#             jobname = lipidpart+'_'+'ulti'+sys.argv[4]
#             jobfile = job[0]
#             start_res  = int(job[2])
#             end_res = int(job[3])
#             jobfile_name = str(start_res)+'-'+str(end_res)+'_'+jobname
#             jobscript_name = 'exec_energycalc'+str(jobfile_name)+'.py'
#             write_jobscript(jobscript_name, lipidpart, start_res, end_res)
#             write_submitfile('submit.sh', jobfile_name)
#             cmd = ['sbatch', '-J', systemdir[2:]+'_'+str(start_res)+'-'+str(end_res)+'_'+jobname, 'submit.sh','python3', jobscript_name]
#             proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#             out, err = proc.communicate()
#             print(out.decode(), err.decode())
#         os.chdir(startdir)
#===============================================================================

def check_and_write():
    ''' Check if all energy files exist and write table with all energies '''
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
                '\nmyenergystate = gmx.Energy(mysystem,"'+lipidpart+'")'
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
            write_submitfile('submit.sh', jobfilename, mem='32G')
            cmd = ['sbatch', '-J', jobfilename, 'submit.sh','python3', scriptfilename]
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = proc.communicate()
            print(out.decode(), err.decode())
        os.chdir(startdir)


