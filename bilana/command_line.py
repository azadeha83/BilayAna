import os
import subprocess
from .common import write_submitfile
from .common import get_minmaxdiv
from .systeminfo import SysInfo



def submit_energycalcs(systemname, temperature, jobname, lipidpart, *args,
    inputfilename="inputfile",
    neighborfile="neighbor_info",
    startdivisor=80,
    overwrite=True,
    **kwargs,):
    ''' Divide energyruns into smaller parts for faster computation '''
    complete_name = './{}_{}'.format(systemname, temperature)
    os.chdir(complete_name)
    mysystem = SysInfo(inputfilename)
    systemsize = mysystem.number_of_lipids
    divisor = get_minmaxdiv(startdivisor, systemsize)
    if divisor % 1 != 0:
        raise ValueError("divisor must be int")
    resids_to_calculate = mysystem.MOLRANGE
    lipids_per_part = systemsize//divisor
    for jobpart in range(divisor):
        list_of_res = resids_to_calculate[jobpart*lipids_per_part:(jobpart+1)*lipids_per_part]
        jobfile_name = str(jobpart)+'_'+jobname
        jobscript_name = 'exec_energycalc'+str(jobfile_name)+'.py'
        with open(jobscript_name, "w") as jobf:
            print(
                '\nimport os, sys'
                '\nfrom bilana.analysis.energy import Energy'
                '\nenergy_instance = Energy("{0}", overwrite="{1}", inputfilename="{2}", neighborfilename="{3}")'
                '\nenergy_instance.run_calculation(resids={4})'
                '\nos.remove(sys.argv[0])'.format(lipidpart, overwrite, inputfilename, neighborfile, list_of_res),
                file=jobf)
        write_submitfile('submit.sh', jobfile_name)
        cmd = ['sbatch', '-J', complete_name[2:]+str(jobpart)+'_'+jobname, 'submit.sh','python3', jobscript_name]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = proc.communicate()
        print(out.decode(), err.decode())


def initialize_system(systemname, temperature, jobname, *args,
    inputfilename="inputfile",
    refatoms="P",
    **kwargs):
    ''' Creates all core files like neighbor_info, resindex_all, scd_distribution '''
    complete_systemname = './{}_{}'.format(systemname, temperature)
    os.chdir(complete_systemname)
    scriptfilename = 'exec'+complete_systemname[2:]+jobname+'.py'
    jobfilename = complete_systemname[2:]+jobname
    with open(scriptfilename, 'w') as scriptf:
        print(\
            'import os, sys'
            '\nimport bilana'
            '\nfrom bilana import analysis'
            '\nfrom bilana.analysis.neighbors import Neighbors'
            '\nfrom bilana.analysis.order import Order'
            '\nneib_inst = Neighbors(inputfilename="{0}")'
            '\norder_inst = Order(inputfilename="{0}")'
            '\nsysinfo_inst = bilana.SysInfo(inputfilename="{0}")'
            '\nneib_inst.determine_neighbors(refatoms="{1}", overwrite=True)'
            '\nneib_inst.create_indexfile()'
            '\norder_inst.create_orderfile()'
            '\nanalysis.lateraldistribution.write_neighbortype_distr(sysinfo_inst)'
            '\nanalysis.create_leaflet_assignment_file(sysinfo_inst)'
            '\nos.remove(sys.argv[0])'.format(inputfilename, refatoms),
            file=scriptf)
        write_submitfile('submit.sh', jobfilename, mem='16G', prio=False)
        cmd = ['sbatch', '-J', jobfilename, 'submit.sh','python3', scriptfilename]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = proc.communicate()
        print(out.decode(), err.decode())

def calc_scd(systemname, temperature, jobname, *args,
    inputfilename="inputfile",
    **kwargs,):
    ''' Only calculate scd distribution and write to scd_distribution.dat '''
    complete_systemname = './{}_{}'.format(systemname, temperature)
    os.chdir(complete_systemname)
    scriptfilename = 'exec'+complete_systemname[2:]+jobname+'.py'
    jobfilename = complete_systemname[2:]+jobname
    with open(scriptfilename, 'w') as scriptf:
        print(
            'import os, sys'
            '\nfrom bilana.analysis.order import Order'
            '\nOrder(inputfilename="{0}").create_scdfile()'
            '\nos.remove(sys.argv[0])'.format(inputfilename),
            file=scriptf)
        write_submitfile('submit.sh', jobfilename, mem='16G', prio=False)
        cmd = ['sbatch', '-J', jobfilename, 'submit.sh','python3', scriptfilename]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = proc.communicate()
        print(out.decode(), err.decode())

def check_and_write(systemname, temperature, jobname, lipidpart, *args,
    overwrite=True,
    startdivisor=80,
    neighborfilename="neighbor_info",
    inputfilename="inputfile",
    energyfilename="all_energies.dat",
    scdfilename="scd_distribution.dat",
    **kwargs,):
    ''' Check if all energy files exist and write table with all energies '''
    raise NotImplementedError("Not yet implemented")
    complete_systemname = './{}_{}'.format(systemname, temperature)
    os.chdir(complete_systemname)
    scriptfilename = 'exec'+complete_systemname[2:]+jobname+'.py'
    jobfilename = complete_systemname[2:]+jobname
    if overwrite:
        oflag = "--overwrite"
    else:
        oflag = ""
    with open(scriptfilename, 'w') as scriptf:
        print(\
            'import os, sys'
            '\nimport subprocess'
            '\nfrom bilana.analysis.energy import Energy'
            '\nfrom bilana.files.eofs import EofScd'
            '\nenergy_instance = Energy("{0}", overwrite="{1}", inputfilename="{2}", neighborfilename="{3}")'
            '\nenergy_instance.run_calculation(resids={4})'
            '\nif energy_instance.check_exist_xvgs():'
            '\n    energy_instance.write_energyfile()'
            '\n    eofs = EofScd("{0}", inputfilename="{2}", energyfilename="{4}", scdfilename="{5}")'
            '\n    eofs.create_eofscdfile()'
            '\nelse:'
            '\n    print("Submitting missing energy calculations")'
            '\n    os.chdir("{0}")'
            '\n    cmd = ["bilana", "energy", "-f", "{0}", "-i", "{2}", -N, "{3}"'
            '\n        "{6}", --divisor, "{7}", -J, "en" ]'
            '\n    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)'
            '\n    out, err = proc.communicate()'
            '\n    proc.wait()'
            '\n    proc.stdout.close()'
            '\n    proc.stderr.close()'
            '\n    print(out.decode(), err.decode())'
            '\n    os.chdir("{0}")'
            '\nos.remove(sys.argv[0])'.format(lipidpart, overwrite,
                inputfilename, neighborfilename, energyfilename, scdfilename, oflag, startdivisor),
            file=scriptf)
        write_submitfile('submit.sh', jobfilename, mem='16G')
        cmd = ['sbatch', '-J', jobfilename, 'submit.sh','python3', scriptfilename]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = proc.communicate()
        print(out.decode(), err.decode())

def write_eofscd(systemname, temperature, jobname, lipidpart, *args,
    inputfilename="inputfile",
    neighborfilename="neighbor_info",
    energyfilename="all_energies.dat",
    scdfilename="scd_distribution.dat",
    **kwargs,):
    ''' Check if all energy files exist and write table with all energies '''
    complete_systemname = './{}_{}'.format(systemname, temperature)
    os.chdir(complete_systemname)
    scriptfilename = 'exec'+complete_systemname[2:]+jobname+'.py'
    jobfilename = complete_systemname[2:]+jobname
    with open(scriptfilename, 'w') as scriptf:
        print(\
            'import os, sys'
            '\nfrom bilana.analysis.energy import Energy'
            '\nfrom bilana.files.eofs import EofScd'
            '\nif energy_instance.check_exist_xvgs():'
            '\n    eofs = EofScd("{1}", inputfilename="{2}", energyfilename="{5}", scdfilename="{4}")'
            '\n    eofs.create_eofscdfile(neighbor_filename="{3}")'
            '\nelse:'
            '\n    raise ValueError("There are .edr files missing.")'
            '\nos.remove(sys.argv[0])'.format(lipidpart, inputfilename, neighborfilename,  scdfilename, energyfilename),
            file=scriptf)
        write_submitfile('submit.sh', jobfilename, mem='16G', prio=True)
        cmd = ['sbatch', '-J', jobfilename, 'submit.sh','python3', scriptfilename]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = proc.communicate()
        print(out.decode(), err.decode())

def write_nofscd(systemname, temperature, jobname, *args,
    inputfilename="inputfile",
    neighborfilename="neighbor_info",
    scdfilename="scd_distribution.dat",
    **kwargs,):
    ''' Check if all energy files exist and write table with all energies '''
    complete_systemname = './{}_{}'.format(systemname, temperature)
    os.chdir(complete_systemname)
    scriptfilename = 'exec'+complete_systemname[2:]+jobname+'.py'
    jobfilename = complete_systemname[2:]+jobname
    with open(scriptfilename, 'w') as scriptf:
        print(\
            '\nimport os, sys'
            '\nfrom bilana.files.nofs import NofScd'
            '\nnofs = NofScd(inputfilename="{}")'
            '\nnofs.create_NofScd_input(scdfilename="{}", neighborfilename="{}")'
            '\nos.remove(sys.argv[0])'.format(inputfilename, scdfilename, neighborfilename),
            file=scriptf)
        write_submitfile('submit.sh', jobfilename, mem='8G', prio=True)
        cmd = ['sbatch', '-J', jobfilename, 'submit.sh','python3', scriptfilename]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = proc.communicate()
        print(out.decode(), err.decode())
