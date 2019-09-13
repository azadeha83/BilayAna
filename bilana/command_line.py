'''
    This module stores functions that write script files for various calculations
    They are used in __main__.py
    All functions must use the same keyword parameter names and the *args and **kwargs parameters as input
'''
import os
import subprocess
from .common import write_submitfile
from .common import get_minmaxdiv
from .systeminfo import SysInfo

def submit_energycalcs(systemname, temperature, jobname, lipidpart, *args,
    inputfilename="inputfile",
    neighborfile="neighbor_info",
    startdivisor=80,
    overwrite=False,
    cores=2,
    dry=False,
    **kwargs,):
    ''' Divide energyruns into smaller parts for faster computation and submit those runs '''
    complete_name = './{}_{}'.format(systemname, temperature)
    os.chdir(complete_name)
    mysystem = SysInfo(inputfilename)
    systemsize = mysystem.number_of_lipids
    divisor = get_minmaxdiv(startdivisor, systemsize)
    if divisor % 1 != 0:
        raise ValueError("divisor must be int")
    resids_to_calculate = mysystem.MOLRANGE
    lipids_per_part = systemsize//divisor
    print("System and temperature:", systemname, temperature)
    print("Will overwrite:", overwrite)
    print("Lipids per job:", lipids_per_part)
    for jobpart in range(divisor):
        list_of_res = resids_to_calculate[jobpart*lipids_per_part:(jobpart+1)*lipids_per_part]
        jobfile_name = str(jobpart)+'_'+jobname
        jobscript_name = 'exec_energycalc'+str(jobfile_name)+'.py'
        with open(jobscript_name, "w") as jobf:
            print(
                '\nimport os, sys'
                '\nfrom bilana.analysis.energy import Energy'
                '\nenergy_instance = Energy("{0}", overwrite={1}, inputfilename="{2}", neighborfilename="{3}")'
                '\nenergy_instance.info()'
                '\nenergy_instance.run_calculation(resids={4})'
                '\nos.remove(sys.argv[0])'.format(lipidpart, overwrite, inputfilename, neighborfile, list_of_res),
                file=jobf)
        if not dry:
            write_submitfile('submit.sh', jobfile_name, ncores=cores)
            cmd = ['sbatch', '-J', complete_name[2:]+str(jobpart)+'_'+jobname, 'submit.sh','python3', jobscript_name]
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = proc.communicate()
            print(out.decode(), err.decode())

def initialize_system(systemname, temperature, jobname, *args,
    inputfilename="inputfile",
    refatoms="name P O3",
    cores=16,
    prio=False,
    dry=False,
    mem='32G',
    **kwargs):
    ''' Creates all core files like neighbor_info, resindex_all, scd_distribution '''
    complete_systemname = './{}_{}'.format(systemname, temperature)
    os.chdir(complete_systemname)
    scriptfilename = 'exec'+complete_systemname[2:]+jobname+'.py'
    jobfilename = complete_systemname[2:]+jobname
    with open(scriptfilename, 'w') as scriptf:
        print(\
            'import os, sys, gc'
            '\nimport bilana'
            '\nfrom bilana import SysInfo'
            '\nfrom bilana import analysis'
            '\nfrom bilana.analysis.neighbors import Neighbors'
            '\nfrom bilana.analysis.order import Order'
            '\nneib_inst = Neighbors(inputfilename="{0}")'
            '\nneib_inst.info()'
            '\nneib_inst.determine_neighbors(refatoms="{1}", overwrite=True)'
            '\ngc.collect()'
            '\nneib_inst.create_indexfile()'
            '\nanalysis.lateraldistribution.write_neighbortype_distr(SysInfo(inputfilename="{0}"))'
            '\nanalysis.leaflets.create_leaflet_assignment_file(SysInfo(inputfilename="{0}"))'
            '\nanalysis.order.calc_tilt(SysInfo(inputfilename="{0}"))'
            .format(inputfilename, refatoms),
            file=scriptf)
        if not dry:
            write_submitfile('submit.sh', jobfilename, mem=mem, ncores=cores, prio=prio)
            cmd = ['sbatch', '-J', jobfilename, 'submit.sh','python3', scriptfilename]
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = proc.communicate()
            print(out.decode(), err.decode())

def calc_scd(systemname, temperature, jobname, *args,
    inputfilename="inputfile",
    dry=False,
    **kwargs,):
    ''' Only calculate scd distribution and write to scd_distribution.dat '''
    complete_systemname = './{}_{}'.format(systemname, temperature)
    os.chdir(complete_systemname)
    scriptfilename = 'exec'+complete_systemname[2:]+jobname+'.py'
    jobfilename = complete_systemname[2:]+jobname
    with open(scriptfilename, 'w') as scriptf:
        print(
            'import os, sys, gc'
            '\nfrom bilana.systeminfo import SysInfo'
            '\nfrom bilana.analysis.order import Order, calc_tilt'
            '\ncalc_tilt(SysInfo())'
            '\ngc.collect()'
            '\nOrder(inputfilename="{0}").create_orderfile()'
            '\nos.remove(sys.argv[0])'.format(inputfilename),
            file=scriptf)
        if not dry:
            write_submitfile('submit.sh', jobfilename, mem='100G', ncores=16, prio=False)
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
    dry=False,
    **kwargs,):
    ''' Check if all energy files exist and write table with all energies '''
    complete_systemname = './{}_{}'.format(systemname, temperature)
    os.chdir(complete_systemname)
    scriptfilename = 'exec'+complete_systemname[2:]+jobname+'.py'
    jobfilename = complete_systemname[2:]+jobname
    if overwrite:
        oflag = "--overwrite"
    else:
        oflag = ""
    with open(scriptfilename, 'w') as scriptf:
        print(
            'import os, sys'
            '\nimport subprocess'
            '\nfrom bilana.analysis.energy import Energy'
            '\nfrom bilana.files.eofs import EofScd'
            '\nenergy_instance = Energy("{0}", overwrite="{1}", inputfilename="{2}", neighborfilename="{3}")'
            '\nenergy_instance.info()'
            '\nif energy_instance.check_exist_xvgs(check_len=energy_instance.universe.trajectory[-1].time):'
            '\n    energy_instance.write_energyfile()'
            '\n    eofs = EofScd("{0}", inputfilename="{2}", energyfilename="{4}", scdfilename="{5}")'
            '\n    eofs.create_eofscdfile()'.format(lipidpart, overwrite,
                inputfilename, neighborfilename, energyfilename, scdfilename),
            file=scriptf)
        if not dry:
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
    dry=False,
    **kwargs,):
    ''' Write eofscd file from table containing all interaction energies '''
    complete_systemname = './{}_{}'.format(systemname, temperature)
    os.chdir(complete_systemname)
    scriptfilename = 'exec'+complete_systemname[2:]+jobname+'.py'
    jobfilename = complete_systemname[2:]+jobname
    with open(scriptfilename, 'w') as scriptf:
        print(\
            'import os, sys'
            '\nfrom bilana.analysis.energy import Energy'
            '\nfrom bilana.files.eofs import EofScd'
            '\nenergy_instance = Energy("{0}", inputfilename="{1}", neighborfilename="{2}")'
            '\nif energy_instance.check_exist_xvgs(check_len=energy_instance.t_end):'
            '\n    eofs = EofScd("{0}", inputfilename="{1}", energyfilename="{4}", scdfilename="{3}")'
            '\n    eofs.create_eofscdfile()'
            '\nelse:'
            '\n    raise ValueError("There are .edr files missing.")'
            '\nos.remove(sys.argv[0])'.format(lipidpart, inputfilename, neighborfilename,  scdfilename, energyfilename),
            file=scriptf)
        if not dry:
            write_submitfile('submit.sh', jobfilename, mem='16G', prio=True)
            cmd = ['sbatch', '-J', jobfilename, 'submit.sh','python3', scriptfilename]
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = proc.communicate()
            print(out.decode(), err.decode())

def write_nofscd(systemname, temperature, jobname, *args,
    inputfilename="inputfile",
    neighborfilename="neighbor_info",
    scdfilename="scd_distribution.dat",
    dry=False,
    **kwargs,):
    ''' Write nofscd file from files containing the order parameter distribution and the neighbor mapping of the system '''
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
        if not dry:
            write_submitfile('submit.sh', jobfilename, mem='8G', prio=True)
            cmd = ['sbatch', '-J', jobfilename, 'submit.sh','python3', scriptfilename]
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = proc.communicate()
            print(out.decode(), err.decode())

def submit_missing_energycalculation(res, part, systemname, temperature):
    jobfilename = "en{}.py".format(res)
    jobname = "{}_{}_res{}".format(systemname, temperature, res)
    with open(jobfilename, "w") as sf:
        print('import os, sys'
            '\nfrom bilana.analysis.energy import Energy'
            '\nenergy_instance = Energy("{}", overwrite=True, inputfilename="inputfile", neighborfilename="neighbor_info")'
            '\nenergy_instance.info()'
            '\nenergy_instance.run_calculation(resids=[{}])'
            '\nos.remove(sys.argv[0])'.format(part, res), file=sf)
    write_submitfile('submit.sh', jobname, mem='8G')
    cmd = ['sbatch', '-J', jobname, 'submit.sh','python3', jobfilename]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = proc.communicate()
    print(out.decode(), err.decode())