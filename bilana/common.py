'''
    Definitions that are used by multiple modules are stored here
'''
import os, sys
import numpy as np
import subprocess
from multiprocessing import Pool

def find_executable(executable, path=None):
    """
    FROM:
    # https://gist.github.com/4368898
    # Public domain code by anatoly techtonik <techtonik@gmail.com>
    # AKA Linux `which` and Windows `where`
    
    Find if 'executable' can be run. Looks for it in 'path'
    (string that lists directories separated by 'os.pathsep';
    defaults to os.environ['PATH']). Checks for all executable
    extensions. Returns full path or None if no command is found.
    """
    if path is None:
        path = os.environ['PATH']
    paths = path.split(os.pathsep)
    extlist = ['']
    if os.name == 'os2':
        (base, ext) = os.path.splitext(executable)
        # executable files on OS/2 can have an arbitrary extension, but
        # .exe is automatically appended if no dot is present in the name
        if not ext:
            executable = executable + ".exe"
    elif sys.platform == 'win32':
        pathext = os.environ['PATHEXT'].lower().split(os.pathsep)
        (base, ext) = os.path.splitext(executable)
        if ext.lower() not in pathext:
            extlist = pathext
    for ext in extlist:
        execname = executable + ext
        if os.path.isfile(execname):
            return execname
        else:
            for p in paths:
                f = os.path.join(p, execname)
                if os.path.isfile(f):
                    return f
    else:
        return None


def loop_to_pool(func, inp, maxtasknum=1000):
    ''' Map inparray to func, if inp is iterable use pool.starmap '''
    with Pool(len(os.sched_getaffinity(0)), maxtasksperchild=maxtasknum) as pool:
        if isinstance(inp, tuple) or isinstance(inp, list):
            data_outputs = pool.starmap(func, inp)
        else:
            data_outputs = pool.map(func, inp)
    pool.close()
    pool.join()
    return data_outputs

def exec_gromacs(cmd, inp_str=None):
    '''
        Execute Gromacs commands.
        cmd variable stores the actual terminal command (can be anything,
        but here is specifically for gromacs)
        cmd must be a list with items being the strings (no whitespaces) of the command
            "['gmx cmd', '-f', 'tprfile', '-e', 'en.edr']"
        if a command needs user input (STDIN, like in gmx make_ndx) inp_str can be parsed
        as string
    '''
    if inp_str is None:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        out, err = proc.communicate()
    else:
        try:
            inp_str = inp_str.encode()
        except AttributeError:
            pass
        proc = subprocess.Popen(cmd, stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = proc.communicate(inp_str)
    proc.wait()
    proc.stdout.close()
    proc.stderr.close()
    #if proc.returncode != 0:
    if proc.returncode == "noval":
        print("CODE", proc.returncode)
        if proc.returncode == 132:
            raise ChildProcessError("Core was dumped. This is probably due to an incompatible gromacs version")
        try:
            err = err.decode()
            print(out)
            print(err)
        except UnboundLocalError:
            pass
        raise ChildProcessError(
            'Failed to execute command "{}" and input "{}"'\
            .format(' '.join(cmd), inp_str))
    return out.decode(), err.decode()

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

def write_submitfile(submitout, jobname, ncores=2, mem='4G', prio=False, queue=None):
    username = "aalaviza"
    hostname = "bagheera"
    if hostname == "bagheera":
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
                  '\n#SBATCH --time=48:00:00'
                  '\n#SBATCH --ntasks=1'
                  '\n#SBATCH --nodes=1'
                  '\n#SBATCH --cpus-per-task={ncores}'
                  '\n#SBATCH --mem={mem}'
                  '\n#SBATCH --exclude=kaa-[72,73,76,77,78,82,12,13]'
                  '\n#SBATCH --constraint="avx|avx2|fma4"'
                  '\nexport OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK'\
                  '\nsrun $@'.format(**locals()), file=sfile)

    elif hostname == 'r06m01':
        if queue is None:
            queue='hims,q0heuer,normal'
        with open(submitout, "w") as sfile:
            print('#!/bin/bash'
                  '\n#SBATCH -A q0heuer'
                  '\n#SBATCH -p {queue}' 
                  '\n#SBATCH --output={jobname}.out'
                  '\n#SBATCH --mail-type=fail'
                  '\n#SBATCH --mail-user={username}@wwu.de'
                  '\n#SBATCH --time=48:00:00'
                  '\n#SBATCH --ntasks=1'
                  '\n#SBATCH --nodes=1'
                  '\n#SBATCH --cpus-per-task={ncores}'
                  '\n#SBATCH --mem={mem}'
                  '\nexport OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK'\
                  '\nsrun $@'.format(**locals()), file=sfile)


def make_movie_from_images(moviename, picture_name, framerate=40):
    ''' Create an .mpg movie from picture files
        <picture_name> must be a string describing the name of the pictures
        example:  %02d_voro.png  for pictures like 01_voro.png 02_voro.png ...

    '''
    cmd = ["ffmpeg", "-framerate", framerate,
    "-i", picture_name,
    "-c:v", "mpeg2video",
    "-pix_fmt", "yuv420p",
    moviename]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = proc.communicate()
    proc.wait()
    proc.stdout.close()
    proc.stderr.close()
    print(out.decode())
    print(err.decode())


def rotate_2d(vec, angle_rad):
    mat = np.matrix([[np.cos(angle_rad), np.sin(angle_rad) ], [-1*np.sin(angle_rad), np.cos(angle_rad) ]])
    nvec = np.array( mat.dot(vec) )[0]
    return nvec

def angle_clockwise(A, B):
    inner_ang = np.dot(A, B) / np.linalg.norm(A) * np.linalg.norm(B)
    det = np.linalg.det(np.matrix([A, B]))
    if det < 0:
        return inner_ang
    else:
        return 2*np.pi - inner_ang

GMXNAME = find_executable("gmx")

