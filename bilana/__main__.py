import sys
import argparse
from . import log
from . import command_line as cmd

commandline_modules = ["initialize", "energy", "assemble_energies", "nofscd", "eofscd", "order"]

LOGGER = log.LOGGER

PARSER = argparse.ArgumentParser()

# Positional arguments
PARSER.add_argument("mode", help="Sets mode on which function should be executed. Can be on of {}".format(commandline_modules))


# Non optional parameters
PARSER.add_argument('-f', action="store", metavar='systembasename', required=True, help="Bilayer descriptor like dppc_chol20")
PARSER.add_argument('-T', action="store", metavar='temperature',    required=True, help="Simulation temperature in K", type=int)
PARSER.add_argument('-J', action="store", metavar='jobname',        required=True, help="Job basename")


# Optional files
PARSER.add_argument('-N', '--neighborfile', action="store", nargs='?', metavar='neighbor_info', required=False, default="neighbor_info",        help="Name of file that stores neighbor information")
PARSER.add_argument('-E', '--energyfile',   action="store", nargs='?', metavar='all_energies.dat',    required=False, default="all_energies.dat",     help="Name of file that stores lipid interaction data")
PARSER.add_argument('-S', '--scdfile',      action="store", nargs='?', metavar='scd_distribution.dat',       required=False, default="scd_distribution.dat", help="Name of file that stores order parameter distribution")
PARSER.add_argument('-i', '--inputfile',    action="store", nargs='?', metavar='inputfile',     required=False, default="inputfile",            help="Name of file that stores name of the analysis inputfule")

# Arguments only for energy calculations
PARSER.add_argument('-p',         action="store", metavar="complete", required=False, default="complete", help="Sets which part of the lipid should be taken into account for interaction energy calculation")
PARSER.add_argument('--divisor',  action="store", metavar="#divisor", required=False, default=80, type=int, help="Sets the number by which the system should be divided for the energy calculation")

# On/Off flags
PARSER.add_argument('--overwrite', action="store_true", required=False, help="If this flag is set all files will be overwritten")
PARSER.add_argument('--debug', action="store_true", help="Sets logger to debug mode. This will print out a lot.")
PARSER.add_argument('--dryrun', action="store_true", help="If set, jobscripts are not submitted.")

# Arbitrary flags
PARSER.add_argument('--arbitrary', nargs='*', help="Store kwargs that is not yet listed in other arguments. Input like key:val")

ARGS = PARSER.parse_args()

kwargs = {}
if ARGS.arbitrary is not None:
    for string in ARGS.arbitrary:
        key, val = string.split(':')
        kwargs[key] = val

if ARGS.debug:
    LOGGER.setLevel("DEBUG")
LOGGER.debug("Arguments of argparse: %s", ARGS)

COMMAND = {
    "initialize":cmd.initialize_system,
    "energy":cmd.submit_energycalcs,
    "assemble_energies":cmd.check_and_write,
    "nofscd":cmd.write_nofscd,
    "eofscd":cmd.write_eofscd,
    "order":cmd.calc_scd,
}

if ARGS.mode not in commandline_modules:
    LOGGER.error("Command not found. Choose one of %s", commandline_modules)
    sys.exit()

def run_cmd():
    COMMAND[ARGS.mode](
        ARGS.f, ARGS.T, ARGS.J, ARGS.p,
        neighborfilename=ARGS.neighborfile,
        scdfilename=ARGS.scdfile,
        inputfilename=ARGS.inputfile,
        energyfilename=ARGS.energyfile,
        overwrite=ARGS.overwrite,
        startdivisor=ARGS.divisor,
        dry=ARGS.dryrun,
        **kwargs,
        )

run_cmd()