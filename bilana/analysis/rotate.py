        
import argparse
import MDAnalysis as mda

parser = argparse.ArgumentParser(
        description='Rotates protein around a selected axis.')
parser.add_argument( "str_file", help="structure file [pdb, gro]")
parser.add_argument( "axis",     help="x|y|z or x,y,z where x,y,z are 0 or 1")
parser.add_argument( "deg", type=float, help="deg")
parser.add_argument( "out_file", help="output fig file [pdb, gro]")
args = parser.parse_args()
rot = { "x": [1,0,0],
        "y": [0,1,0],
        "z": [0,0,1]
}
#----------------------------------------
axis = args.axis.lower()
if axis.find(",")>-1:
    rotc = [int(v) for v in axis.split(",")]
elif axis in "x y z".split():
    rotc = rot[axis]
else:
    raise ValueError("Axis should be x, y, or z.")
u = mda.Universe(args.str_file)
s = u.select_atoms("all")
s.rotateby( args.deg, rot[axis])
s.write(args.out_file)