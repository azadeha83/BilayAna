import numpy as np
import re


class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value


def calculate_geometriccenter(self,coordinateinput):
    geocenter=[0,0,0]
    for atomcoords in coordinateinput:
        for dimension in atomcoords:
            geocenter[atomcoords.index(dimension)]+=dimension
    geocenter=[x/len(coordinateinput) for x in geocenter]
    return np.array(geocenter)

class PDB_format():
    #                    Record   Serial   atomname    resn    resid        x       y           z    Rest
    pdb_string_format = '{: <4}  {: >5} {: >2}{: <2}{}{: >4}{}{: >4}{}   {: >8.3f}{: >8.3f}{: >8.3f}{}'#{: >6.2f}{: >6.2f}'
    pdb_pattern = r"""        
            ^
            ATOM              (?# Record Type)
    \s+
            (\d+)             (?# Serial numbers  Grp 1)
    \s*
            (\w+?)([\d,\w,\']+)       (?# Atom names;     Grp 2+3)
    \s+
            ([\w,\d]+)        (?# Residue name    Grp 4)
    \s+
            (\d+)             (?# Resid           Grp 5)
    \s+
            (-?\d+\.\d+)      (?# X               Grp 6-8)
    \s*
            (-?\d+\.\d+)      (?# Y)
    \s*
            (-?\d+\.\d+)      (?# Z)
            (.*)              (?# Remainder)              
                """
    pdb_pattern = ''.join(pdb_pattern.split()) # Remove all whitespaces in above defined pattern
    regexp = re.compile(pdb_pattern)


class GRO_format():
    ''' Provides regexp for .gro structure file format
    Entries are ordered as follows
    #grp     Name
    1       resid
    2       resname
    3       atom name
    4       atom number
    5-7     x y z
    8       velocities (rest)
    '''
    #                          resid resnm atmnm atmnr  x           y      z         vx      vy        vz
    gro_string_format = '{: >5}{: <5}{: >5}{: >5}{: >8.3f}{: >8.3f}{: >8.3f}{: >8.4f}{: >8.4f}{: >8.4f} '
    gro_pattern = r"""
        ([\s,\d]{5})           (?# Resid       Grp 1)
        ([\s,\w]{5})           (?# Resname     Grp 2)
        ([\w,\s,\d,\']{5})           (?# Atom name   Grp 3)
        ([\s,\d]{5})           (?# Atom number Grp 4)
    \s*
            (-?\d+\.\d+)      (?# X               Grp 5-7)
    \s*
            (-?\d+\.\d+)      (?# Y)
    \s*
            (-?\d+\.\d+)      (?# Z)
    \s*
            (.*)            (?# The rest velocities Grp 8)
    """
    gro_pattern = ''.join(gro_pattern.split())
    regexp = re.compile(gro_pattern)
    gro_box_pattern = r'\s*(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)$'
    regexp_box = re.compile(gro_box_pattern)
