''' Contains information about lipid molecules '''
from .. import log

LOGGER = log.LOGGER


# Glycerolpart, also carbonylpart of FA is included
GLYCATM = ['GL1','GL2', ]

HEADATM = {
    'PC':['PO4', 'NC3'],
          }

TAILCARBS = {
    'DP':[['C1A', 'C2A', 'C3A', 'C4A'],                 #16:0
          ['C1B', 'C2B', 'C3B', 'C4B']],                #16:0
    'DI':[['C1A', 'D2A', 'D3A', 'C4A'],                              #14:0
          ['C1B', 'D2B', 'D3B', 'C4B']],
    'DU':[['C1A', 'D2A', 'D3A', 'C4A'],                              #14:0
          ['C1B', 'D2B', 'D3B', 'C4B']],                             #14:0                             #14:0
            }

SCD_TAIL_ATOMS_OF = {\
    'DP':[TAILCARBS['DP'][0][::2], TAILCARBS['DP'][1][::2]],
    'DI':[TAILCARBS['DI'][0][::2], TAILCARBS['DI'][1][::2]],
    'DU':[TAILCARBS['DU'][0][::2], TAILCARBS['DU'][1][::2]],
    'CHOL':[['C1', 'C2']],
    'ERG':[['C1', 'C2']],
                    }

HEAD_ATOMS_OF = {\
    'PC':HEADATM['PC']+GLYCATM,
    
    'CHOL':['ROH', 'R1', 'R2', 'R3', 'R4', 'R5'],
    #'CHOL':['ROH'],
    #'ch1m':['all'],
    #'CHIM':['all'],
    'ERG':['ROH', 'R1', 'R2', 'R3', 'R4', 'R5'],
                }

TAIL_ATOMS_OF = {\
    'DP':[TAILCARBS['DP'][0], TAILCARBS['DP'][1],],
    'DI':[TAILCARBS['DI'][0], TAILCARBS['DI'][1],],

    'CHOL':[['C1', 'C2']],
    'ERG':[['C1', 'C2']],
                }

CENTRAL_ATOM_OF = {
    'PC':'P',
    'CHOL':'ROH',
    'ERG':'ROH',
                 }

INCLUDED_TAILS = ['DP', 'DI', 'DU']
INCLUDED_HEADS = ['PC']
STEROLS        = ['CHOL', 'ERG']
PROTEINS       = ['VAL', 'GLY', 'ALA', 'ILE', 'LEU', 'CYS', 'ARG', 'HSD']
SOLVENTS       = ["W", "TIP3", "SOL", "CL", "POT", "NA"]
#SHORTESTCHAIN  = len(TAILCARBS['DM'])



def is_sterol(lipid):
    if lipid in STEROLS:
        return True
    else:
        return False
def is_protein(lipid):
    if lipid in PROTEINS:
        return True
    else:
        return False

def central_atom_of(lipid):
    ''' Return the central atom of specified lipid as string '''
    LOGGER.debug("Lipid: %s", lipid)
    if is_sterol(lipid):
        return CENTRAL_ATOM_OF[lipid]
    elif is_protein(lipid):
        return 'N'
    else:
        head = lipid[-2:]
        return CENTRAL_ATOM_OF[head]

def head_atoms_of(lipid):
    ''' Returns a list of all head atoms for specified lipid species '''
    if is_sterol(lipid):
        return HEAD_ATOMS_OF[lipid]
    else:
        head = lipid[-2:]
        return HEAD_ATOMS_OF[head]

def tail_atoms_of(lipid):
    ''' Returns a list of all tail atoms for specified lipid species '''
    if is_sterol(lipid):
        LOGGER.warning("WARNING: No hydrogens are included")
        return TAIL_ATOMS_OF[lipid]
    else:
        tail = lipid[:-2]
        return TAIL_ATOMS_OF[tail]

def tailcarbons_of(lipid):
    ''' Returns only carbon atoms of tail '''
    if is_sterol(lipid):
        return TAIL_ATOMS_OF[lipid]
    else:
        tail = lipid[:-2]
        return TAILCARBS[tail]

def scd_tail_atoms_of(lipid):
    ''' Returns a list of relevant carbons for calculation of scd '''
    if is_sterol(lipid):
        return SCD_TAIL_ATOMS_OF[lipid]
    elif is_protein(lipid):
        return [['N', 'C']]
    else:
        tail = lipid[:-2]
        return SCD_TAIL_ATOMS_OF[tail]
