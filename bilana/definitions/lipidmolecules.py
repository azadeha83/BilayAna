''' Contains information about lipid molecules '''

from ..systeminfo import SysInfo
from .. import log

LOGGER = log.LOGGER
try:
    with open('force_field','r') as f:
        ls = f.readlines()
        if ls[0].strip() == 'all_atom':
            ff = 'all_atom'
        else:
            ff = 'coarse'
except FileNotFoundError:
    ff = 'all_atom'

# Glycerolpart, also carbonylpart of FA is included
if ff == 'all_atom':
    
    GLYCATM = ['C1', 'O11', 'C2', 'O21', 'C21', 'O22', 'C3', 'O31', 'C31', 'O32', 'HA', 'HB', 'HY', 'HX', 'HS', ]
else:
    GLYCATM = ['GL1','GL2', ]

if ff == 'all_atom':

    HEADATM = {
        'PC':['P', 'O11', 'O12', 'O13', 'O14', 'N', 'C11', 'C12', 'C13', 'C14', 'C15',
            'H11A', 'H11B', 'H12A', 'H12B', 'H13A', 'H13B', 'H13C', 'H14A', 'H14B', 'H14C', 'H15A', 'H15B', 'H15C'],
        'PA':['P', 'O12', 'O13', 'O14', 'H12'],
        'PE':['P', 'O12', 'O13', 'O14', 'C11', 'C12', 'N', 'HN1', 'HN2', 'HN3', 'H11A', 'H11B', 'H12A', 'H12B',],
        'PI':['P', 'O12', 'O13', 'O14', 'C11', 'C12', 'O2', 'C13', 'O3', 'C14', 'O4', 'C15', 'O5', 'C16', 'O6',
            'H1', 'H2', 'HO2', 'H3', 'HO3', 'H4', 'HO4', 'H5', 'HO5', 'H6', 'HO6'],
        'PS':['P', 'O12', 'O13', 'O14', 'C11', 'C12', 'C13', 'O13A', 'O13B', 'N', 'H11A', 'H11B', 'H12A', 'HN1', 'HN2', 'HN3'],
        }
else:
    HEADATM = {
        'PC':['PO4', 'NC3'],
        'PS':['PO4', 'CNO'],
            }

if ff == 'all_atom':

    TAILCARBS = {
    'DP':[['C22', 'C23', 'C24', 'C25', 'C26', 'C27', 'C28', 'C29', 'C210', 'C211', 'C212', 'C213', 'C214', 'C215', 'C216'],                  #16:0
          ['C32', 'C33', 'C34', 'C35', 'C36', 'C37', 'C38', 'C39', 'C310', 'C311', 'C312', 'C313', 'C314', 'C315', 'C316']],                 #16:0
    'DM':[['C22', 'C23', 'C24', 'C25', 'C26', 'C27', 'C28', 'C29', 'C210', 'C211', 'C212', 'C213', 'C214',],                                 #14:0
          ['C32', 'C33', 'C34', 'C35', 'C36', 'C37', 'C38', 'C39', 'C310', 'C311', 'C312', 'C313', 'C314',]],                                #14:0
    'DS':[['C22', 'C23', 'C24', 'C25', 'C26', 'C27', 'C28', 'C29', 'C210', 'C211', 'C212', 'C213', 'C214', 'C215', 'C216', 'C217', 'C218'],  #18:0
          ['C32', 'C33', 'C34', 'C35', 'C36', 'C37', 'C38', 'C39', 'C310', 'C311', 'C312', 'C313', 'C314', 'C315', 'C316', 'C317', 'C318']], #18:0
    'DU':[['C22', 'C23', 'C24', 'C25', 'C26', 'C27', 'C28', 'C29', 'C210', 'C211', 'C212', 'C213', 'C214', 'C215', 'C216', 'C217', 'C218'],  #18:2
          ['C32', 'C33', 'C34', 'C35', 'C36', 'C37', 'C38', 'C39', 'C310', 'C311', 'C312', 'C313', 'C314', 'C315', 'C316', 'C317', 'C318']], #18:2
    'DLI':[['C22', 'C23', 'C24', 'C25', 'C26', 'C27', 'C28', 'C29', 'C210', 'C211', 'C212', 'C213', 'C214', 'C215', 'C216', 'C217', 'C218'],  #18:2
          ['C32', 'C33', 'C34', 'C35', 'C36', 'C37', 'C38', 'C39', 'C310', 'C311', 'C312', 'C313', 'C314', 'C315', 'C316', 'C317', 'C318']], #18:2
    'DY':[['C22', 'C23', 'C24', 'C25', 'C26', 'C27', 'C28', 'C29', 'C210', 'C211', 'C212', 'C213', 'C214', 'C215', 'C216'],                  #16:1
          ['C32', 'C33', 'C34', 'C35', 'C36', 'C37', 'C38', 'C39', 'C310', 'C311', 'C312', 'C313', 'C314', 'C315', 'C316']],                 #16:1
    'PL':[['C22', 'C23', 'C24', 'C25', 'C26', 'C27', 'C28', 'C29', 'C210', 'C211', 'C212', 'C213', 'C214', 'C215', 'C216', 'C217', 'C218'],  #18:2
          ['C32', 'C33', 'C34', 'C35', 'C36', 'C37', 'C38', 'C39', 'C310', 'C311', 'C312', 'C313', 'C314', 'C315', 'C316']],                 #16:0
    'PO':[['C22', 'C23', 'C24', 'C25', 'C26', 'C27', 'C28', 'C29', 'C210', 'C211', 'C212', 'C213', 'C214', 'C215', 'C216', 'C217', 'C218'],  #18:1
          ['C32', 'C33', 'C34', 'C35', 'C36', 'C37', 'C38', 'C39', 'C310', 'C311', 'C312', 'C313', 'C314', 'C315', 'C316']],                 #16:0
    'DO':[['C22', 'C23', 'C24', 'C25', 'C26', 'C27', 'C28', 'C29', 'C210', 'C211', 'C212', 'C213', 'C214', 'C215', 'C216', 'C217', 'C218'],  #18:1
          ['C32', 'C33', 'C34', 'C35', 'C36', 'C37', 'C38', 'C39', 'C310', 'C311', 'C312', 'C313', 'C314', 'C315', 'C316', 'C317', 'C318']], #18:1
    }
else:
    TAILCARBS = {
        'DP':[['C1A', 'C2A', 'C3A', 'C4A'],                 #16:0
            ['C1B', 'C2B', 'C3B', 'C4B']],                #16:0
        'DI':[['C1A', 'D2A', 'D3A', 'C4A'],                              #14:0
            ['C1B', 'D2B', 'D3B', 'C4B']],
        'DU':[['C1A', 'D2A', 'D3A', 'C4A'],                              #14:0
            ['C1B', 'D2B', 'D3B', 'C4B']],
        'DO':[['C1A', 'D2A', 'C3A', 'C4A'],                              #14:0
            ['C1B', 'D2B', 'C3B', 'C4B']], 
        'PO':[['C1A', 'D2A', 'C3A', 'C4A'],                              #14:0
            ['C1B', 'C2B', 'C3B', 'C4B']],                                            #14:0                             #14:0
                }

if ff == 'all_atom':
    
    TAILHYDR = {
        'DP':[['H2R', 'H2S','H3R', 'H3S', 'H4R', 'H4S', 'H5R', 'H5S', 'H6R', 'H6S', 'H7R', 'H7S', 'H8R', 'H8S', 'H9R', 'H9S', 'H10R', 'H10S', 'H11R', 'H11S', 'H12R', 'H12S', 'H13R', 'H13S', 'H14R', 'H14S', 'H15R', 'H15S','H16R', 'H16S', 'H16T'],
            ['H2X', 'H2Y','H3X', 'H3Y', 'H4X', 'H4Y', 'H5X', 'H5Y', 'H6X', 'H6Y', 'H7X', 'H7Y', 'H8X', 'H8Y', 'H9X', 'H9Y', 'H10X', 'H10Y', 'H11X', 'H11Y', 'H12X', 'H12Y', 'H13X', 'H13Y', 'H14X', 'H14Y', 'H15X', 'H15Y','H16X', 'H16Y', 'H16Z']],
        'DM':[['H2R', 'H2S','H3R', 'H3S', 'H4R', 'H4S', 'H5R', 'H5S', 'H6R', 'H6S','H7R','H7S','H8R',
                    'H8S', 'H9R', 'H9S', 'H10R', 'H10S','H11R', 'H11S', 'H12R', 'H12S','H13R', 'H13S',
                    'H14R', 'H14S', 'H14T'],
                ['H2X', 'H2Y','H3X', 'H3Y','H4X', 'H4Y','H5X', 'H5Y', 'H6X', 'H6Y','H7X', 'H7Y','H8X',
                    'H8Y','H9X', 'H9Y','H10X', 'H10Y','H11X', 'H11Y','H12X', 'H12Y', 'H13X', 'H13Y','H14X',
                    'H14Y', 'H14Z']],
        'DS':[['H2R', 'H2S','H3R', 'H3S', 'H4R', 'H4S', 'H5R', 'H5S', 'H6R', 'H6S','H7R','H7S','H8R',
                    'H8S', 'H9R', 'H9S', 'H10R', 'H10S','H11R', 'H11S', 'H12R', 'H12S','H13R', 'H13S',
                    'H14R', 'H14S', 'H15R', 'H15S','H16R', 'H16S', 'H17R', 'H17S', 'H118R', 'H18S', 'H18T'],
                ['H2X', 'H2Y','H3X', 'H3Y','H4X', 'H4Y','H5X', 'H5Y', 'H6X', 'H6Y','H7X', 'H7Y','H8X',
                    'H8Y','H9X', 'H9Y','H10X', 'H10Y','H11X', 'H11Y','H12X', 'H12Y', 'H13X', 'H13Y','H14X',
                    'H14Y','H15X', 'H15Y','H16X', 'H16Y','H17X', 'H17Y', 'H18X', 'H18Y', 'H18Z']],
        'DU':[['H2R', 'H2S', 'H3R', 'H3S', 'H4R', 'H4S', 'H5R', 'H5S', 'H6R', 'H6S', 'H7R', 'H7S','H8R', 'H8S', 'H9R', 'H10R', 'H11R', 'H11S', 'H12R', 'H13R', 'H14R', 'H14S', 'H15R', 'H15S', 'H16R', 'H16S', 'H17R', 'H17S', 'H18R', 'H18S', 'H18T'],
            ['H2X', 'H2Y', 'H3X', 'H3Y', 'H4X', 'H4Y', 'H5X', 'H5Y', 'H6X', 'H6Y', 'H7X', 'H7Y','H8X', 'H8Y', 'H9X', 'H10X', 'H11X', 'H11Y', 'H12X', 'H13X', 'H14X', 'H14Y', 'H15X', 'H15Y', 'H16X', 'H16Y', 'H17X', 'H17Y', 'H18X', 'H18Y', 'H18Z']],
        
        'DLI':[['H2R', 'H2S', 'H3R', 'H3S', 'H4R', 'H4S', 'H5R', 'H5S', 'H6R', 'H6S', 'H7R', 'H7S','H8R', 'H8S', 'H9R', 'H10R', 'H11R', 'H11S', 'H12R', 'H13R', 'H14R', 'H14S', 'H15R', 'H15S', 'H16R', 'H16S', 'H17R', 'H17S', 'H18R', 'H18S', 'H18T'],
            ['H2X', 'H2Y', 'H3X', 'H3Y', 'H4X', 'H4Y', 'H5X', 'H5Y', 'H6X', 'H6Y', 'H7X', 'H7Y','H8X', 'H8Y', 'H9X', 'H10X', 'H11X', 'H11Y', 'H12X', 'H13X', 'H14X', 'H14Y', 'H15X', 'H15Y', 'H16X', 'H16Y', 'H17X', 'H17Y', 'H18X', 'H18Y', 'H18Z']],
        
        'DY':[['H2R', 'H2S','H3R', 'H3S', 'H4R', 'H4S', 'H5R', 'H5S', 'H6R', 'H6S','H7R','H7S','H8R',
                    'H8S', 'H9R', 'H10R', 'H11R', 'H11S', 'H12R', 'H12S','H13R', 'H13S',
                    'H14R', 'H14S', 'H15R', 'H15S','H16R', 'H16S', 'H16T'],
                ['H2X', 'H2Y','H3X', 'H3Y','H4X', 'H4Y','H5X', 'H5Y', 'H6X', 'H6Y','H7X', 'H7Y','H8X',
                    'H8Y','H9X', 'H10X', 'H11X', 'H11Y','H12X', 'H12Y', 'H13X', 'H13Y','H14X',
                    'H14Y','H15X', 'H15Y','H16X', 'H16Y','H16Z']],
        'PL':[['H2R', 'H2S','H3R', 'H3S', 'H4R', 'H4S', 'H5R', 'H5S', 'H6R', 'H6S','H7R','H7S','H8R',
                    'H8S', 'H9R', 'H10R', 'H11R', 'H11S', 'H12R', 'H13R',
                    'H14R', 'H14S', 'H15R', 'H15S','H16R', 'H16S', 'H17X', 'H17Y', 'H18X', 'H18Y', 'H18Z'],
                ['H2X', 'H2Y','H3X', 'H3Y','H4X', 'H4Y','H5X', 'H5Y', 'H6X', 'H6Y','H7X', 'H7Y','H8X',
                    'H8Y','H9X', 'H9Y','H10X', 'H10Y','H11X', 'H11Y','H12X', 'H12Y', 'H13X', 'H13Y','H14X',
                    'H14Y','H15X', 'H15Y','H16X', 'H16Y','H16Z']],
        'PO':[['H2R', 'H2S','H3R', 'H3S', 'H4R', 'H4S', 'H5R', 'H5S', 'H6R', 'H6S','H7R','H7S','H8R',
                    'H8S', 'H91', 'H101', 'H11R', 'H11S', 'H12R', 'H12S','H13R', 'H13S',
                    'H14R', 'H14S', 'H15R', 'H15S','H16R', 'H16S', 'H17R', 'H17S', 'H118R', 'H18S', 'H18T'],
                ['H2X', 'H2Y','H3X', 'H3Y','H4X', 'H4Y','H5X', 'H5Y', 'H6X', 'H6Y','H7X', 'H7Y','H8X',
                    'H8Y','H9X', 'H9Y','H10X', 'H10Y', 'H11X', 'H11Y','H12X', 'H12Y', 'H13X', 'H13Y','H14X',
                    'H14Y','H15X', 'H15Y','H16X', 'H16Y','H16Z']],
        'DO':[['H2R', 'H2S','H3R', 'H3S', 'H4R', 'H4S', 'H5R', 'H5S', 'H6R', 'H6S','H7R','H7S','H8R',
                    'H8S', 'H9R', 'H10R', 'H11R', 'H11S', 'H12R', 'H12S','H13R', 'H13S',
                    'H14R', 'H14S', 'H15R', 'H15S','H16R', 'H16S', 'H17R', 'H17S', 'H118R', 'H18S', 'H18T'],
                ['H2X', 'H2Y','H3X', 'H3Y','H4X', 'H4Y','H5X', 'H5Y', 'H6X', 'H6Y','H7X', 'H7Y','H8X',
                    'H8Y','H9X', 'H10X', 'H11X', 'H11Y','H12X', 'H12Y', 'H13X', 'H13Y','H14X',
                    'H14Y','H15X', 'H15Y','H16X', 'H16Y','H17X', 'H17Y', 'H18X', 'H18Y', 'H18Z']],
        }

if ff == 'all_atom':

    SCD_TAIL_ATOMS_OF = {\
    'DP':[TAILCARBS['DP'][0][::2], TAILCARBS['DP'][1][::2]],
    'DU':[['C22', 'C24', 'C26', 'C28', 'C211', 'C214', 'C216', 'C218'], # Double bonds between 9-10, 12-12
            ['C32', 'C34', 'C36', 'C38', 'C311', 'C314', 'C316', 'C318']],
    'DLI':[['C22', 'C24', 'C26', 'C28', 'C211', 'C214', 'C216', 'C218'], # Double bonds between 9-10, 12-12
            ['C32', 'C34', 'C36', 'C38', 'C311', 'C314', 'C316', 'C318']],
    'DM':[TAILCARBS['DM'][0][::2], TAILCARBS['DM'][1][::2]],
    'DS':[TAILCARBS['DS'][0][::2], TAILCARBS['DS'][1][::2]],
    'DY':[['C22', 'C24', 'C26', 'C28', 'C211', 'C213', 'C215', ],   # Double bonds between 9-10
            ['C32', 'C34', 'C36', 'C38', 'C311', 'C313', 'C315',]],
    'PL':[['C22', 'C24', 'C26', 'C28', 'C211', 'C214', 'C216', 'C218'], TAILCARBS['PL'][1][::2]],
    'PO':[['C22', 'C24', 'C26', 'C28', 'C211', 'C213', 'C215', 'C217'], TAILCARBS['PO'][1][::2]],
    'DO':[['C22', 'C24', 'C26', 'C28', 'C211', 'C213', 'C215', 'C217'],  # Double bonds between 9-10
            ['C32', 'C34', 'C36', 'C38', 'C311', 'C313', 'C315', 'C317']],
    'CHL1':[['C3', 'C17']],
    #'CHL1':[['C20', 'C25']],
    'ch1m':[['C20', 'C12']],
    'CHIM':[['C20', 'C12']],
    'CH0M':[['C20', 'C12']],
    'ERG':[['C3', 'C17']],
    #'ERG':[['C20', 'C25']],
                    }
else:
    SCD_TAIL_ATOMS_OF = {\
        'DP':[TAILCARBS['DP'][0][::1], TAILCARBS['DP'][1][::1]],
        'DI':[TAILCARBS['DI'][0][::1], TAILCARBS['DI'][1][::1]],
        'DU':[TAILCARBS['DU'][0][::1], TAILCARBS['DU'][1][::1]],
        'DO':[TAILCARBS['DO'][0][::1], TAILCARBS['DO'][1][::1]],
        'PO':[TAILCARBS['PO'][0][::1], TAILCARBS['PO'][1][::1]],
        'CHOL':[['R1', 'R5']],
        'ERGO':[['R1', 'R5']],
                        }

if ff == 'all_atom':

    HEAD_ATOMS_OF = {\
        'PC':HEADATM['PC']+GLYCATM,
        'PE':HEADATM['PE']+GLYCATM,
        'PS':HEADATM['PS']+GLYCATM,
        'PI':HEADATM['PI']+GLYCATM,
        'PA':HEADATM['PA']+GLYCATM,
        'CHL1':['O3', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C10', 'C11', 'C12', 'C13', 'C14', 'C15', 'C16', 'C17', 'C18', 'C19'],
        #'ch1m':['all'],
        #'CHIM':['all'],
        'ERG':['O3', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C10', 'C11', 'C12', 'C13', 'C14', 'C15', 'C16', 'C17', 'C18', 'C19'],
        'CHIM':['C20', "C1", "C2", "N1", "N2"],
        'CH0M':['C20', "C1", "C2", "N1", "N2"],                }

else:
    HEAD_ATOMS_OF = {\
        'PC':HEADATM['PC']+GLYCATM,
        'PS':HEADATM['PS']+GLYCATM,
        
        'CHOL':['ROH', 'R1', 'R2', 'R3', 'R4', 'R5'],
        #'CHOL':['ROH'],
        #'ch1m':['all'],
        #'CHIM':['all'],
        'ERGO':['ROH', 'R1', 'R2', 'R3', 'R4', 'R5'],
                    }

if ff == 'all_atom':

    TAIL_ATOMS_OF = {\
        'DP':[TAILCARBS['DP'][0], TAILHYDR['DP'][0], TAILCARBS['DP'][1], TAILHYDR['DP'][1],],
        'DM':[TAILCARBS['DM'][0], TAILHYDR['DM'][0], TAILCARBS['DM'][1], TAILHYDR['DM'][1],],
        'DS':[TAILCARBS['DS'][0], TAILHYDR['DS'][0], TAILCARBS['DS'][1], TAILHYDR['DS'][1],],
        'DO':[TAILCARBS['DO'][0], TAILHYDR['DO'][0], TAILCARBS['DO'][1], TAILHYDR['DO'][1],],
        'DY':[TAILCARBS['DY'][0], TAILHYDR['DY'][0], TAILCARBS['DY'][1], TAILHYDR['DY'][1],],
        'DU':[TAILCARBS['DU'][0], TAILHYDR['DU'][0], TAILCARBS['DU'][1], TAILHYDR['DU'][1],],
        'DLI':[TAILCARBS['DLI'][0], TAILHYDR['DLI'][0], TAILCARBS['DLI'][1], TAILHYDR['DLI'][1],],
        'PO':[TAILCARBS['PO'][0], TAILHYDR['PO'][0], TAILCARBS['PO'][1], TAILHYDR['PO'][1],],
        'PL':[TAILCARBS['PL'][0], TAILHYDR['PL'][0], TAILCARBS['PL'][1], TAILHYDR['PL'][1],],
        'CHL1':[['C13', 'C14', 'C15', 'C16', 'C17', 'C20', 'C22', 'C23', 'C24', 'C25']],
        'ERG':[['C13', 'C14', 'C15', 'C16', 'C17', 'C20', 'C22', 'C23', 'C24', 'C25', 'C28']],
        'CHIM':[['C12', 'C22', 'C23', 'C24', 'C25', 'C26', 'C27', 'C28', 'C29', 'C30']],
        'CH0M':[['C12', 'C22', 'C23', 'C24', 'C25', 'C26', 'C27', 'C28', 'C29', 'C30']],
                    }
else:
    TAIL_ATOMS_OF = {\
        'DP':[TAILCARBS['DP'][0], TAILCARBS['DP'][1],],
        'DI':[TAILCARBS['DI'][0], TAILCARBS['DI'][1],],
        'DO':[TAILCARBS['DO'][0], TAILCARBS['DO'][1],],
        'PO':[TAILCARBS['PO'][0], TAILCARBS['PO'][1],],

        'CHOL':[['C1', 'C2']],
        'ERGO':[['C1', 'C2']],
                }

if ff == 'all_atom':

    CENTRAL_ATOM_OF = {
        'PC':'P',
        'PE':'P',
        'PS':'P',
        'PI':'P',
        'PA':'P',
        'CHL1':'O3',
        'CHIM':'C20',
        'CH0M':'C20',
        'ch1m':'C20',
        'ERG':'O3',
        'WSC1':'N',
                    }
else:
    CENTRAL_ATOM_OF = {
        'PC':'P',
        'CHOL':'ROH',
        'ERGO':'ROH',
                    }

if ff == 'all_atom':

    INCLUDED_TAILS = ['DP', 'DM', 'DS', 'DO', 'DY', 'DU', 'PO', 'PL']
    INCLUDED_HEADS = ['PC', 'PE', 'PS', 'PI', 'PA']
    STEROLS        = ['CHL1', 'CHIM', 'CH0M', 'ERG', 'ch1m']
    PROTEINS       = ['VAL', 'GLY', 'ALA', 'ILE', 'LEU', 'CYS', 'ARG', 'HSD', 'ASN','HIS']
    SOLVENTS       = ["TIP3", "SOL", "CL", "POT", "NA"]
    SHORTESTCHAIN  = len(TAILCARBS['DM'])
else:
    INCLUDED_TAILS = ['DP', 'DI', 'DU', 'DO', 'PO']
    INCLUDED_HEADS = ['PC']
    STEROLS        = ['CHOL', 'ERGO',]
    PROTEINS       = ['VAL', 'GLY', 'ALA', 'ILE', 'LEU', 'CYS', 'ARG', 'HSD', 'ASN','HIS']
    SOLVENTS       = ["W", "TIP3", "SOL", "CL", "POT", "NA"]
    SHORTESTCHAIN  = len(TAILCARBS['DP'])



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
