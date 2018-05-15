''' Contains information about lipid molecules '''

# Head
headheteroat = ['C11', 'C12', 'C13', 'C14', 'C15', 'N', 'P', 'O11', 'O12', 'O13', 'O14', 'C1', 'C2', 'O21', 'C21', 'O22', 'C3', 'O31', 'O32', 'C31']
headhydrog = ['H13A', 'H13B', 'H13C', 'H14A', 'H14B', 'H14C', 'H15A', 'H15B', 'H15C', 'H12A', 'H12B', 'H11A', 'H11B', 'HA', 'HB', 'HS', 'HX', 'HY']
# Tail
tailcarbons_of = {
    'DPPC': [['C22', 'C23', 'C24', 'C25', 'C26', 'C27', 'C28','C29', 'C210','C211', 'C212','C213', 'C214', 'C215', 'C216'],
             ['C32', 'C33', 'C34', 'C35', 'C36', 'C37', 'C38','C39', 'C310','C311', 'C312','C313', 'C314', 'C315', 'C316']],
    'DUPC': [['C22', 'C23', 'C24', 'C25', 'C26', 'C27', 'C28','C29', 'C210','C211', 'C212','C213', 'C214', 'C215', 'C216', 'C217','C218'],
             ['C32', 'C33', 'C34', 'C35', 'C36', 'C37', 'C38','C39', 'C310','C311', 'C312','C313', 'C314', 'C315', 'C316', 'C317','C318']],
                  }

tailhydr_of = {
    'DPPC': [['H2R', 'H2S','H3R', 'H3S', 'H4R', 'H4S', 'H5R', 'H5S', 'H6R', 'H6S','H7R','H7S','H8R',
                'H8S', 'H9R', 'H9S', 'H10R', 'H10S','H11R', 'H11S', 'H12R', 'H12S','H13R', 'H13S',
                'H14R', 'H14S', 'H15R', 'H15S','H16R', 'H16S', 'H16T'],
             ['H2X', 'H2Y','H3X', 'H3Y','H4X', 'H4Y','H5X', 'H5Y', 'H6X', 'H6Y','H7X', 'H7Y','H8X',
                'H8Y','H9X', 'H9Y','H10X', 'H10Y','H11X', 'H11Y','H12X', 'H12Y', 'H13X', 'H13Y','H14X',
                'H14Y','H15X', 'H15Y','H16X', 'H16Y','H16Z']],
    'DUPC': [['H2R', 'H2S','H3R', 'H3S', 'H4R', 'H4S', 'H5R', 'H5S', 'H6R', 'H6S','H7R','H7S','H8R',
                'H8S', 'H9R', 'H9S', 'H10R', 'H10S','H11R', 'H11S', 'H12R', 'H12S', 'H13R', 'H13S',
                'H14R', 'H14S', 'H15R', 'H15S','H16R', 'H16S', 'H17X', 'H17Y', 'H18X', 'H18Y', 'H18Z'],
             ['H2X', 'H2Y','H3X', 'H3Y','H4X', 'H4Y','H5X', 'H5Y', 'H6X', 'H6Y','H7X', 'H7Y','H8X',
                'H8Y', 'H9X', 'H9Y', 'H10X', 'H10Y', 'H11X', 'H11Y', 'H12X', 'H12Y', 'H13X', 'H13Y',
                'H14X', 'H14Y', 'H15X', 'H15Y', 'H16X', 'H16Y', 'H17X', 'H17Y', 'H18X', 'H18Y', 'H18Z']],
                }
# DUPC
#tailcarbonsdupc = [tailcarbonsdppc[0]+['C217','C218'], tailcarbonsdppc[1]+['C317','C318']]
#tailhydrdupc = [tailhydrdppc[0].copy()+['H17X','H17Y','H18X','H18Y','H18Z'], tailhydrdppc[1].copy()+['H17X','H17Y','H18X','H18Y','H18Z']]
#tailhydrdupc[0].remove('H16T')
#tailhydrdupc[1].remove('H16Z')
#
scd_tail_atoms_of = {\
    'DPPC':[tailcarbons_of['DPPC'][0][::2], tailcarbons_of['DPPC'][1][::2]],
    'DUPC':[['C22', 'C24', 'C26', 'C28', 'C211', 'C214', 'C216', 'C218'], # Double bonds between 9-10, 12-13
            ['C32', 'C34', 'C36', 'C38', 'C311', 'C314', 'C316', 'C318']],
    'CHL1':[['C3', 'C17']],
    'CHIM':[['C20', 'C12']],
                    }

head_atoms_of = {\
    'DPPC':headheteroat+headhydrog,
    'DUPC':headheteroat+headhydrog,
    'CHL1':['all'],
    'CHIM':['all'],
                }

tail_atoms_of = {\
    'DPPC':[tailcarbons_of['DPPC'][0], tailhydr_of['DPPC'][0], tailcarbons_of['DPPC'][1], tailhydr_of['DPPC'][1]],
    'DUPC':[tailcarbons_of['DUPC'][0], tailhydr_of['DUPC'][0], tailcarbons_of['DUPC'][1], tailhydr_of['DUPC'][1]],
    'CHL1':['all'],
    'CHIM':['all']
                }

central_atom_of = {
    'DPPC':'P',
    'CHL1':'O3',
    'CHIM':'C20',
    'DUPC':'P',
    'DOPC':'P',
                }

described_molecules = ['DPPC', 'DUPC', 'CHL1', 'CHIM']
sterols = ['CHL1', 'CHIM']
shortestchain = len(tailcarbons_of['DPPC'])


