''' Contains information about lipid molecules '''

##### For energy groups
headhydr=['H13A', 'H13B', 'H13C', 'H14A', 'H14B', 'H14C','H15A', 'H15B', 'H15C', 'H12A', 'H12B','H11A', 'H11B','HA', 'HB', 'HS', 'HX', 'HY']
tailcarbonsdppc=[['C22','C23','C24','C25','C26','C27', 'C28','C29', 'C210','C211', 'C212','C213', 'C214','C215', 'C216'],\
            ['C32','C33','C34','C35','C36','C37', 'C38','C39', 'C310','C311', 'C312','C313', 'C314','C315', 'C316']]
tailhydrdppc=[['H2R', 'H2S','H3R', 'H3S', 'H4R', 'H4S', 'H5R', 'H5S', 'H6R', 'H6S','H7R','H7S','H8R', 'H8S', 'H9R', 'H9S', 'H10R', 'H10S','H11R', 'H11S',\
   'H12R', 'H12S','H13R', 'H13S','H14R', 'H14S', 'H15R', 'H15S','H16R', 'H16S', 'H16T'],\
          ['H2X', 'H2Y','H3X', 'H3Y','H4X', 'H4Y','H5X', 'H5Y', 'H6X', 'H6Y','H7X', 'H7Y','H8X', 'H8Y','H9X', 'H9Y','H10X', 'H10Y','H11X', 'H11Y',\
  'H12X', 'H12Y', 'H13X', 'H13Y','H14X', 'H14Y','H15X', 'H15Y','H16X', 'H16Y','H16Z']]
tailcarbonsdupc=[tailcarbonsdppc[0]+['C217','C218'],tailcarbonsdppc[1]+['C317','C318']]
tailhydrdupc=[tailhydrdppc[0].copy()+['H17X','H17Y','H18X','H18Y','H18Z'],tailhydrdppc[1].copy()+['H17X','H17Y','H18X','H18Y','H18Z']]
tailhydrdupc[0].remove('H16T')
tailhydrdupc[1].remove('H16Z')

scd_tail_atoms_of={\
    'DPPC':[['C22','C24','C26', 'C28', 'C210', 'C212', 'C214', 'C216'],\
        ['C32', 'C34', 'C36','C38', 'C310', 'C312', 'C314', 'C316']],\
    'CHL1':[['C3','C17']],\
    'DUPC':[['C21','C23','C25', 'C27','C29','C211', 'C213', 'C215','C217'],\
        ['C31','C33','C35', 'C37','C39','C311', 'C313', 'C315','C317']]}
head_atoms_of={\
    'DPPC':['C11','C12','C13','C14','C15','N','P','O11','O12','O13','O14','C1','C2','O21','C21','O22','C3','O31','O32','C31']+headhydr,\
    'DUPC':['C11','C12','C13','C14','C15','N','P','O11','O12','O13','O14','C1','C2','O21','C21','O22','C3','O31','O32','C31']+headhydr,\
    'CHL1':['all']}
tail_atoms_of={\
    'DPPC':[tailcarbonsdppc[0],tailhydrdppc[0],tailcarbonsdppc[1],tailhydrdppc[1]],\
    'DUPC':[tailcarbonsdupc[0],tailhydrdupc[0],tailcarbonsdupc[1],tailhydrdupc[1]],\
    'CHL1':['all']}
central_atom_of={\
    'DPPC':'P',\
    'CHL1':'O3',\
    'DUPC':'P',\
    'DOPC':'P'}